#!/usr/bin/env python3
"""Bacteria Atlas: paste bacterial species names into a text file, get a static website.

    python build_site.py species.txt                 # writes ./site/
    python build_site.py species.txt --out docs --title "Bacteria Atlas" --email me@lab.org

Input: one species per line (blank lines and #comments are ignored). Messy lines such as
"SKIPPED Cutibacterium acnes: no output yet" or "Escherichia_coli" are cleaned automatically.

For every species the script fetches, live from the internet (standard library only):
  * NCBI Taxonomy      -> taxid, lineage, former names, genetic code
  * NCBI Datasets      -> reference genome size, GC%, assembly level and accession, assembly count
  * NCBI PubMed        -> number of papers mentioning the species
  * Wikipedia          -> species article (or genus article) summary + link + lead picture
  * Wikimedia Commons  -> a picture of the species when its Wikipedia article has none
and builds one HTML page per species plus an index page with a search box and a genome-size chart.
Pictures are downloaded into site/images/ with credit and a link to the file page.
Everything fetched is cached in ./.bacteria_cache so re-runs are fast (use --refresh to refetch).
If the network is down or a species is unknown, its page is still built with whatever was found.
"""
import argparse
import colorsys
import datetime
import html
import json
import pathlib
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

HERE = pathlib.Path(__file__).parent
CACHE = HERE / ".bacteria_cache"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
DATASETS = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/taxon/"
WIKI_SUMMARY = "https://en.wikipedia.org/api/rest_v1/page/summary/"
COMMONS_API = "https://commons.wikimedia.org/w/api.php"
BACTERIA_WORDS = re.compile(r"bacteri|archae|microb|prokaryot|cyanobacter|proteobacter|actinomycet|firmicute", re.I)
RANKS = ["domain", "kingdom", "phylum", "class", "order", "family", "genus"]

EMAIL = ""
_last_ncbi = [0.0]


# ---------------------------------------------------------------- networking

def user_agent():
    return f"BacteriaAtlas/1.0 (static site generator; {EMAIL or 'no contact given'})"


PLACEHOLDERS = {"missing", "not applicable", "not collected", "not available", "not provided", "unknown", "na", "n/a",
                "none", "restricted access", "not determined", "not recorded", "other", "unspecified", "not known"}
NET_ERRORS = [0]  # failed lookups (not 404s); a species with any is not cached


def http_get(url, retries=3, timeout=30):
    """Return response bytes; None on 404 (genuinely absent) or after repeated failures."""
    for attempt in range(retries + 1):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": user_agent()})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.read()
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return None
            # 429 = rate limited: wait as long as the server asks, then try again
            try:
                wait = float(e.headers.get("Retry-After", ""))
            except (TypeError, ValueError):
                wait = 4.0 * (attempt + 1) if e.code in (429, 503) else 1.5 * (attempt + 1)
            time.sleep(min(wait, 30))
        except (urllib.error.URLError, TimeoutError, OSError):
            time.sleep(1.5 * (attempt + 1))
    NET_ERRORS[0] += 1
    return None


def get_json(url):
    data = http_get(url)
    if not data:
        return None
    try:
        return json.loads(data.decode("utf-8"))
    except ValueError:
        return None


def ncbi_url(path, **params):
    """NCBI asks for at most ~3 requests per second without an API key."""
    wait = 0.4 - (time.time() - _last_ncbi[0])
    if wait > 0:
        time.sleep(wait)
    _last_ncbi[0] = time.time()
    params["tool"] = "bacteria_atlas"
    if EMAIL:
        params["email"] = EMAIL
    return EUTILS + path + "?" + urllib.parse.urlencode(params)


# ------------------------------------------------------------- input parsing

def parse_names(text):
    """Pull 'Genus species' pairs out of free-form text, keeping order, dropping duplicates."""
    names, seen = [], set()
    for chunk in re.split(r"[\n;,]", text):
        chunk = chunk.strip().replace("_", " ")
        if not chunk or chunk.startswith("#"):
            continue
        chunk = re.sub(r"^(SKIPPED|DONE|FAILED)\s+", "", chunk, flags=re.I).split(":")[0].strip()
        m = re.match(r"^\[?([A-Z][a-z]+)\]?\s+([a-z][a-z.\-]+)", chunk)
        if m:
            name = f"{m.group(1)} {m.group(2)}"
            if name.lower() not in seen:
                seen.add(name.lower())
                names.append(name)
    return names


# --------------------------------------------------------------------- NCBI

def ncbi_taxonomy(name):
    out = {"taxid": None, "ncbi_name": None, "rank": None, "lineage": {}, "synonyms": [], "genetic_code": None}
    for term in (f'"{name}"[Scientific Name]', f'"{name}"[All Names]'):
        d = get_json(ncbi_url("esearch.fcgi", db="taxonomy", term=term, retmode="json"))
        ids = ((d or {}).get("esearchresult") or {}).get("idlist") or []
        if ids:
            out["taxid"] = int(ids[0])
            break
    if not out["taxid"]:
        return out
    raw = http_get(ncbi_url("efetch.fcgi", db="taxonomy", id=out["taxid"], retmode="xml"))
    if not raw:
        return out
    try:
        tax = ET.fromstring(raw).find("Taxon")
    except ET.ParseError:
        return out
    if tax is None:
        return out
    out["ncbi_name"] = tax.findtext("ScientificName")
    out["rank"] = tax.findtext("Rank")
    out["genetic_code"] = tax.findtext("GeneticCode/GCName")
    for t in tax.findall("LineageEx/Taxon"):
        rank = t.findtext("Rank")
        if rank == "superkingdom":
            rank = "domain"
        if rank in RANKS:
            out["lineage"][rank] = t.findtext("ScientificName")
    for tag in ("Synonym", "GenbankSynonym"):
        for s in tax.findall(f"OtherNames/{tag}"):
            s = (s.text or "").strip()
            if s and not re.search(r"\d", s) and s != out["ncbi_name"] and s not in out["synonyms"]:
                out["synonyms"].append(s)
    out["synonyms"] = out["synonyms"][:6]
    return out


def ncbi_genome(taxid):
    """Reference (else RefSeq, else any) assembly stats plus the total assembly count."""
    out = {}
    base = f"{DATASETS}{taxid}/dataset_report?page_size=1"
    for flt in ("&filters.reference_only=true", "&filters.assembly_source=refseq", ""):
        d = get_json(base + flt)
        reports = (d or {}).get("reports") or []
        if reports:
            r = reports[0]
            st = r.get("assembly_stats") or {}
            info = r.get("assembly_info") or {}
            try:
                out["size_mb"] = round(int(st["total_sequence_length"]) / 1e6, 2)
            except (KeyError, ValueError, TypeError):
                pass
            if st.get("gc_percent") is not None:
                out["gc"] = st["gc_percent"]
            out["level"] = info.get("assembly_level")
            out["accession"] = r.get("accession")
            out["strain"] = ((r.get("organism") or {}).get("infraspecific_names") or {}).get("strain")
            out["reference"] = flt.endswith("true")
            break
    # Sample records of up to 40 assemblies tell where strains of this species are usually collected.
    d = get_json(f"{DATASETS}{taxid}/dataset_report?page_size=40")
    if d and d.get("total_count") is not None:
        out["n_assemblies"] = d["total_count"]
    found = {"isolation_source": {}, "host": {}}
    for r in (d or {}).get("reports") or []:
        attrs = {a.get("name"): a.get("value") for a in (((r.get("assembly_info") or {}).get("biosample") or {}).get("attributes") or [])}
        for key, counts in found.items():
            v = (attrs.get(key) or "").strip().lower()
            if v and v not in PLACEHOLDERS and not re.fullmatch(r"[\W\d_]*", v):
                counts[v] = counts.get(v, 0) + 1
    top = []
    for key, n in (("isolation_source", 3), ("host", 1)):
        top += sorted(found[key], key=lambda v: -found[key][v])[:n]
    if top:
        out["found_in"] = ", ".join(dict.fromkeys(top))
    merged = {}
    for counts in found.values():
        for v, n in counts.items():
            merged[v] = merged.get(v, 0) + n
    out["sources"] = dict(sorted(merged.items(), key=lambda kv: -kv[1])[:20])  # raw evidence for environment()
    return out


def pubmed_text(name, retmax=30):
    """Titles + abstracts of the most relevant papers: the habitat evidence of last resort."""
    d = get_json(ncbi_url("esearch.fcgi", db="pubmed", term=f'"{name}"[Title/Abstract]', retmax=retmax,
                          sort="relevance", retmode="json"))
    ids = ((d or {}).get("esearchresult") or {}).get("idlist") or []
    if not ids:
        return ""
    raw = http_get(ncbi_url("efetch.fcgi", db="pubmed", id=",".join(ids), retmode="xml"))
    if not raw:
        return ""
    try:
        root = ET.fromstring(raw)
    except ET.ParseError:
        return ""
    parts = []
    for art in root.iter("PubmedArticle"):
        for tag in ("ArticleTitle", "AbstractText"):
            for el in art.iter(tag):
                parts.append("".join(el.itertext()))
    return " ".join(parts)[:40000]


def pubmed_count(name):
    d = get_json(ncbi_url("esearch.fcgi", db="pubmed", term=f'"{name}"[Title/Abstract]', retmode="json", rettype="count"))
    try:
        return int(d["esearchresult"]["count"])
    except (TypeError, KeyError, ValueError):
        return None


# ---------------------------------------------------------- Wikipedia / Commons

def wiki_search(query, must_start):
    """Wikipedia full-text search; return the first title starting with `must_start`, else None."""
    params = dict(action="query", list="search", srsearch=query, srlimit=5, format="json")
    d = get_json("https://en.wikipedia.org/w/api.php?" + urllib.parse.urlencode(params))
    for hit in ((d or {}).get("query") or {}).get("search") or []:
        if hit.get("title", "").lower().startswith(must_start.lower()):
            return hit["title"]
    return None


def wiki_fulltext(title, limit=12000):
    """Plain text of the whole article (habitat/ecology sections live below the lead paragraph)."""
    time.sleep(0.25)
    params = dict(action="query", prop="extracts", explaintext=1, exsectionformat="plain", redirects=1,
                  titles=title, format="json")
    d = get_json("https://en.wikipedia.org/w/api.php?" + urllib.parse.urlencode(params))
    pages = ((d or {}).get("query") or {}).get("pages") or {}
    text = next(iter(pages.values()), {}).get("extract", "") if pages else ""
    return re.split(r"\n(?:References|External links|See also|Further reading|Notes)\n", text)[0][:limit]


def wiki_summary(title):
    time.sleep(0.25)  # be polite to Wikipedia
    d = get_json(WIKI_SUMMARY + urllib.parse.quote(title.replace(" ", "_"), safe=""))
    if not d or d.get("type") == "disambiguation" or not d.get("extract"):
        return None
    if not BACTERIA_WORDS.search(d["extract"] + " " + (d.get("description") or "")):
        return None  # a page with that title, but not about a microbe
    thumb = (d.get("thumbnail") or {}).get("source")
    return {
        "title": d.get("title"),
        "extract": d["extract"],
        "url": (d.get("content_urls") or {}).get("desktop", {}).get("page")
               or "https://en.wikipedia.org/wiki/" + urllib.parse.quote(title.replace(" ", "_")),
        "thumb": thumb,
    }


def find_wikipedia(name, synonyms, genus):
    """Return (kind, summary): kind is 'species', 'genus' or None."""
    for t in [name] + synonyms:
        s = wiki_summary(t)
        if s:
            title_words = (s["title"] or "").lower().split()
            if len(title_words) >= 2:
                return "species", s
            if title_words and title_words[0] == genus.lower():
                return "genus", s
    s = wiki_summary(genus)
    if not s:  # last resort: let Wikipedia's own search find the species, then the genus
        for q, start in ((name, genus), (f"{genus} bacteria genus", genus)):
            title = wiki_search(q, start)
            s = wiki_summary(title) if title else None
            if s:
                words = (s["title"] or "").lower().split()
                return ("species" if len(words) >= 2 else "genus"), s
    return ("genus", s) if s else (None, None)


def strip_tags(s):
    return re.sub(r"\s+", " ", re.sub(r"<[^>]+>", "", s or "")).strip()


def commons_image(names):
    """Search Wikimedia Commons for a file whose title contains the exact species name."""
    for q in names:
        params = dict(action="query", generator="search", gsrnamespace=6, gsrsearch=f'"{q}"', gsrlimit=10,
                      prop="imageinfo", iiprop="url|extmetadata|mime", iiurlwidth=640, format="json")
        d = get_json(COMMONS_API + "?" + urllib.parse.urlencode(params))
        pages = ((d or {}).get("query") or {}).get("pages") or {}
        for p in sorted(pages.values(), key=lambda p: p.get("index", 0)):
            ii = (p.get("imageinfo") or [{}])[0]
            if ii.get("mime") not in ("image/jpeg", "image/png", "image/gif", "image/svg+xml"):
                continue
            if q.lower() not in p.get("title", "").lower().replace("_", " "):
                continue
            meta = ii.get("extmetadata") or {}
            return {
                "url": ii.get("thumburl") or ii.get("url"),
                "page": ii.get("descriptionurl"),
                "credit": strip_tags((meta.get("Artist") or {}).get("value")),
                "license": strip_tags((meta.get("LicenseShortName") or {}).get("value")),
            }
    return None


def file_page_from_thumb(url):
    """Best-effort link to the Commons/Wikipedia file page behind an upload.wikimedia.org URL."""
    if not url:
        return None
    parts = urllib.parse.urlparse(url).path.split("/")
    if "thumb" in parts:
        name = parts[parts.index("thumb") + 3] if len(parts) > parts.index("thumb") + 3 else None
    else:
        name = parts[-1]
    if not name:
        return None
    host = "en.wikipedia.org" if "/wikipedia/en/" in url else "commons.wikimedia.org"
    return f"https://{host}/wiki/File:{name}"


def choose_image(name, synonyms, kind, wiki, genus_wiki):
    if kind == "species" and wiki and wiki.get("thumb"):
        return {"url": wiki["thumb"], "page": file_page_from_thumb(wiki["thumb"]), "credit": "",
                "license": "see file page", "source": "Wikipedia article", "scope": "species"}
    c = commons_image([name] + synonyms[:2])
    if c:
        c.update(source="Wikimedia Commons", scope="species")
        return c
    g = wiki if kind == "genus" else genus_wiki
    if g and g.get("thumb"):
        return {"url": g["thumb"], "page": file_page_from_thumb(g["thumb"]), "credit": "",
                "license": "see file page", "source": f"Wikipedia article '{g['title']}'", "scope": "genus"}
    return None


# ------------------------------------------------------------------- fetching

def slug(name):
    return re.sub(r"[^a-z0-9]+", "_", name.lower()).strip("_")


def fetch_species(name, refresh=False):
    cache = CACHE / f"{slug(name)}.json"
    if cache.exists() and not refresh:
        try:
            return json.loads(cache.read_text(encoding="utf-8"))
        except ValueError:
            pass
    genus = name.split()[0]
    errors_before = NET_ERRORS[0]
    tax = ncbi_taxonomy(name)
    genome = ncbi_genome(tax["taxid"]) if tax["taxid"] else {}
    syn = tax["synonyms"]
    kind, wiki = find_wikipedia(name, syn, genus)
    genus_wiki = None
    if kind == "species":
        genus_wiki = wiki_summary(genus)
    d = {
        "name": name, "genus": genus, "tax": tax, "genome": genome, "pubmed": pubmed_count(name),
        "wiki_kind": kind, "wiki": wiki,
        "wiki_full": wiki_fulltext(wiki["title"]) if wiki else "",
        "genus_full": wiki_fulltext(genus, 6000) if kind == "species" else "",
        "pubmed_text": pubmed_text(name),
        "image": choose_image(name, syn, kind, wiki, genus_wiki),
        "fetched": datetime.date.today().isoformat(),
    }
    if NET_ERRORS[0] == errors_before:  # never cache a result that had failed lookups
        CACHE.mkdir(exist_ok=True)
        cache.write_text(json.dumps(d, indent=1), encoding="utf-8")
    else:
        print("    warning: some lookups failed (network/rate limit); not cached, re-run to retry", flush=True)
    return d


def download_image(d, out_dir):
    """Download the chosen picture to out_dir/images; return its relative path or None."""
    img = d.get("image")
    if not img or not img.get("url"):
        return None
    ext = pathlib.Path(urllib.parse.urlparse(img["url"]).path).suffix.lower()
    if ext not in (".jpg", ".jpeg", ".png", ".gif", ".webp"):
        ext = ".jpg" if ext in ("", ".svg") and "svg" not in img["url"].lower() else ".png"
    rel = f"images/{slug(d['name'])}{ext}"
    dest = out_dir / rel
    if dest.exists() and dest.stat().st_size > 0:
        return rel
    # Wikipedia thumbnails are small (~330px); ask for a larger rendition first. Wikimedia only
    # serves standard thumbnail widths (e.g. 500, 960), so 640 would be rejected.
    candidates = [re.sub(r"/\d+px-", "/960px-", img["url"]), img["url"]]
    for u in dict.fromkeys(candidates):
        data = http_get(u)
        if data:
            dest.parent.mkdir(parents=True, exist_ok=True)
            dest.write_bytes(data)
            return rel
    return None


# ------------------------------------------------------------------- HTML

CSS = """
:root{--bg:#f7f8f5;--fg:#1d2420;--card:#fff;--muted:#5d6a62;--accent:#2f7d5a;--line:#dde3dd}
@media (prefers-color-scheme:dark){:root{--bg:#14181a;--fg:#e6ece8;--card:#1d2326;--muted:#9aa8a0;--accent:#6cc29b;--line:#2c3539}}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--fg);font:16px/1.55 system-ui,sans-serif}
main{max-width:900px;margin:0 auto;padding:24px 16px 64px}
h1{margin:.2em 0 .1em}h1 i,.card i{font-style:italic}
a{color:var(--accent)}
.sub{color:var(--muted);margin-bottom:16px}
.grid{display:grid;grid-template-columns:repeat(auto-fill,minmax(230px,1fr));gap:12px}
.card{background:var(--card);border:1px solid var(--line);border-radius:10px;overflow:hidden;text-decoration:none;color:inherit;display:block}
.card:hover{border-color:var(--accent)}
.card .pic{aspect-ratio:4/3;background:var(--line);overflow:hidden}
.card .pic img,.card .pic svg{width:100%;height:100%;object-fit:cover;display:block}
.card .txt{padding:10px 12px}
.card b{display:block;margin-bottom:4px}
.tag.env{background:var(--accent);color:var(--bg)}
.tag{display:inline-block;font-size:.78rem;padding:1px 8px;border-radius:99px;background:var(--line);color:var(--muted);margin:0 4px 4px 0}
input[type=search]{width:100%;padding:10px 12px;margin:0 0 14px;border:1px solid var(--line);border-radius:8px;background:var(--card);color:var(--fg);font:inherit}
table{width:100%;border-collapse:collapse;background:var(--card);border:1px solid var(--line);border-radius:10px;margin:8px 0 16px}
td,th{text-align:left;padding:8px 12px;border-bottom:1px solid var(--line);vertical-align:top}
th{width:34%;color:var(--muted);font-weight:500}
.note{border-left:4px solid var(--accent);background:var(--card);padding:10px 14px;margin:18px 0;border-radius:4px}
nav{margin-bottom:8px;font-size:.9rem;display:flex;justify-content:space-between;gap:8px}
.hero{margin:0 0 14px 20px;float:right;width:min(42%,360px)}
.hero img,.hero svg{width:100%;max-height:420px;object-fit:contain;display:block;border:1px solid var(--line);border-radius:10px;background:var(--card)}
.clear{clear:both}
@media (max-width:640px){.hero{float:none;width:100%;margin:0 0 14px}}
.hero figcaption{font-size:.82rem;color:var(--muted);margin-top:4px}
.links{margin:0 0 16px}
.links a,.links span{display:inline-block;margin:0 8px 8px 0;padding:5px 14px;border:1px solid var(--line);border-radius:99px;text-decoration:none;background:var(--card)}
.links a.wiki{border-color:var(--accent);font-weight:600}
.links span{border-style:dashed;color:var(--muted)}
.bar{display:grid;grid-template-columns:200px 1fr 70px;gap:8px;align-items:center;font-size:.9rem;margin:3px 0}
.bar i{font-style:italic;overflow:hidden;text-overflow:ellipsis;white-space:nowrap}
.bar div{height:12px;border-radius:6px;background:var(--accent)}
.small{font-size:.82rem;color:var(--muted)}
"""


def e(v):
    return html.escape(str(v))


def hue_for(name):
    return (sum(ord(c) * (i + 3) for i, c in enumerate(name)) % 360) / 360.0


def placeholder_svg(name):
    """Generated stand-in when no real picture exists: rods coloured from the name."""
    r, g, b = colorsys.hls_to_rgb(hue_for(name), 0.55, 0.45)
    c = f"#{int(r * 255):02x}{int(g * 255):02x}{int(b * 255):02x}"
    rods = "".join(
        f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{c}" stroke-width="26" stroke-linecap="round"/>'
        f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="#fff" stroke-opacity=".28" stroke-width="12" stroke-linecap="round"/>'
        for x1, y1, x2, y2 in [(90, 70, 190, 90), (250, 150, 340, 120), (120, 190, 220, 175), (400, 60, 470, 95), (430, 190, 520, 205)])
    return (f'<svg viewBox="0 0 600 260" role="img" aria-label="Placeholder for {e(name)}">'
            f'<rect width="600" height="260" fill="#e9eee9"/>{rods}'
            f'<text x="12" y="248" font-size="13" fill="#5d6a62" font-family="system-ui,sans-serif">'
            f'no picture found online: generic schematic</text></svg>')


def page(title, body, script=""):
    return (f'<!doctype html><html lang="en"><head><meta charset="utf-8">'
            f'<meta name="viewport" content="width=device-width,initial-scale=1">'
            f'<title>{e(title)}</title><style>{CSS}</style></head><body><main>{body}</main>{script}</body></html>')


def wiki_chip(d):
    w = d.get("wiki")
    if d["wiki_kind"] == "species" and w:
        return f'<a class="wiki" href="{e(w["url"])}">Wikipedia article</a>'
    if d["wiki_kind"] == "genus" and w:
        return f'<a class="wiki" href="{e(w["url"])}">Wikipedia (genus: {e(w["title"])})</a>'
    return '<span>No Wikipedia article found</span>'


def link_chips(d):
    t, q = d["tax"], urllib.parse.quote_plus(d["name"])
    chips = [wiki_chip(d)]
    if t.get("taxid"):
        chips.append(f'<a href="https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id={t["taxid"]}">NCBI Taxonomy</a>')
        chips.append(f'<a href="https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon={t["taxid"]}">NCBI genomes</a>')
    else:
        chips.append(f'<a href="https://www.ncbi.nlm.nih.gov/taxonomy/?term={q}">Search NCBI Taxonomy</a>')
    chips.append(f'<a href="https://pubmed.ncbi.nlm.nih.gov/?term={q}">PubMed</a>')
    chips.append(f'<a href="https://lpsn.dsmz.de/search?word={q}">LPSN</a>')
    chips.append(f'<a href="https://bacdive.dsmz.de/search?search={q}">BacDive</a>')
    chips.append(f'<a href="https://scholar.google.com/scholar?q={q}">Scholar</a>')
    return "".join(chips)


def hero_html(d, rel):
    img = d.get("image")
    if rel and img:
        scope = ("Representative image from the genus article, not necessarily this species. "
                 if img["scope"] == "genus" else "")
        bits = [scope + f'Source: {e(img["source"])}']
        if img.get("credit"):
            bits.append(f'credit: {e(img["credit"])}')
        if img.get("license"):
            bits.append(f'licence: {e(img["license"])}')
        if img.get("page"):
            bits.append(f'<a href="{e(img["page"])}">file page</a>')
        return (f'<figure class="hero"><img src="{e(rel)}" alt="{e(d["name"])}">'
                f'<figcaption>{" | ".join(bits)}</figcaption></figure>')
    return f'<figure class="hero">{placeholder_svg(d["name"])}</figure>'


GRAM_POSITIVE_PHYLA = {"Actinomycetota", "Bacillota", "Firmicutes", "Actinobacteria"}
GRAM_NEGATIVE_PHYLA = {
    "Pseudomonadota", "Proteobacteria", "Bacteroidota", "Cyanobacteriota", "Spirochaetota", "Chlamydiota",
    "Campylobacterota", "Fusobacteriota", "Planctomycetota", "Verrucomicrobiota", "Acidobacteriota",
    "Thermotogota", "Aquificota", "Myxococcota", "Desulfobacterota", "Nitrospirota", "Synergistota",
    "Deinococcota", "Gemmatimonadota", "Rhodothermota", "Chlorobiota", "Fibrobacterota", "Elusimicrobiota",
}


def gram_stain(d):
    """'Gram-positive' / 'Gram-negative' from the Wikipedia text if it says so, else inferred from the phylum."""
    w = d.get("wiki") or {}
    text = (w.get("extract") or "").lower()
    pos, neg = re.search(r"gram[- ]positive", text), re.search(r"gram[- ]negative", text)
    if pos and not neg:
        return "Gram-positive"
    if neg and not pos:
        return "Gram-negative"
    lin = d["tax"]["lineage"]
    if lin.get("domain") == "Archaea":
        return "Archaea"
    if lin.get("class") == "Negativicutes":
        return "Gram-negative"
    if lin.get("phylum") == "Mycoplasmatota":
        return "No cell wall"
    if lin.get("phylum") in GRAM_POSITIVE_PHYLA:
        return "Gram-positive"
    if lin.get("phylum") in GRAM_NEGATIVE_PHYLA:
        return "Gram-negative"
    return None


# Habitat keywords per category. Skin/gut/mouth/clinical vocabulary all counts as "Human".
ENVIRONMENTS = [
    ("Human", r"\b(?:humans?|patients?|clinical|sputum|skin|acne|dermat\w*|sebaceous|wounds?|blood\w*|urin\w*|gut|intestin\w*|oral|mouth|nasal|infections?|infect\w*|pathogen\w*|disease\w*|hospital\w*|homo sapiens|commensals?|opportunistic)\b"),
    ("Animal", r"\b(?:animals?|dogs?|cats?|cattle|cows?|pigs?|swine|poultry|chickens?|livestock|mammals?|fish|insects?|ticks?|bees?|rodents?|mice|mouse|bovine|ruminants?|sheep|horses?)\b"),
    ("Soil", r"\b(?:soils?|rhizosphere|compost|mud|desert|forest|peat|agricultur\w*)\b"),
    ("Plant", r"\b(?:plants?|roots?|nodul\w*|leaf|leaves|legumes?|seeds?|tobacco|rice|wheat|phyllosphere|endophyt\w*|crops?)\b"),
    ("Freshwater", r"\b(?:freshwater|rivers?|lakes?|ponds?|streams?|reservoirs?|groundwater|drinking water|water)\b"),
    ("Marine", r"\b(?:marine|sea|seawater|oceans?|coastal|estuar\w*|tidal|saltwater|sediments?|reefs?)\b"),
    ("Extreme (hot, cold, salty)", r"\b(?:hot springs?|thermal|hydrothermal|thermophil\w*|geothermal|hypersaline|saline|salt\w*|glaciers?|antarctic\w*|arctic|polar|cryoconite|ice|permafrost|acidic|alkaline)\b"),
    ("Airborne", r"\b(?:air|airborne|aerosols?|atmospher\w*|dust)\b"),
    ("Food & industry", r"\b(?:food\w*|dairy|milk|cheese|ferment\w*|wine|beer|yogurt|probiotics?|industrial|biofilms?|wastewater|sewage|sludge|paper)\b"),
]
# how much each evidence source counts; each source is first turned into shares (summing to 1)
# so a long PubMed text cannot drown out a short but precise NCBI isolation-source list
ENV_WEIGHTS = {"NCBI strains": 3.0, "Wikipedia lead": 2.0, "Wikipedia article": 1.5, "PubMed": 1.5, "Genus article": 0.5}


def _shares(text):
    hits = {label: len(re.findall(p, (text or "").lower())) for label, p in ENVIRONMENTS}
    total = sum(hits.values())
    return {k: v / total for k, v in hits.items() if v} if total else {}


def environment_detail(d):
    """(label, evidence) for the main habitat, weighing NCBI strain sources, Wikipedia and PubMed."""
    w = d.get("wiki") or {}
    strains = " ; ".join(f"{v} ;" * n for v, n in (d["genome"].get("sources") or {}).items())
    texts = {
        "NCBI strains": strains or (d["genome"].get("found_in") or ""),
        "Wikipedia lead": (w.get("extract") or ""),
        "Wikipedia article": d.get("wiki_full"),
        "PubMed": d.get("pubmed_text"),
        "Genus article": d.get("genus_full"),
    }
    dists = {src: _shares(t) for src, t in texts.items()}
    dists = {s: dist for s, dist in dists.items() if dist}
    if not dists:
        return None, ""
    used = sum(ENV_WEIGHTS[s] for s in dists)
    score = {c: sum(ENV_WEIGHTS[s] * dist.get(c, 0) for s, dist in dists.items()) / used for c, _ in ENVIRONMENTS}
    best = max(score, key=score.get)
    if score[best] < 0.2:
        return None, ""
    evidence = " &middot; ".join(f"{s} {dist.get(best, 0):.0%}" for s, dist in dists.items())
    return best, evidence


def environment(d):
    return environment_detail(d)[0]


def facts_rows(d):
    t, g = d["tax"], d["genome"]
    rows = []
    if gram_stain(d):
        rows.append(("Gram stain", e(gram_stain(d))))
    env, evidence = environment_detail(d)
    if env:
        rows.append(("Typical environment", e(env) + f' <span class="small">(best guess; share of evidence pointing to it: {evidence})</span>'))
    if t.get("taxid"):
        rows.append(("NCBI taxid", f'{t["taxid"]} ({e(t.get("rank") or "")})'))
    for r in RANKS:
        if t["lineage"].get(r):
            rows.append((r.capitalize(), e(t["lineage"][r])))
    if t.get("synonyms"):
        rows.append(("Also known as", ", ".join(f"<i>{e(s)}</i>" for s in t["synonyms"])))
    if g.get("found_in"):
        rows.append(("Found in", e(g["found_in"]) + ' <span class="small">(most common collection sources among sequenced strains in NCBI)</span>'))
    if g.get("size_mb"):
        kind ="reference genome" if g.get("reference") else "representative assembly"
        rows.append(("Genome size", f'{g["size_mb"]} Mb ({kind}, {e(g.get("level") or "?")}'
                                    f'{", strain " + e(g["strain"]) if g.get("strain") else ""})'))
    if g.get("gc") is not None:
        rows.append(("GC content", f'{g["gc"]}%'))
    if g.get("accession"):
        rows.append(("Assembly", f'<a href="https://www.ncbi.nlm.nih.gov/datasets/genome/{e(g["accession"])}/">{e(g["accession"])}</a>'))
    if g.get("n_assemblies") is not None:
        rows.append(("Genome assemblies in NCBI", e(g["n_assemblies"])))
    if t.get("genetic_code"):
        rows.append(("Genetic code", e(t["genetic_code"])))
    if d.get("pubmed") is not None:
        q = urllib.parse.quote_plus(f'"{d["name"]}"[Title/Abstract]')
        rows.append(("PubMed papers", f'<a href="https://pubmed.ncbi.nlm.nih.gov/?term={q}">{d["pubmed"]:,}</a> mention it in title/abstract'))
    return rows


def species_page(d, rel, prev, nxt, site_title):
    t = d["tax"]
    lin = " &middot; ".join(e(t["lineage"][r]) for r in ("phylum", "class", "family") if t["lineage"].get(r))
    trs = "".join(f"<tr><th>{k}</th><td>{v}</td></tr>" for k, v in facts_rows(d))
    w = d.get("wiki")
    if w:
        label = "article" if d["wiki_kind"] == "species" else f'genus article "{w["title"]}"'
        about = (f'<h2>About</h2><p>{e(w["extract"])}</p>'
                 f'<p class="small">Text from the Wikipedia {e(label)} (CC BY-SA). '
                 f'<a href="{e(w["url"])}">Read the full page &rarr;</a></p>')
    else:
        about = '<h2>About</h2><p class="small">No Wikipedia article was found for this species or its genus.</p>'
    if not t.get("taxid"):
        about += ('<div class="note">NCBI Taxonomy has no entry under this exact name. '
                  'Check the spelling or search for a synonym.</div>')
    navl = f'<a href="{slug(prev)}.html">&larr; <i>{e(prev)}</i></a>' if prev else "<span></span>"
    navr = f'<a href="{slug(nxt)}.html"><i>{e(nxt)}</i> &rarr;</a>' if nxt else "<span></span>"
    body = (f'<nav><a href="index.html">&larr; {e(site_title)}</a></nav>'
            f'{hero_html(d, rel)}'
            f'<h1><i>{e(d["name"])}</i></h1><div class="sub">{lin or "classification unavailable"}</div>'
            f'<div class="links">{link_chips(d)}</div>'
            f'{about}<div class="clear"></div><h2>Key facts</h2><table>{trs}</table>'
            f'<nav>{navl}{navr}</nav>'
            f'<p class="small">Fetched {e(d["fetched"])} from NCBI, Wikipedia and Wikimedia Commons. Verify before citing.</p>')
    return page(f'{d["name"]} | {site_title}', body)


def genome_chart(items):
    rows = sorted([(d["name"], d["genome"]["size_mb"]) for d in items if d["genome"].get("size_mb")], key=lambda r: -r[1])
    if not rows:
        return ""
    top = rows[0][1]
    bars = "".join(
        f'<div class="bar"><i>{e(n)}</i><div style="width:{max(2, 100 * s / top):.1f}%"></div><span>{s} Mb</span></div>'
        for n, s in rows)
    return f'<h2>Genome sizes</h2>{bars}<p class="small">NCBI reference (else representative) assembly.</p>'


def index_page(items, rels, site_title):
    cards = ""
    for d in items:
        rel = rels[d["name"]]
        pic = f'<img src="{e(rel)}" alt="" loading="lazy">' if rel else placeholder_svg(d["name"])
        lin = d["tax"]["lineage"]
        labels = [("Genus", d["genus"]), ("Family", lin.get("family")), (None, gram_stain(d)), (None, environment(d))]
        tags = "".join(f'<span class="tag{" env" if i == 3 else ""}">{e(v)}</span>'
                       for i, (_, v) in enumerate(labels) if v)
        search = " ".join(str(v) for _, v in labels if v) + " " + d["name"]
        cards += (f'<a class="card" data-q="{e(search.lower())}" href="{slug(d["name"])}.html">'
                  f'<div class="pic">{pic}</div><div class="txt"><b><i>{e(d["name"])}</i></b>{tags}</div></a>')
    script = ('<script>document.getElementById("q").addEventListener("input",function(){var v=this.value.toLowerCase();'
              'document.querySelectorAll(".card").forEach(function(c){c.style.display=c.dataset.q.indexOf(v)<0?"none":""})})</script>')
    body = (f'<h1>{e(site_title)}</h1><div class="sub">{len(items)} species. Pictures, summaries and genome statistics '
            f'fetched on {datetime.date.today().isoformat()} from NCBI, Wikipedia and Wikimedia Commons.</div>'
            f'<input type="search" id="q" placeholder="Filter species...">'
            f'<div class="grid">{cards}</div>{genome_chart(items)}'
            f'<p class="small">Pictures keep their original licences; see the credit under each picture. '
            f'Wikipedia text is CC BY-SA. Verify facts against primary sources before citing.</p>')
    return page(site_title, body, script)


# ------------------------------------------------------------------- main

def main():
    global EMAIL
    ap = argparse.ArgumentParser(description="Build a static website from a list of bacterial species.")
    ap.add_argument("names_file", help="text file with one species per line")
    ap.add_argument("--out", default="site", help="output folder (default: site)")
    ap.add_argument("--title", default="Bacteria Atlas")
    ap.add_argument("--email", default="", help="contact email, sent to NCBI/Wikimedia as they request")
    ap.add_argument("--refresh", action="store_true", help="ignore the cache and refetch everything")
    ap.add_argument("--no-images", action="store_true", help="skip downloading pictures")
    args = ap.parse_args()
    EMAIL = args.email

    names = parse_names(pathlib.Path(args.names_file).read_text(encoding="utf-8", errors="replace"))
    if not names:
        sys.exit("No 'Genus species' names found in the input file.")
    out = pathlib.Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    items, rels = [], {}
    for i, name in enumerate(names, 1):
        print(f"[{i}/{len(names)}] {name}", flush=True)
        d = fetch_species(name, args.refresh)
        rels[name] = None if args.no_images else download_image(d, out)
        print(f"    taxid={d['tax'].get('taxid')} wikipedia={d['wiki_kind']} "
              f"image={'yes' if rels[name] else 'no'} genome={d['genome'].get('size_mb')}Mb", flush=True)
        items.append(d)

    (out / "index.html").write_text(index_page(items, rels, args.title), encoding="utf-8")
    for i, d in enumerate(items):
        prev = items[i - 1]["name"] if i else None
        nxt = items[i + 1]["name"] if i + 1 < len(items) else None
        (out / f"{slug(d['name'])}.html").write_text(species_page(d, rels[d["name"]], prev, nxt, args.title), encoding="utf-8")
    print(f"Wrote {len(items) + 1} pages to {out.resolve()}")


if __name__ == "__main__":
    main()
