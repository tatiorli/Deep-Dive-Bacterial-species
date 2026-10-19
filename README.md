# Bacteria Atlas - Build easily bacterial focused websites to highlight species of interest and deep dive into them in seconds

Paste bacterial species names into a text file and get a static website: one page per species
plus an index with a search box and a genome-size chart. Standard-library Python only.

```
python build_site.py species.txt
python build_site.py species.txt --out docs --title "Bacteria Atlas" --email me@lab.org
```

Open `site/index.html`. Needs internet access when it runs. Results are cached in `.bacteria_cache/`
(`--refresh` refetches, `--no-images` skips pictures).

## What each page contains
- Lead picture (Wikipedia article image, else Wikimedia Commons, else the genus article's image, else a generic schematic) with credit and licence link
- Visible link buttons: Wikipedia, NCBI Taxonomy, NCBI genomes, PubMed, LPSN, BacDive, Scholar
- Wikipedia summary, lineage, former names, reference genome size/GC/accession, assembly count, PubMed paper count

## Input format
One `Genus species` per line; `#` comments and blank lines are skipped; log-style lines such as
`SKIPPED Cutibacterium acnes: no output yet` and `Escherichia_coli` are cleaned automatically.

## Sources
NCBI Taxonomy, NCBI Datasets, PubMed (E-utilities), Wikipedia, Wikimedia Commons. Pictures keep their
own licences; Wikipedia text is CC BY-SA. Verify facts before citing.
