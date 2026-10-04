# Perennial ryegrass root atlas - web viewer

Single-nucleus RNA-seq atlas of the *Lolium perenne* root (76,684 nuclei, 15 cell states; Barhoney and Fagerlin, full N and 25% less N). This repository is the static web page served by GitHub Pages: a cell-by-cell root schematic coloured by cell state, the cell-state UMAP, and a gene viewer that draws expression on the UMAP in the browser for every gene detected in ≥10 nuclei (35,145 genes; keyword search on gene descriptions).

- `index.html` - the page (no build step, no server; plain HTML/JS)
- `schematic.svg` - root schematic (longitudinal tip + mature cross-section)
- `data/` - 20,000-nucleus display subsample: `umap.png`, `labels.png`, `meta.json`, `genes.json` (gene index), `chunk_<k>.png` (128 genes × 20,000 nuclei, 8-bit, log1p CP10k scaled to each gene's 99th percentile)
- `98_website_atlas_page.py` - builder, run on the processed atlas (code and data: see the paper's Data availability)

## Licence and citation

CC BY 4.0 (see `LICENSE`). Any use of the viewer, the schematic or the display data must cite the article:

> Kovi, M. R. *et al.* A constitutively active cortical nitrate-uptake cell state underlies nitrogen-use efficiency in perennial ryegrass. *bioRxiv* (2026). DOI to be added.

A machine-readable citation is in `CITATION.cff` (GitHub shows a "Cite this repository" button from it). Until the article is posted, the data are unpublished: please contact the author before reuse.

Sustainable Plant Group, NMBU - https://sustainplantgroup.com
