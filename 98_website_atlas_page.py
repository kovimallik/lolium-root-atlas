#!/usr/bin/env python
"""Static root-atlas page for sustainplantgroup.com (GitHub Pages; no server, no file > 100 MB).

Writes results/browser/website_root_atlas/:
  index.html                 the page (schematic + cell-state UMAP + gene box -> feature plot drawn in the browser)
  schematic.svg              cell-by-cell root schematic (figures/Root_schematic_Lp.svg)
  data/umap.png              2 x N RGB PNG: row 0 = x, row 1 = y as (R = high byte, G = low byte) of 0..65535-scaled UMAP, 20,000-nucleus subsample
  data/labels.png            1 x N grey PNG: cell-state index
  data/libs.png              1 x N grey PNG: library index into meta.json['libs'] (genotype, N level) into data/meta.json["states"]
  data/meta.json             states (name, colour, n nuclei), counts
  data/genes.json            gene index: id -> [chunk, row, max, pct nuclei, description, role]
  data/chunk_<k>.png         32 x N grey PNG: genes x nuclei, uint8 (log1p CP10k scaled to the gene's 99th percentile), fetched on demand
Genes: every gene detected in >= MIN_NUCLEI nuclei (expression chunks) + the rest in the index only; roles annotate the NUE panel (ref/nue_panel_v2.tsv), the 40 strongest markers of each Leiden cluster (results/clean/markers_leiden1
aggregated by state) and the 80 most expressed tau >= 0.6 specific genes of each state (sn_label_specificity.tsv), plus a few named genes.
Usage: 98_website_atlas_page.py
"""
import json, os, sys
import numpy as np, pandas as pd
sys.path.insert(0, "/home/kovimallik/Desktop/Lp_NUE_spatial/04_manuscript/figures/scripts")
import style as S

OUT = f"{S.A}/results/browser/website_root_atlas"; os.makedirs(f"{OUT}/data", exist_ok=True)
N_SUB, CHUNK, SEED = 20000, 128, 0
MIN_NUCLEI = 10            # genes detected in fewer nuclei are listed in the index but have no expression chunk
rng = np.random.default_rng(SEED)

obs = S.load_obs(["label_from_audit", "sample"]); U = S.load_umap()
grp = S.merge_other(obs["group"].astype(str).values)
states = [g for g in S.GROUP_ORDER_ANATOMICAL]
idx = np.sort(rng.choice(len(obs), N_SUB, replace=False))
Us = U[idx]; lab = np.array([states.index(g) for g in grp[idx]], dtype=np.uint8)
lo, hi = Us.min(0), Us.max(0)
Uq = np.round((Us - lo) / (hi - lo) * 65535).astype(np.uint16)
# PNGs rather than raw binaries: static hosts and the artifact viewer serve images; the page reads them back via canvas
from PIL import Image
um = np.zeros((2, N_SUB, 3), dtype=np.uint8)
um[0, :, 0], um[0, :, 1] = Uq[:, 0] >> 8, Uq[:, 0] & 255
um[1, :, 0], um[1, :, 1] = Uq[:, 1] >> 8, Uq[:, 1] & 255
Image.fromarray(um, "RGB").save(f"{OUT}/data/umap.png", optimize=True)
Image.fromarray(lab[None, :], "L").save(f"{OUT}/data/labels.png", optimize=True)
samp = obs["sample"].astype(str).values[idx]
Image.fromarray(np.array([S.SAMPLE_ORDER.index(x) for x in samp], dtype=np.uint8)[None, :], "L").save(f"{OUT}/data/libs.png", optimize=True)

# ---- gene set
ann = S.gene_annot()
panel = pd.read_csv(f"{S.A}/ref/nue_panel_v2.tsv", sep="\t")
genes = {g: f"NUE panel: {f}" for g, f in zip(panel.LOC, panel.fam)}
spec = pd.read_csv(f"{S.A}/results/bulk_timecourse/sn_label_specificity.tsv", sep="\t")
spec = spec[(spec.tau >= 0.6) & (spec.top_cpm >= 5)].sort_values("top_cpm", ascending=False).groupby("top_label").head(80)
for g, l in zip(spec.gene, spec.top_label):
    genes.setdefault(g, f"specific to {S.canon_group(l)}")
mk = pd.read_csv(f"{S.A}/results/clean/markers_leiden1.tsv", sep="\t")
cl2grp = obs.groupby("label_from_audit", observed=True).size()  # not needed; markers are per leiden cluster
mk = mk[(mk.pvals_adj < 1e-10) & (mk.logfoldchanges > 1.5)].sort_values("scores", ascending=False)
for g in mk.groupby("group").head(40).names:
    genes.setdefault(g, "cluster marker")
ish = pd.read_csv("/home/kovimallik/Desktop/Lp_NUE_spatial/01_design/ISH_probe_sequences/ISH_probe_targets.tsv", sep="\t")
for g, t in zip(ish.gene, ish.probe_target): genes.setdefault(g, f"ISH target: {t.split('_')[0]}")
for g in ["LOC127335392", "LOC127307575", "LOC127326287", "LOC127305389", "LOC127343520", "LOC127308872", "LOC127297367", "LOC127297370", "LOC127297389"]:
    genes.setdefault(g, "named gene")
var = S.load_var_names(); vset = set(var)
roles = {g: r for g, r in genes.items() if g in vset}
X = S.load_counts_csr(); totals = np.asarray(X.sum(axis=1)).ravel()
nnz = np.diff(X.tocsc().indptr); nnz_by_gene = dict(zip(var, nnz))
glist = sorted(g for g, n in zip(var, nnz) if n >= MIN_NUCLEI)
genes = {g: roles.get(g, "") for g in glist}
rare = [g for g, n in zip(var, nnz) if n < MIN_NUCLEI]
print(len(glist), "genes")

# ---- expression chunks
meta_genes = {}
for k in range(0, len(glist), CHUNK):
    gs = glist[k:k + CHUNK]
    E, _, _ = S.gene_expression(gs, X=X, totals=totals)
    block = np.zeros((len(gs), N_SUB), dtype=np.uint8)
    for i, g in enumerate(gs):
        v = E[g].values[idx]
        vmax = float(np.percentile(v[v > 0], 99)) if (v > 0).any() else 1.0
        block[i] = np.clip(np.round(v / vmax * 255), 0, 255).astype(np.uint8)
        d = str(ann.desc.get(g, ""))
        meta_genes[g] = [k // CHUNK, i, round(vmax, 3), round(float((v > 0).mean()) * 100, 1), d[:90], genes[g]]
    Image.fromarray(block, "L").save(f"{OUT}/data/chunk_{k // CHUNK}.png", optimize=True)
    if (k // CHUNK) % 10 == 0: print("chunk", k // CHUNK, flush=True)

counts = pd.Series(grp).value_counts()
meta = dict(n=N_SUB, n_total=int(len(obs)), umap_range=[lo.tolist(), hi.tolist()], chunk=CHUNK, n_genes=len(glist), n_rare=len(rare),
            states=[dict(name=s, label=S.GROUP_LEGEND[s], color=S.GROUP_COLORS[s], n=int(counts.get(s, 0))) for s in states],
            libs=[dict(id=l, geno=S.GENOTYPE_NAMES[l.split("_")[1]], trt=S.TREATMENT_NAMES[l.split("_")[2]]) for l in S.SAMPLE_ORDER],
            geno_colors={S.GENOTYPE_NAMES[k]: v for k, v in S.GENOTYPE_COLORS.items()})
json.dump(meta, open(f"{OUT}/data/meta.json", "w"), separators=(",", ":"))
# gene index: id -> [chunk, row, max, pct, description, role]; rare genes -> [-1, nnz, 0, 0, description, ""]
for g in rare: meta_genes[g] = [-1, int(nnz_by_gene[g]), 0, 0, str(ann.desc.get(g, ""))[:90], ""]
json.dump(meta_genes, open(f"{OUT}/data/genes.json", "w"), separators=(",", ":"))
svg = open("/home/kovimallik/Desktop/Lp_NUE_spatial/04_manuscript/figures/Root_schematic_Lp.svg").read()
import re as _re; open(f"{OUT}/schematic.svg", "w").write(_re.sub(r"<!DOCTYPE[^>]*>\s*", "", svg, flags=_re.S))  # no DTD: static hosts and the artifact viewer refuse it
print("wrote", OUT, "size MB", sum(os.path.getsize(os.path.join(dp, f)) for dp, _, fs in os.walk(OUT) for f in fs) / 1e6)
