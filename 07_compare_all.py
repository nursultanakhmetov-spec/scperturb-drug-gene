import anndata as ad
import numpy as np
import pandas as pd

P_THRESH = 0.05          # Replogle energy-test threshold for "knockdown has an effect"

# ---- my assumption: target class from sciplex3 -> concrete genes (edit freely) ----
TOKEN2GENES = {
    "ALK": ["ALK"], "c-Met": ["MET"], "Androgen Receptor": ["AR"], "Aromatase": ["CYP19A1"],
    "Aurora Kinase": ["AURKA", "AURKB"], "Bcr-Abl": ["BCR", "ABL1"], "FLT3": ["FLT3"],
    "JAK": ["JAK1", "JAK2", "JAK3", "TYK2"], "c-RET": ["RET"],
    "FGFR": ["FGFR1", "FGFR2", "FGFR3"], "CDK": ["CDK1", "CDK2", "CDK4", "CDK6", "CDK7", "CDK9"],
    "VEGFR": ["KDR", "FLT1", "FLT4"], "ROCK": ["ROCK1", "ROCK2"],
    "Bcl-2": ["BCL2", "BCL2L1", "MCL1"], "c-Kit": ["KIT"], "Src": ["SRC", "LYN"],
    "Gamma-secretase": ["PSEN1", "PSENEN", "NCSTN", "APH1A"], "CCR": ["CCR5"],
    "COX": ["PTGS1", "PTGS2"], "CSF-1R": ["CSF1R"], "PDGFR": ["PDGFRA", "PDGFRB"],
    "DNA Methyltransferase": ["DNMT1", "DNMT3A", "DNMT3B"], "Dehydrogenase": ["ALDH2", "ALDH1A1"],
    "Dopamine Receptor": ["DRD2", "DRD3"], "EGFR": ["EGFR"], "HER2": ["ERBB2"],
    "HDAC": ["HDAC1", "HDAC2", "HDAC3", "HDAC6", "HDAC8"],
    "Epigenetic Reader Domain": ["BRD4", "BRD2", "BRD3", "EED"],
    "Estrogen/progestogen Receptor": ["ESR1"], "FAAH": ["FAAH"], "FAK": ["PTK2"],
    "HIF": ["EGLN1"], "HSP (e.g. HSP90)": ["HSP90AA1", "HSP90AB1"],
    "Histone Acetyltransferase": ["EP300", "KAT2B"],
    "Histone Demethylase": ["KDM1A", "KDM6A", "KDM6B"],
    "Histone Methyltransferase": ["EZH2", "EHMT1", "EHMT2", "KMT5A"],
    "IDH1": ["IDH1"], "IGF-1R": ["IGF1R"], "LPA Receptor": ["LPAR1"], "Lipoxygenase": ["ALOX5"],
    "MAO": ["MAOA", "MAOB", "KDM1A"], "MEK": ["MAP2K1", "MAP2K2"],
    "Microtubule Associated": ["TUBB", "TUBA1B"], "PARP": ["PARP1", "PARP2"],
    "Raf": ["RAF1", "BRAF"], "PI3K": ["PIK3CA", "PIK3CB", "PIK3CD"],
    "PKC": ["PRKCA", "PRKCB", "PRKCD"], "PLK": ["PLK1"], "Pim": ["PIM1", "PIM2"],
    "STAT": ["STAT3"], "Survivin": ["BIRC5"], "TGF-beta/Smad": ["TGFBR1"],
    "Telomerase": ["TERT"], "Tie-2": ["TEK"], "Topoisomerase": ["TOP1", "TOP2A", "TOP2B"],
    "Wnt/beta-catenin": ["TNKS", "TNKS2"], "mTOR": ["MTOR", "RPTOR"],
    "Trk receptor": ["NTRK1"], "DNA/RNA Synthesis": ["TYMS"],
}
# classes left unmapped on purpose: agonists/activators or no single protein target
# (AMPK, Glucocorticoid Receptor, Sirtuin, E3 Ligase, p53, TNF-alpha, DNA alkylator,
#  Autophagy, GABA Receptor, Beta Amyloid, NF-kB, Nrf2, MT Receptor, Histamine Receptor, Others)

def expected_genes(target):
    genes = []
    for tok in str(target).split(","):
        genes += TOKEN2GENES.get(tok.strip(), [])
    return sorted(set(genes))

def zrows(M):
    M = M - M.mean(axis=1, keepdims=True)
    s = M.std(axis=1, keepdims=True)
    s[s == 0] = 1.0
    return M / s

out = open("compare_all_summary.txt", "w", encoding="utf-8")
def p(*args):
    print(*args)
    print(*args, file=out)

# ================= Replogle =================
rep = ad.read_h5ad("K562_gwps_normalized_bulk_01.h5ad")
R = pd.DataFrame(np.asarray(rep.X), index=rep.obs_names, columns=rep.var_names)
R = R.replace([np.inf, -np.inf], np.nan)
R = R.loc[:, R.notna().all(axis=0)]

is_ctrl = R.index.str.contains("non-targeting")
pval = rep.obs["energy_test_p_value"].reindex(R.index)
keep = (pval < P_THRESH).fillna(False).values & ~is_ctrl
p("Replogle rows:", R.shape[0], "| non-targeting:", int(is_ctrl.sum()),
  f"| energy-test p < {P_THRESH}:", int(keep.sum()))
if keep.sum() < 500:                                  # fallback if the p-value filter is too strict
    strength = R.std(axis=1).values
    strength[is_ctrl] = -1
    keep = np.zeros(len(R), dtype=bool)
    keep[np.argsort(-strength)[:3000]] = True
    p("too few rows passed -> fallback: 3000 strongest knockdown profiles")
R = R.loc[keep]
rep_gene = np.array([s.split("_")[1] for s in R.index])
p("Replogle rows kept:", R.shape[0], "| unique genes:", len(set(rep_gene)))

# ================= sciplex3 K562 =================
fp = pd.read_csv("sciplex3_k562_fingerprints.csv.gz", index_col=0)
meta = pd.read_csv("sciplex3_k562_meta.csv", index_col=0).loc[fp.index]
var = pd.read_csv("sciplex3_var.csv", index_col=0)
sym2ens = var.loc[~var.index.duplicated(), "ensembl_id"].astype(str).str.split(".").str[0]
fp = fp.loc[:, fp.columns.isin(sym2ens.index)]
fp.columns = sym2ens.loc[fp.columns].values
fp = fp.loc[:, ~fp.columns.duplicated()]

common = R.columns.intersection(fp.columns)
p("common genes used:", len(common))

# ---- improvement 1: remove shared background on the drug side (per dose, median over drugs)
D = fp[common].copy()
for dose in sorted(meta["dose_nM"].unique()):
    idx = meta.index[meta["dose_nM"] == dose]
    D.loc[idx] = D.loc[idx] - D.loc[idx].median(axis=0)
meta["effect_corrected"] = np.linalg.norm(D.values, axis=1)

# ---- same idea on the gene side: center every gene across kept knockdowns
Rc = R[common] - R[common].mean(axis=0)

# ================= compare everything =================
Rz = zrows(Rc.values)
Dz = zrows(D.values)
corr = Dz @ Rz.T / len(common)                         # drugs x knockdowns
order = np.argsort(-corr, axis=1)
rank = np.empty_like(order)
rows = np.arange(corr.shape[0])[:, None]
rank[rows, order] = np.arange(1, corr.shape[1] + 1)
n = corr.shape[1]

records = []
for i, name in enumerate(D.index):
    exp = expected_genes(meta.loc[name, "target"])
    hits = np.where(np.isin(rep_gene, exp))[0]
    rec = {"name": name, "drug": meta.loc[name, "drug"], "dose_nM": meta.loc[name, "dose_nM"],
           "target": meta.loc[name, "target"], "n_cells": meta.loc[name, "n_cells"],
           "effect_corrected": round(float(meta.loc[name, "effect_corrected"]), 2),
           "expected_genes": ";".join(exp), "n_expected_rows": len(hits),
           "top5_matches": ", ".join(f"{rep_gene[j]} ({corr[i, j]:+.2f})" for j in order[i, :5])}
    if len(hits):
        best = hits[np.argmin(rank[i, hits])]
        rec.update(best_gene=rep_gene[best], best_rank=int(rank[i, best]),
                   best_pct=round(100 * rank[i, best] / n, 2), best_r=round(float(corr[i, best]), 3))
    records.append(rec)
T = pd.DataFrame(records).set_index("name")
T.to_csv("compare_all_table.csv")

# ================= summary =================
def chance(k, q):                                      # P(best of k random rows lands in top q)
    return 1 - (1 - q) ** k

p("\n=== HIT RATES (target gene among best matches) ===")
p("drugs with at least one expected gene present in kept Replogle rows are scored")
for dose in sorted(T["dose_nM"].unique()):
    S = T[(T["dose_nM"] == dose) & T["best_rank"].notna()]
    act_cut = T.loc[T["dose_nM"] == dose, "effect_corrected"].quantile(0.75)
    for label, sub in [("all scored drugs", S), ("strongest 25% drugs", S[S["effect_corrected"] >= act_cut])]:
        if len(sub) == 0:
            continue
        p(f"dose {dose:>7g} nM | {label:20s} n = {len(sub):3d} | "
          f"top 1%: {100 * (sub['best_pct'] <= 1).mean():5.1f}% (chance {100 * chance(sub['n_expected_rows'], 0.01).mean():4.1f}%) | "
          f"top 5%: {100 * (sub['best_pct'] <= 5).mean():5.1f}% (chance {100 * chance(sub['n_expected_rows'], 0.05).mean():4.1f}%)")

p("\n=== BY TARGET CLASS, dose 1000 nM (classes with >= 3 scored drugs) ===")
S = T[(T["dose_nM"] == 1000) & T["best_rank"].notna()]
g = S.groupby("target").agg(n=("best_pct", "size"), median_pct=("best_pct", "median"),
                            in_top5=("best_pct", lambda x: int((x <= 5).sum())))
p(g[g["n"] >= 3].sort_values("median_pct").to_string())

p("\n=== ALL HITS IN TOP 1% (any dose) ===")
H = T[T["best_pct"] <= 1].sort_values("best_pct")
p(H[["target", "best_gene", "best_rank", "best_r", "effect_corrected"]].to_string())

p("\n=== SANITY ROWS ===")
for part in ["Nilotinib", "Dasatinib", "Temsirolimus", "Trametinib", "Ruxolitinib (INCB018424)|", "Panobinostat"]:
    sub = T[T.index.str.contains(part, regex=False)]
    p(sub[["best_gene", "best_rank", "best_pct", "top5_matches"]].to_string())

n_unmapped = T.loc[T["dose_nM"] == 1000, "expected_genes"].eq("").sum()
n_absent = ((T["dose_nM"] == 1000) & T["expected_genes"].ne("") & T["best_rank"].isna()).sum()
p(f"\ndrugs without a mapped gene target: {n_unmapped} | mapped but target not in kept Replogle rows: {n_absent}")
out.close()