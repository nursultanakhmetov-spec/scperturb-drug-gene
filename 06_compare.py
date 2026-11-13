import anndata as ad
import numpy as np
import pandas as pd

# ---------- Replogle: gene knockdown profiles ----------
rep = ad.read_h5ad("K562_gwps_normalized_bulk_01.h5ad")
R = pd.DataFrame(np.asarray(rep.X), index=rep.obs_names, columns=rep.var_names)
R = R.replace([np.inf, -np.inf], np.nan)
R = R.loc[:, R.notna().all(axis=0)]                      # drop genes with inf/NaN
rep_gene = np.array([s.split("_")[1] for s in R.index])  # target gene of each row

# ---------- sciplex3 K562: drug fingerprints ----------
fp = pd.read_csv("sciplex3_k562_fingerprints.csv.gz", index_col=0)
meta = pd.read_csv("sciplex3_k562_meta.csv", index_col=0)
var = pd.read_csv("sciplex3_var.csv", index_col=0)

sym2ens = var.loc[~var.index.duplicated(), "ensembl_id"].astype(str).str.split(".").str[0]
fp = fp.loc[:, fp.columns.isin(sym2ens.index)]
fp.columns = sym2ens.loc[fp.columns].values
fp = fp.loc[:, ~fp.columns.duplicated()]

common = R.columns.intersection(fp.columns)

def zrows(M):
    M = M - M.mean(axis=1, keepdims=True)
    return M / M.std(axis=1, keepdims=True)

Rz = zrows(R[common].values)

# ---------- queries: (label, fingerprint vector, expected target genes) ----------
queries = []

def add_class(cls, genes, dose=10000):
    idx = meta.index[(meta["target"] == cls) & (meta["dose_nM"] == dose)]
    queries.append((f"CLASS {cls} (mean of {len(idx)} drugs, {dose} nM)",
                    fp.loc[idx, common].mean(axis=0).values, genes))

def add_drug(name_part, genes, doses=(1000, 10000)):
    for dose in doses:
        idx = meta.index[meta["drug"].str.contains(name_part, regex=False) & (meta["dose_nM"] == dose)]
        for i in idx:
            queries.append((i, fp.loc[i, common].values, genes))

add_drug("Nilotinib", ["BCR", "ABL1"])
add_drug("Dasatinib", ["BCR", "ABL1"])
add_class("HDAC", ["HDAC1", "HDAC2", "HDAC3"])
add_class("Aurora Kinase", ["AURKA", "AURKB"])
add_class("JAK", ["JAK1", "JAK2"])                 # expected negative control
add_drug("JQ1", ["BRD4"])
add_drug("Trametinib", ["MAP2K1", "MAP2K2"])
add_drug("Temsirolimus", ["MTOR"])
add_drug("Flavopiridol", ["CDK9", "CDK1"])
add_drug("Rigosertib", ["PLK1"])
add_drug("Tanespimycin", ["HSP90AB1"])
add_drug("Tazemetostat", ["EZH2"])
add_drug("Decitabine", ["DNMT1"])

# ---------- compare ----------
with open("compare_check.txt", "w", encoding="utf-8") as out:
    def p(*args):
        print(*args)
        print(*args, file=out)

    p("Replogle perturbations:", R.shape[0], "| genes after cleaning:", R.shape[1])
    p("sciplex3 fingerprint genes:", fp.shape[1], "| common genes used:", len(common))

    n = len(rep_gene)
    for label, vec, expected in queries:
        v = zrows(vec[None, :])[0]
        corr = Rz @ v / len(common)                 # Pearson r with every knockdown
        order = np.argsort(-corr)
        rank = np.empty(n, dtype=int)
        rank[order] = np.arange(1, n + 1)

        p("\n" + "=" * 70)
        p(label)
        for g in expected:
            hits = np.where(rep_gene == g)[0]
            for h in hits:
                p(f"  expected {g:9s} r = {corr[h]:+.3f}   rank {rank[h]:5d} of {n}"
                  f"   (top {100 * rank[h] / n:.1f}%)")
            if len(hits) == 0:
                p(f"  expected {g:9s} not found in Replogle")
        top = [f"{rep_gene[i]} ({corr[i]:+.2f})" for i in order[:10]]
        p("  top 10 matches:", ", ".join(top))