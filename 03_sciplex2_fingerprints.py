import anndata as ad, scanpy as sc, pandas as pd, numpy as np

a = ad.read_h5ad("SrivatsanTrapnell2020_sciplex2.h5ad")
a.var_names_make_unique()
a.obs["dose"] = pd.to_numeric(a.obs["dose_value"], errors="coerce")

# убрать неприписанные клетки ("control") и ненадёжно приписанные
ok = (a.obs["perturbation"].astype(str) != "control") & a.obs["dose"].notna() \
     & (a.obs["qval_demultiplexing"] <= 0.05)
a = a[ok.values].copy()
print("клеток после фильтра:", a.n_obs)

sc.pp.filter_genes(a, min_cells=50)
sc.pp.normalize_total(a, target_sum=1e4)
sc.pp.log1p(a)

drug = a.obs["perturbation"].astype(str).values
dose = a.obs["dose"].values
rows, names = [], []
for d in sorted(set(drug)):
    base = np.asarray(a.X[(drug == d) & (dose == 0)].mean(axis=0)).ravel()   # растворитель этого препарата
    for x in sorted(set(dose[drug == d])):
        if x == 0: continue
        m = np.asarray(a.X[(drug == d) & (dose == x)].mean(axis=0)).ravel()
        rows.append(m - base); names.append(f"{d}_{x:g}")

fp = pd.DataFrame(rows, index=names, columns=a.var_names)
fp.to_csv("sciplex2_fingerprints.csv.gz")
print("отпечатков:", fp.shape[0], "генов:", fp.shape[1])

check = [g for g in ["ANGPTL4", "GDF15", "CDKN1A", "MDM2"] if g in fp.columns]
with open("sciplex2_check.txt", "w", encoding="utf-8") as out:
    t = fp[check].round(3).to_string()
    print(t); print(t, file=out)