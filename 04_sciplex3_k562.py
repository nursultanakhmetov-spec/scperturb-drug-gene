import anndata as ad, scanpy as sc, pandas as pd, numpy as np

F = "SrivatsanTrapnell2020_sciplex3.h5ad"
a = ad.read_h5ad(F, backed="r")
mask = (a.obs["cell_line"].astype(str) == "K562").values
k = a[mask].to_memory()                 # в память берём только K562
a.file.close()
k.var_names_make_unique()
print("K562 клеток:", k.n_obs)

sc.pp.filter_genes(k, min_cells=200)
sc.pp.normalize_total(k, target_sum=1e4)
sc.pp.log1p(k)
print("генов после фильтра:", k.n_vars)

drug = k.obs["perturbation"].astype(str).values
dose = k.obs["dose_value"].astype(float).values
target = dict(zip(drug, k.obs["target"].astype(str).values))

base = np.asarray(k.X[drug == "control"].mean(axis=0)).ravel()
rows, meta = [], []
for d in sorted(set(drug) - {"control"}):
    for x in sorted(set(dose[drug == d])):
        sel = (drug == d) & (dose == x)
        rows.append(np.asarray(k.X[sel].mean(axis=0)).ravel() - base)
        meta.append((f"{d}|{x:g}", d, x, target[d], int(sel.sum())))

meta = pd.DataFrame(meta, columns=["name", "drug", "dose_nM", "target", "n_cells"]).set_index("name")
fp = pd.DataFrame(rows, index=meta.index, columns=k.var_names)
fp.to_csv("sciplex3_k562_fingerprints.csv.gz")
meta["effect"] = np.linalg.norm(fp.values, axis=1).round(2)      # сила отпечатка
meta.to_csv("sciplex3_k562_meta.csv")

with open("sciplex3_k562_check.txt", "w", encoding="utf-8") as out:
    def p(*args):
        print(*args); print(*args, file=out)

    top = meta[meta["dose_nM"] == 10000]
    p("=== 20 самых сильных препаратов на 10 000 нМ ===")
    p(top.sort_values("effect", ascending=False).head(20)[["target", "n_cells", "effect"]].to_string())
    p("\n=== медиана силы эффекта по дозам ===")
    p(meta.groupby("dose_nM")["effect"].median().to_string())
    p("\n=== ингибиторы Bcr-Abl по дозам ===")
    p(meta[meta["target"].str.contains("Bcr-Abl")][["dose_nM", "target", "n_cells", "effect"]].to_string())

    # похожи ли препараты одного класса друг на друга
    f = fp.loc[top.index]
    genes = f.var().sort_values().index[-2000:]
    c = pd.DataFrame(np.corrcoef(f[genes].values), index=top.index, columns=top.index)
    p("\n=== средняя корреляция отпечатков: внутри класса / с остальными ===")
    for cls in ["HDAC", "JAK", "Aurora Kinase", "PARP", "Sirtuin"]:
        inn = top.index[top["target"] == cls]; oth = top.index[top["target"] != cls]
        w = c.loc[inn, inn].values[np.triu_indices(len(inn), 1)].mean()
        b = c.loc[inn, oth].values.mean()
        p(f"{cls:15s} препаратов {len(inn):2d}   внутри {w:.2f}   с остальными {b:.2f}")