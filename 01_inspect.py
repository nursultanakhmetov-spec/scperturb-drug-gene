import anndata as ad

FILES = {
    "nadig":    "NadigOConner2024_hepg2.h5ad",
    "norman":   "NormanWeissman2019_filtered.h5ad",
    "sciplex2": "SrivatsanTrapnell2020_sciplex2.h5ad",
    "sciplex3": "SrivatsanTrapnell2020_sciplex3.h5ad",
}

with open("inspect_output.txt", "w", encoding="utf-8") as out:
    def p(*args):
        print(*args)
        print(*args, file=out)

    for name, f in FILES.items():
        a = ad.read_h5ad(f, backed="r")          # матрица в память не грузится
        p("\n" + "=" * 70)
        p(name, f, "клеток x генов:", a.shape)
        p("obs:", list(a.obs.columns))
        p("var:", list(a.var.columns))
        p("layers:", list(a.layers.keys()))
        p("obsm:", list(a.obsm.keys()))
        for c in a.obs.columns:
            n = a.obs[c].nunique()
            if n < 300:
                p(f"  [{c}] уникальных: {n}")
                p("   ", a.obs[c].value_counts().head(10).to_dict())
            else:
                p(f"  [{c}] уникальных: {n} (примеры: {list(a.obs[c].unique()[:5])})")
        p(a.var.head(5).to_string())
        a.obs.to_csv(f"{name}_obs.csv.gz")       # таблица клеток
        a.var.to_csv(f"{name}_var.csv")          # таблица генов
        a.file.close()
        out.flush()                              # результат по файлу сразу пишется на диск