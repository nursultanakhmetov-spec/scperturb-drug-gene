import anndata as ad
import numpy as np

a = ad.read_h5ad("K562_gwps_normalized_bulk_01.h5ad")

with open("replogle_inspect.txt", "w", encoding="utf-8") as out:
    def p(*args):
        print(*args)
        print(*args, file=out)

    p("shape (perturbations x genes):", a.shape)
    p("obs columns:", list(a.obs.columns))
    p("var columns:", list(a.var.columns))
    p("layers:", list(a.layers.keys()))

    p("\nfirst rows of obs:")
    p(a.obs.head(5).to_string())

    p("\nfirst rows of var:")
    p(a.var.head(5).to_string())

    p("\nexample row names:", list(a.obs_names[:5]))
    p("rows containing 'non-targeting':", sum("non-targeting" in s for s in a.obs_names))

    X = np.asarray(a.X)
    p("X values (min, median, max):", float(X.min()), float(np.median(X)), float(X.max()))

    # check that the genes we want to compare against are present as targets
    wanted = ["ABL1", "BCR", "HDAC1", "HDAC2", "HDAC3", "AURKA", "AURKB",
              "PARP1", "JAK1", "JAK2", "BRD4", "EZH2", "DNMT1", "MTOR",
              "PLK1", "CDK1", "CDK9", "HSP90AB1", "MAP2K1", "MAP2K2"]
    p("\ntarget genes found in row names:")
    for g in wanted:
        hits = [s for s in a.obs_names if f"_{g}_" in s]
        p(f"  {g}: {len(hits)} row(s)", hits[:2])