"""Optional vector figure; no dependency on plotting to run the audit."""
from itertools import product
import numpy as np


def core_heatmap(out, records, preview=None):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    models = ["cq", "ca", "uq", "ua", "nq", "na"]
    labels = ["Censored QUAIDS + demographics", "Censored AIDS + demographics",
              "Uncensored QUAIDS + demographics", "Uncensored AIDS + demographics",
              "Uncensored QUAIDS, no demographics", "Uncensored AIDS, no demographics"]
    columns = list(product(("nls", "fgnls", "ifgnls"), ("gn", "lm"), ("zero", "linear")))
    byid = {r["id"]: r for r in records}
    data = np.full((6, 12), 3, dtype=int)
    names = {"PASS": 0, "NONCONVERGED": 1, "FAIL": 2}
    for i, model in enumerate(models):
        for j, (method, algorithm, start) in enumerate(columns):
            case = byid.get(f"core-{model}-{method}-{algorithm}-{start}")
            if case:
                data[i, j] = names.get(case["status"], 3)
    matplotlib.rcParams["svg.hashsalt"] = "pyquaidsce-small4-160"
    fig, ax = plt.subplots(figsize=(13.5, 4.8))
    ax.imshow(data, cmap=ListedColormap(["#d9efdf", "#ffe3ab", "#f4c4c4", "#eeeeee"]), vmin=0, vmax=3, aspect="auto")
    ax.set_yticks(range(6), labels)
    ax.set_xticks(range(12), [f"{m.upper()}\n{a.upper()} / {s}" for m, a, s in columns], fontsize=9)
    for i in range(6):
        for j in range(12):
            ax.text(j, i, ["Pass", "No conv.", "Fail", "Not run"][data[i, j]],
                    ha="center", va="center", fontsize=9, color="#253237")
    ax.set_xticks(np.arange(-.5, 12), minor=True)
    ax.set_yticks(np.arange(-.5, 6), minor=True)
    ax.grid(which="minor", color="white", linewidth=2)
    ax.tick_params(which="both", length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_title("pyquaidsce 1.6.0 · small4 · 72 core cases", loc="left", pad=20, fontsize=14, weight="bold")
    fig.text(.015, .015, "Default numerical controls. No convergence is preserved as a failed audit outcome.\n"
             "Pass covers the case's recorded assertions; see separate inference and interface findings.", fontsize=9, color="#455a64")
    fig.tight_layout(rect=(0, .08, 1, 1))
    fig.savefig(out / "core-matrix.svg", metadata={"Date": None}, bbox_inches="tight")
    svg = out / "core-matrix.svg"
    svg.write_text("\n".join(line.rstrip() for line in svg.read_text().splitlines()) + "\n", encoding="utf-8")
    if preview is not None:
        fig.savefig(preview, dpi=140, bbox_inches="tight")
    plt.close(fig)
