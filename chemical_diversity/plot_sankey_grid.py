"""Sankey diagrams of the scaffold chemotypes retained, lost and new along the AL rounds.

One panel per target (columns) and bias (rows).
Run compute_chemotypes.py first.
"""
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.path import Path as MplPath

# Define paths, targets, biases and colors for the Sankey diagrams
OUT = Path(__file__).resolve().parent / "results"
TARGETS = {"IL1b": "IL-1β", "IL2": "IL-2", "IL36g": "IL-36γ"}
BIASES = {"unbias": "Unbias", "bias1": "Bias 1", "bias2": "Bias 2"}
C = {"init": "#7B84D1", "retained": "#E8837A", "new1": "#59C3A0", "new2": "#F4B46E", "lost1": "#5CC8D1",
     "lost2": "#F07FC0", "gray": "#CCCCCC", "ret_rib1": "#A8B8E8", "ret_rib2": "#8FD4A0"}
BAR_W = 0.045  # half width of the node bars
GAP = 0.03     # vertical gap between stacked nodes

# Functions to draw the Sankey diagrams
def ribbon(ax, x0, x1, src, dst, color, alpha):
    """Bezier band joining the source (bottom, top) interval to the destination one."""
    cx = (x0 + x1) / 2
    verts = [(x0, src[0]), (cx, src[0]), (cx, dst[0]), (x1, dst[0]),
             (x1, dst[1]), (cx, dst[1]), (cx, src[1]), (x0, src[1]), (x0, src[0])]
    codes = [MplPath.MOVETO] + [MplPath.CURVE4] * 3 + [MplPath.LINETO] + [MplPath.CURVE4] * 3 + [MplPath.CLOSEPOLY]
    ax.add_patch(mpatches.PathPatch(MplPath(verts, codes), facecolor=color, edgecolor="none", alpha=alpha))


def draw_sankey(ax, r, scale):
    h = lambda v: v * scale
    ret1, lost1, new1 = r["Initial_to_Pred1_Retained"], r["Initial_to_Pred1_Lost"], r["Pred1_New_Scaffolds"]
    ret2, lost2, new2 = r["Pred1_to_Pred2_Retained"], r["Pred1_to_Pred2_Lost"], r["Pred2_New_Scaffolds"]
    init_top = h(r["Initial_Scaffolds"])

    # Nodes stacked bottom to top: lost, new, retained
    def stack(lost, new, ret):
        return {"lost": (0, h(lost)), "new": (h(lost) + GAP, h(lost + new) + GAP),
                "retained": (h(lost + new) + 2 * GAP, h(lost + new + ret) + 2 * GAP)}
    n1, n2 = stack(lost1, new1, ret1), stack(lost2, new2, ret2)

    bars = [(0, (0, init_top), C["init"]), (1, n1["lost"], C["lost1"]), (1, n1["new"], C["new1"]),
            (1, n1["retained"], C["retained"]), (2, n2["lost"], C["lost2"]), (2, n2["new"], C["new2"]),
            (2, n2["retained"], C["retained"])]
    for x, (bottom, top), color in bars:
        ax.fill_betweenx([bottom, top], x - BAR_W, x + BAR_W, color=color, lw=0, zorder=4)

    ribbon(ax, BAR_W, 1 - BAR_W, (init_top - h(ret1), init_top), n1["retained"], C["ret_rib1"], 0.7)
    ribbon(ax, BAR_W, 1 - BAR_W, (init_top - h(ret1 + new1), init_top - h(ret1)), n1["new"], C["gray"], 0.35)
    ribbon(ax, BAR_W, 1 - BAR_W, (0, h(lost1)), n1["lost"], C["gray"], 0.35)
    ribbon(ax, 1 + BAR_W, 2 - BAR_W, n1["retained"], n2["retained"], C["ret_rib2"], 0.7)
    ribbon(ax, 1 + BAR_W, 2 - BAR_W, (n1["new"][1] - h(new2), n1["new"][1]), n2["new"], C["new2"], 0.45)
    ribbon(ax, 1 + BAR_W, 2 - BAR_W, (n1["new"][0], n1["new"][0] + h(lost2)), n2["lost"], C["gray"], 0.35)

    # Total counts on top and percentages next to each node
    for x, nodes, total in [(0, {"": (0, init_top)}, r["Initial_Scaffolds"]),
                            (1, n1, r["Pred1_Scaffolds"]), (2, n2, r["Pred2_Scaffolds"])]:
        top = max(t for _, t in nodes.values())
        ax.text(x, top + 0.02, f"{total:,}", ha="center", va="bottom", fontsize=6.5, fontweight="bold")
    labels = {1: [("retained", "Retained", "Initial_to_Pred1_Retained_%"), ("new", "New", "Pred1_New_Scaffolds_%"),
                  ("lost", "Lost", "Initial_to_Pred1_Lost_%")],
              2: [("retained", "Retained", "Pred1_to_Pred2_Retained_%"), ("new", "New", "Pred2_New_Scaffolds_%"),
                  ("lost", "Lost", "Pred1_to_Pred2_Lost_%")]}
    for x, nodes in [(1, n1), (2, n2)]:
        for node, label, col in labels[x]:
            bottom, top = nodes[node]
            if top - bottom > 0.01 or node == "lost":
                ax.text(x + BAR_W + 0.03, (bottom + top) / 2, f"{label} ({r[col]:.1f}%)",
                        va="center", fontsize=5.5, color="#555")

    ax.set_xlim(-0.2, 2.6)
    ax.set_ylim(-0.05, 1.8)
    ax.set_xticks([0, 1, 2])
    ax.set_xticklabels(["Initial", "Pred 1", "Pred 2"], fontsize=7)
    ax.set_yticks([])
    ax.tick_params(axis="x", length=0)
    for side in ["top", "right", "left"]:
        ax.spines[side].set_visible(False)

# Load data (previously computed with compute_chemotypes.py)
df = pd.read_csv(OUT / "scaffold_chemotype_evolution.csv")
scale = 1 / df["Initial_Scaffolds"].max()

# Create figure
fig, axes = plt.subplots(3, 3, figsize=(10, 8), gridspec_kw={"hspace": 0.2, "wspace": 0.2})
for i, (bias, bias_label) in enumerate(BIASES.items()):
    for j, (target, target_label) in enumerate(TARGETS.items()):
        ax = axes[i, j]
        row = df[(df["Target"] == target) & (df["Bias"] == bias)]
        if row.empty:
            ax.set_visible(False)
            continue
        draw_sankey(ax, row.iloc[0], scale)
        if j == 0 or not axes[i, 0].get_visible():
            ax.set_ylabel(bias_label, fontsize=9, color="#444", style="italic")
        if i == 0:
            ax.set_title(target_label, fontsize=11, fontweight="bold")

legend = [("init", "Initial scaffolds"), ("retained", "Retained"), ("new1", "New (Pred 1)"),
          ("new2", "New (Pred 2)"), ("lost1", "Lost after Pred 1"), ("lost2", "Lost after Pred 2")]
fig.legend(handles=[mpatches.Patch(facecolor=C[k], label=l) for k, l in legend], loc="upper center", ncol=6,
           bbox_to_anchor=(0.5, 0.975), fontsize=8, frameon=False)

fig.savefig(OUT / "sankey_grid.png", dpi=600, bbox_inches="tight")
fig.savefig(OUT / "sankey_grid.pdf", bbox_inches="tight")
