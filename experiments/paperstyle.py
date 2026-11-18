"""Shared figure style for the manuscript.

Design intent: restrained, journal-typical, no decoration. Small type, thin
rules, two visible spines, a muted palette that survives greyscale printing,
and panel letters set outside the axes rather than as bold headings.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Muted, print-safe, clearly distinguishable in greyscale by lightness order.
NAVY   = "#27496D"
BRICK  = "#9E3A26"
FOREST = "#3B6E4F"
PLUM   = "#6C5B7B"
OCHRE  = "#B07D2B"
GREY   = "#595959"
LGREY  = "#B4B4B4"

CYCLE = [NAVY, BRICK, FOREST, PLUM, OCHRE, GREY]

METHOD = {
    "MC-ContrastiveVI": BRICK,
    "contrastiveVI":    NAVY,
    "multiGroupVI":     FOREST,
    "scVI":             PLUM,
    "scDisInFact":      OCHRE,
    "PCA(32)":          GREY,
    "PCA":              GREY,
}

BLOCK = {"Shared": NAVY, "Doxorubicin": BRICK, "Cisplatin": FOREST,
         "Background": GREY, "Drug": FOREST}


def apply():
    plt.rcParams.update({
        "figure.dpi": 300, "savefig.dpi": 300, "savefig.bbox": "tight",
        "savefig.pad_inches": 0.02,
        "font.family": "DejaVu Sans", "font.size": 7,
        "axes.labelsize": 7, "axes.titlesize": 7,
        "xtick.labelsize": 6.5, "ytick.labelsize": 6.5,
        "legend.fontsize": 6.5, "legend.frameon": False,
        "legend.handletextpad": 0.4, "legend.columnspacing": 1.1,
        "legend.borderaxespad": 0.2, "legend.handlelength": 1.1,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.linewidth": 0.6, "axes.labelpad": 2.5,
        "axes.prop_cycle": matplotlib.cycler(color=CYCLE),
        "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "xtick.major.size": 2.4, "ytick.major.size": 2.4,
        "xtick.direction": "out", "ytick.direction": "out",
        "xtick.major.pad": 2, "ytick.major.pad": 2,
        "lines.linewidth": 0.9, "lines.markersize": 3,
        "grid.linewidth": 0.4, "grid.color": "#E2E2E2",
        "pdf.fonttype": 42, "ps.fonttype": 42,
    })


def panel(ax, letter, title=None, dx=-0.16, dy=1.06):
    """Panel letter outside the axes; optional plain-weight descriptive title."""
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=8,
            fontweight="bold", va="bottom", ha="left")
    if title:
        ax.set_title(title, loc="left", fontsize=7, pad=4)


def faint_grid(ax, axis="y"):
    ax.grid(axis=axis, lw=0.4, color="#E8E8E8", zorder=0)
    ax.set_axisbelow(True)
