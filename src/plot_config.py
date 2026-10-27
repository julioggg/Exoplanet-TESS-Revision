"""
============================================================
PLOT CONFIGURATION
============================================================

Common plotting style for all figures in the manuscript.

This module centralizes the graphical configuration used
throughout the project in order to obtain publication-quality
figures suitable for scientific journals.

Author:
Bryan Cabrera

Project:
Multi-objective Optimization for Exoplanet Detection
============================================================
"""

import matplotlib.pyplot as plt


# ============================================================
# GLOBAL STYLE
# ============================================================

FIG_WIDTH = 8
FIG_HEIGHT = 6

DPI = 600

TITLE_SIZE = 20
LABEL_SIZE = 15
TICK_SIZE = 12
LEGEND_SIZE = 12

LINE_WIDTH = 2.5

GRID_ALPHA = 0.30

# ============================================================
# SCIENTIFIC COLOR PALETTE
# ============================================================

COLORS = {

    "blue": "#0072B2",

    "orange": "#D55E00",

    "green": "#009E73",

    "purple": "#CC79A7",

    "yellow": "#F0E442",

    "gray": "#666666"

}

def get_color(name):

    """
    Return a scientific color from the project palette.
    """

    return COLORS[name]

def setup_plot():
    """
    Configure matplotlib using a common style for every figure.
    """

    plt.rcParams.update({

        "figure.figsize": (FIG_WIDTH, FIG_HEIGHT),

        "figure.dpi": DPI,

        "savefig.dpi": DPI,

        "font.family": "DejaVu Sans",

        "font.size": LABEL_SIZE,

        "axes.titlesize": TITLE_SIZE,

        "axes.labelsize": LABEL_SIZE,

        "axes.linewidth": 1.2,

        "xtick.labelsize": TICK_SIZE,

        "ytick.labelsize": TICK_SIZE,

        "legend.fontsize": LEGEND_SIZE,

        "grid.alpha": GRID_ALPHA,

        "lines.linewidth": LINE_WIDTH,

        "savefig.bbox": "tight"

    })


def finish_plot():

    """
    Apply common finishing touches before saving.
    """

    plt.grid(True, linestyle="--", alpha=GRID_ALPHA)

    plt.tight_layout()