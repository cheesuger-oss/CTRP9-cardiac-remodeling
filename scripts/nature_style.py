#!/usr/bin/env python3
"""Shared Nature-style plotting configuration (backend: python only).

Implements the nature-figure quick-start:
- Arial/Helvetica sans-serif, base font size 7 pt
- editable text: svg.fonttype='none', pdf.fonttype=42
- spines off top/right, thin axes, no legend frame
- restrained palette (one neutral + one signal + one accent family)
- save_pub(): SVG + PDF + TIFF(600 dpi) + PNG preview(600 dpi)
- panel_label(): lowercase bold panel letters near top-left
"""
import os
os.environ.setdefault("MPLCONFIGDIR", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".mplcache"))
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
    "svg.fonttype": "none",       # editable text in SVG
    "pdf.fonttype": 42,           # editable TrueType text in PDF
    "font.size": 7,
    "axes.titlesize": 8,
    "axes.labelsize": 7.5,
    "xtick.labelsize": 6.5,
    "ytick.labelsize": 6.5,
    "legend.fontsize": 6.5,
    "axes.spines.right": False,
    "axes.spines.top": False,
    "axes.linewidth": 0.8,
    "legend.frameon": False,
    "figure.dpi": 100,
    "savefig.dpi": 600,
})

# Journal column widths (inches) for final-size figures
MM = 1.0 / 25.4
W_SINGLE = 89 * MM   # single column
W_DOUBLE = 183 * MM  # double column

# Restrained palette: neutral / signal / accent
C_CTL = "#4C72B0"   # control / normal (blue)
C_CASE = "#C44E52"  # disease (red)
C_ACC = "#55A868"   # accent (green)
C_PUR = "#8172B3"   # secondary (purple)
C_GRAY = "#7F8C8D"  # neutral gray
GRID = "#D9D9D9"

def panel_label(ax, text, x=0.015, y=0.97, size=8):
    """Lowercase bold panel letter near top-left of the axes."""
    ax.text(x, y, text, transform=ax.transAxes, fontsize=size, fontweight="bold",
            va="top", ha="left", color="black")

def save_pub(fig, filename, dpi=600):
    """Export SVG + PDF + TIFF + PNG preview with editable text."""
    fig.savefig(f"{filename}.svg", bbox_inches="tight")
    fig.savefig(f"{filename}.pdf", bbox_inches="tight")
    fig.savefig(f"{filename}.tiff", dpi=dpi, bbox_inches="tight")
    fig.savefig(f"{filename}.png", dpi=dpi, bbox_inches="tight")
    plt.close(fig)
