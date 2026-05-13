"""Static data: atomic weights, default reagents, perovskite presets, periodic table layout."""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Dict, List


ATOMIC_WEIGHTS: Dict[str, float] = {
    "H": 1.008, "He": 4.003,
    "Li": 6.94, "Be": 9.012, "B": 10.81, "C": 12.011, "N": 14.007, "O": 15.999, "F": 18.998, "Ne": 20.18,
    "Na": 22.99, "Mg": 24.305, "Al": 26.982, "Si": 28.085, "P": 30.974, "S": 32.06, "Cl": 35.45, "Ar": 39.948,
    "K": 39.098, "Ca": 40.078, "Sc": 44.956, "Ti": 47.867, "V": 50.942, "Cr": 51.996, "Mn": 54.938,
    "Fe": 55.845, "Co": 58.933, "Ni": 58.693, "Cu": 63.546, "Zn": 65.38, "Ga": 69.723, "Ge": 72.63,
    "As": 74.922, "Se": 78.971, "Br": 79.904, "Kr": 83.798,
    "Rb": 85.468, "Sr": 87.62, "Y": 88.906, "Zr": 91.224, "Nb": 92.906, "Mo": 95.95,
    "Ru": 101.07, "Rh": 102.91, "Pd": 106.42, "Ag": 107.87, "Cd": 112.41, "In": 114.82, "Sn": 118.71,
    "Sb": 121.76, "Te": 127.6, "I": 126.9, "Xe": 131.29,
    "Cs": 132.91, "Ba": 137.327, "La": 138.91, "Ce": 140.12, "Pr": 140.91, "Nd": 144.24,
    "Sm": 150.36, "Eu": 151.96, "Gd": 157.25, "Tb": 158.93, "Dy": 162.5, "Ho": 164.93, "Er": 167.26,
    "Tm": 168.93, "Yb": 173.05, "Lu": 174.97,
    "Hf": 178.49, "Ta": 180.95, "W": 183.84, "Re": 186.21, "Os": 190.23, "Ir": 192.22, "Pt": 195.08,
    "Au": 196.97, "Hg": 200.59, "Tl": 204.38, "Pb": 207.2, "Bi": 208.98,
}


# Typical site occupancy in perovskites (used for color-coding the periodic table picker)
A_SITE_CATIONS: List[str] = [
    "Na", "K", "Li", "Rb", "Cs", "Ca", "Sr", "Ba", "Pb", "Bi",
    "La", "Ce", "Pr", "Nd", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Y", "Ag",
]

B_SITE_CATIONS: List[str] = [
    "Ti", "Zr", "Hf", "Nb", "Ta", "V", "Fe", "Mn", "Mg", "Zn", "Cu", "Al",
    "Sn", "W", "Mo", "Co", "Ni", "Sc", "Cr", "Ga", "In", "Sb", "Ru", "Rh", "Ir", "Pt",
]


# Periodic table layout: (symbol, row, col). Lanthanides placed at row 8.
PERIODIC_TABLE_LAYOUT: List[tuple] = [
    ("H", 1, 1), ("He", 1, 18),
    ("Li", 2, 1), ("Be", 2, 2), ("B", 2, 13), ("C", 2, 14), ("N", 2, 15), ("O", 2, 16), ("F", 2, 17), ("Ne", 2, 18),
    ("Na", 3, 1), ("Mg", 3, 2), ("Al", 3, 13), ("Si", 3, 14), ("P", 3, 15), ("S", 3, 16), ("Cl", 3, 17), ("Ar", 3, 18),
    ("K", 4, 1), ("Ca", 4, 2), ("Sc", 4, 3), ("Ti", 4, 4), ("V", 4, 5), ("Cr", 4, 6), ("Mn", 4, 7),
    ("Fe", 4, 8), ("Co", 4, 9), ("Ni", 4, 10), ("Cu", 4, 11), ("Zn", 4, 12),
    ("Ga", 4, 13), ("Ge", 4, 14), ("As", 4, 15), ("Se", 4, 16), ("Br", 4, 17), ("Kr", 4, 18),
    ("Rb", 5, 1), ("Sr", 5, 2), ("Y", 5, 3), ("Zr", 5, 4), ("Nb", 5, 5), ("Mo", 5, 6),
    ("Ru", 5, 8), ("Rh", 5, 9), ("Pd", 5, 10), ("Ag", 5, 11), ("Cd", 5, 12),
    ("In", 5, 13), ("Sn", 5, 14), ("Sb", 5, 15), ("Te", 5, 16), ("I", 5, 17), ("Xe", 5, 18),
    ("Cs", 6, 1), ("Ba", 6, 2),
    ("Hf", 6, 4), ("Ta", 6, 5), ("W", 6, 6), ("Re", 6, 7), ("Os", 6, 8), ("Ir", 6, 9), ("Pt", 6, 10),
    ("Au", 6, 11), ("Hg", 6, 12), ("Tl", 6, 13), ("Pb", 6, 14), ("Bi", 6, 15),
    # Lanthanide row (separate row 8)
    ("La", 8, 4), ("Ce", 8, 5), ("Pr", 8, 6), ("Nd", 8, 7), ("Sm", 8, 8), ("Eu", 8, 9), ("Gd", 8, 10),
    ("Tb", 8, 11), ("Dy", 8, 12), ("Ho", 8, 13), ("Er", 8, 14), ("Tm", 8, 15), ("Yb", 8, 16), ("Lu", 8, 17),
]


# Default reagent database — edit/extend freely
DEFAULT_REAGENTS: List[dict] = [
    {"name": "Na2CO3", "element": "Na", "atoms": 2, "mw": 105.99, "purity": 99.9, "notes": "Sodium carbonate"},
    {"name": "Bi2O3", "element": "Bi", "atoms": 2, "mw": 465.96, "purity": 99.9, "notes": "Bismuth oxide"},
    {"name": "BaCO3", "element": "Ba", "atoms": 1, "mw": 197.34, "purity": 99.9, "notes": "Barium carbonate"},
    {"name": "TiO2", "element": "Ti", "atoms": 1, "mw": 79.87, "purity": 99.9, "notes": "Titanium dioxide"},
    {"name": "CaCO3", "element": "Ca", "atoms": 1, "mw": 100.09, "purity": 99.9, "notes": "Calcium carbonate"},
    {"name": "Nb2O5", "element": "Nb", "atoms": 2, "mw": 265.81, "purity": 99.9, "notes": "Niobium pentoxide"},
    {"name": "K2CO3", "element": "K", "atoms": 2, "mw": 138.205, "purity": 99.9, "notes": "Potassium carbonate"},
    {"name": "MgCO3", "element": "Mg", "atoms": 1, "mw": 84.3149, "purity": 99.9, "notes": "Magnesium carbonate"},
    {"name": "SrCO3", "element": "Sr", "atoms": 1, "mw": 147.63, "purity": 99.9, "notes": "Strontium carbonate"},
    {"name": "Fe2O3", "element": "Fe", "atoms": 2, "mw": 159.69, "purity": 99.9, "notes": "Iron(III) oxide"},
    {"name": "ZrO2", "element": "Zr", "atoms": 1, "mw": 123.22, "purity": 99.9, "notes": "Zirconium dioxide"},
    {"name": "PbO", "element": "Pb", "atoms": 1, "mw": 223.2, "purity": 99.9, "notes": "Lead(II) oxide"},
    {"name": "MnO2", "element": "Mn", "atoms": 1, "mw": 86.94, "purity": 99.9, "notes": "Manganese dioxide"},
    {"name": "La2O3", "element": "La", "atoms": 2, "mw": 325.81, "purity": 99.9, "notes": "Lanthanum oxide"},
    {"name": "Y2O3", "element": "Y", "atoms": 2, "mw": 225.81, "purity": 99.9, "notes": "Yttrium oxide"},
    {"name": "Ta2O5", "element": "Ta", "atoms": 2, "mw": 441.89, "purity": 99.9, "notes": "Tantalum pentoxide"},
    {"name": "CuO", "element": "Cu", "atoms": 1, "mw": 79.55, "purity": 99.9, "notes": "Copper(II) oxide"},
    {"name": "ZnO", "element": "Zn", "atoms": 1, "mw": 81.38, "purity": 99.9, "notes": "Zinc oxide"},
]


# Common perovskite presets
PRESETS: List[dict] = [
    {
        "id": "nbt", "name": "NBT", "full": "Na0.5Bi0.5TiO3",
        "end_members": [{"name": "NBT", "fraction": 1.0, "A": {"Na": 0.5, "Bi": 0.5}, "B": {"Ti": 1.0}}],
    },
    {
        "id": "bt", "name": "BT", "full": "BaTiO3",
        "end_members": [{"name": "BaTiO3", "fraction": 1.0, "A": {"Ba": 1.0}, "B": {"Ti": 1.0}}],
    },
    {
        "id": "nbt-bt", "name": "NBT-BT", "full": "0.94 NBT - 0.06 BT",
        "end_members": [
            {"name": "NBT", "fraction": 0.94, "A": {"Na": 0.5, "Bi": 0.5}, "B": {"Ti": 1.0}},
            {"name": "BaTiO3", "fraction": 0.06, "A": {"Ba": 1.0}, "B": {"Ti": 1.0}},
        ],
    },
    {
        "id": "knn", "name": "KNN", "full": "(K0.5Na0.5)NbO3",
        "end_members": [{"name": "KNN", "fraction": 1.0, "A": {"K": 0.5, "Na": 0.5}, "B": {"Nb": 1.0}}],
    },
    {
        "id": "bfbt", "name": "BFBT", "full": "0.7 BiFeO3 - 0.3 BaTiO3",
        "end_members": [
            {"name": "BiFeO3", "fraction": 0.7, "A": {"Bi": 1.0}, "B": {"Fe": 1.0}},
            {"name": "BaTiO3", "fraction": 0.3, "A": {"Ba": 1.0}, "B": {"Ti": 1.0}},
        ],
    },
    {
        "id": "bf", "name": "BiFeO3", "full": "BiFeO3",
        "end_members": [{"name": "BiFeO3", "fraction": 1.0, "A": {"Bi": 1.0}, "B": {"Fe": 1.0}}],
    },
    {
        "id": "pzt", "name": "PZT", "full": "PbZr0.52Ti0.48O3",
        "end_members": [{"name": "PZT", "fraction": 1.0, "A": {"Pb": 1.0}, "B": {"Zr": 0.52, "Ti": 0.48}}],
    },
    {
        "id": "st", "name": "SrTiO3", "full": "SrTiO3",
        "end_members": [{"name": "SrTiO3", "fraction": 1.0, "A": {"Sr": 1.0}, "B": {"Ti": 1.0}}],
    },
]


NON_CATION_ELEMENTS = {"H", "He", "C", "N", "O", "F", "Ne", "Si", "P", "S", "Cl", "Ar",
                       "Ge", "As", "Se", "Br", "Kr", "Te", "I", "Xe"}
