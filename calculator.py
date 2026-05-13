"""Pure calculation engine for perovskite stoichiometry.

No Streamlit dependency here — this is the testable math core.
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Dict, List, Optional
from copy import deepcopy

from data import ATOMIC_WEIGHTS


@dataclass
class EndMember:
    """A single end-member in a (possibly composite) perovskite composition."""
    name: str
    fraction: float
    A: Dict[str, float]   # A-site cations → coefficients (should sum to 1)
    B: Dict[str, float]   # B-site cations → coefficients (should sum to 1)

    def to_dict(self) -> dict:
        return {"name": self.name, "fraction": self.fraction, "A": dict(self.A), "B": dict(self.B)}

    @staticmethod
    def from_dict(d: dict) -> "EndMember":
        return EndMember(name=d["name"], fraction=float(d["fraction"]),
                         A={k: float(v) for k, v in d.get("A", {}).items()},
                         B={k: float(v) for k, v in d.get("B", {}).items()})


@dataclass
class Dopant:
    enabled: bool = False
    site: str = "A"       # "A" or "B"
    cation: str = "Ca"
    level: float = 0.015  # mol fraction (0..1)

    def to_dict(self) -> dict:
        return {"enabled": self.enabled, "site": self.site, "cation": self.cation, "level": self.level}

    @staticmethod
    def from_dict(d: dict) -> "Dopant":
        return Dopant(enabled=bool(d.get("enabled", False)),
                      site=d.get("site", "A"),
                      cation=d.get("cation", "Ca"),
                      level=float(d.get("level", 0.0)))


@dataclass
class Reagent:
    name: str
    element: str
    atoms: int      # atoms of `element` per formula unit of reagent
    mw: float       # molecular weight (g/mol)
    purity: float   # percent (e.g. 99.9)
    notes: str = ""

    def to_dict(self) -> dict:
        return {"name": self.name, "element": self.element, "atoms": self.atoms,
                "mw": self.mw, "purity": self.purity, "notes": self.notes}

    @staticmethod
    def from_dict(d: dict) -> "Reagent":
        return Reagent(name=d["name"], element=d["element"], atoms=int(d["atoms"]),
                       mw=float(d["mw"]), purity=float(d["purity"]), notes=d.get("notes", ""))


@dataclass
class CalculationRow:
    element: str
    on_a_site: bool
    on_b_site: bool
    atomic_weight: float
    coeff: float                # element coefficient in target f.u.
    effective_coeff: float      # coeff * (1 + excess/100)
    excess_pct: float
    reagent: Optional[Reagent]
    mass_target: float          # grams of reagent for target stoichiometry
    mass_with_excess: float     # grams accounting for excess additions
    mass_to_weigh: float        # final mass to weigh on the balance (with purity correction if enabled)


@dataclass
class CalculationResult:
    rows: List[CalculationRow]
    a_coeffs: Dict[str, float]
    b_coeffs: Dict[str, float]
    oxygen_coeff: float
    oxygen_mass: float
    total_mw: float
    moles: float
    total_mass_with_excess: float
    total_mass_to_weigh: float


def build_composition_coeffs(
    end_members: List[EndMember],
    dopant: Optional[Dopant] = None,
) -> tuple[Dict[str, float], Dict[str, float]]:
    """Combine end-members (weighted by mole fraction) and apply dopant substitution.

    Returns (A_site_coeffs, B_site_coeffs) keyed by element symbol.
    """
    A: Dict[str, float] = {}
    B: Dict[str, float] = {}

    for em in end_members:
        f = em.fraction or 0.0
        for el, c in em.A.items():
            A[el] = A.get(el, 0.0) + f * c
        for el, c in em.B.items():
            B[el] = B.get(el, 0.0) + f * c

    if dopant is not None and dopant.enabled and dopant.cation and dopant.level > 0:
        x = dopant.level
        target = A if dopant.site == "A" else B
        # Scale host cations down by (1 - x)
        for el in list(target.keys()):
            target[el] *= (1.0 - x)
        # Add dopant
        target[dopant.cation] = target.get(dopant.cation, 0.0) + x

    return A, B


def calculate(
    batch_size_g: float,
    end_members: List[EndMember],
    dopant: Optional[Dopant],
    reagents_by_element: Dict[str, Reagent],
    excess_percent: Dict[str, float],
    apply_purity: bool = True,
) -> CalculationResult:
    """Run the full stoichiometric calculation.

    Parameters
    ----------
    batch_size_g : target mass of final ceramic product (after firing) in grams.
    end_members  : list of EndMember objects (their fractions should sum to 1).
    dopant       : optional Dopant; if enabled, substitutes on its chosen site.
    reagents_by_element : {element_symbol: Reagent} — which reagent supplies each element.
    excess_percent : {element_symbol: percent} — extra moles of that element to compensate volatilization.
    apply_purity : if True, divide reagent mass by (purity/100) to get the mass to weigh.
    """
    a_coeffs, b_coeffs = build_composition_coeffs(end_members, dopant)

    # Total moles of formula units = batch_size / total_MW
    # Oxygen coefficient: ABO3 perovskite → 3 per formula unit, scaled by sum of end-member fractions
    em_sum = sum(em.fraction for em in end_members)
    oxygen_coeff = 3.0 * em_sum

    elements = sorted(set(list(a_coeffs.keys()) + list(b_coeffs.keys())))
    elements = [el for el in elements if (a_coeffs.get(el, 0) + b_coeffs.get(el, 0)) > 1e-12]

    # Molar mass of target f.u.
    total_metal_mw = 0.0
    for el in elements:
        coeff = a_coeffs.get(el, 0.0) + b_coeffs.get(el, 0.0)
        total_metal_mw += coeff * ATOMIC_WEIGHTS.get(el, 0.0)
    oxygen_mass = oxygen_coeff * ATOMIC_WEIGHTS["O"]
    total_mw = total_metal_mw + oxygen_mass

    moles = batch_size_g / total_mw if total_mw > 0 else 0.0

    rows: List[CalculationRow] = []
    for el in elements:
        coeff = a_coeffs.get(el, 0.0) + b_coeffs.get(el, 0.0)
        excess = excess_percent.get(el, 0.0)
        effective = coeff * (1.0 + excess / 100.0)
        reagent = reagents_by_element.get(el)
        if reagent is not None and reagent.atoms > 0:
            mass_target = moles * coeff * reagent.mw / reagent.atoms
            mass_with_excess = moles * effective * reagent.mw / reagent.atoms
            if apply_purity and reagent.purity > 0:
                mass_to_weigh = mass_with_excess * (100.0 / reagent.purity)
            else:
                mass_to_weigh = mass_with_excess
        else:
            mass_target = 0.0
            mass_with_excess = 0.0
            mass_to_weigh = 0.0

        rows.append(CalculationRow(
            element=el,
            on_a_site=a_coeffs.get(el, 0.0) > 0,
            on_b_site=b_coeffs.get(el, 0.0) > 0,
            atomic_weight=ATOMIC_WEIGHTS.get(el, 0.0),
            coeff=coeff,
            effective_coeff=effective,
            excess_pct=excess,
            reagent=reagent,
            mass_target=mass_target,
            mass_with_excess=mass_with_excess,
            mass_to_weigh=mass_to_weigh,
        ))

    total_mass_with_excess = sum(r.mass_with_excess for r in rows)
    total_mass_to_weigh = sum(r.mass_to_weigh for r in rows)

    return CalculationResult(
        rows=rows,
        a_coeffs=a_coeffs,
        b_coeffs=b_coeffs,
        oxygen_coeff=oxygen_coeff,
        oxygen_mass=oxygen_mass,
        total_mw=total_mw,
        moles=moles,
        total_mass_with_excess=total_mass_with_excess,
        total_mass_to_weigh=total_mass_to_weigh,
    )


def normalize_dict(d: Dict[str, float]) -> Dict[str, float]:
    """Return a copy of d with values rescaled to sum to 1.0. No-op if sum is 0."""
    s = sum(d.values())
    if s <= 0:
        return dict(d)
    return {k: v / s for k, v in d.items()}


def format_composition(end_members: List[EndMember], dopant: Optional[Dopant]) -> str:
    """Human-readable composition string."""
    parts = []
    for em in end_members:
        if em.fraction <= 0:
            continue
        a_str = "".join(f"{el}{'' if c == 1 else c}" for el, c in em.A.items())
        b_str = "".join(f"{el}{'' if c == 1 else c}" for el, c in em.B.items())
        prefix = "" if em.fraction == 1.0 else f"{em.fraction:.4g} "
        parts.append(f"{prefix}({a_str})({b_str})O3")
    s = " - ".join(parts)
    if dopant and dopant.enabled and dopant.cation and dopant.level > 0:
        s += f"  +{dopant.level:.4g} {dopant.cation} on {dopant.site}-site"
    return s
