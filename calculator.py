"""Pure calculation engine for perovskite stoichiometry.

No Streamlit dependency here — this is the testable math core.
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Dict, List, Optional
from copy import deepcopy

from data import ATOMIC_WEIGHTS


# Dopant unit identifiers used in DopantEntry.unit
UNIT_MOL_PCT = "mol_pct"
UNIT_WT_PCT = "wt_pct"


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
class DopantEntry:
    """An additive dopant — weighed in on top of the base batch, not substituted.

    `amount` is the user-supplied number; `unit` selects how to interpret it:
      - "mol_pct": amount mol% of the base composition's formula-unit moles
                   (so 1.0 means +0.01 × moles of dopant cation added).
      - "wt_pct":  amount wt% of the base batch_size_g of dopant reagent
                   (so 1.0 means dopant reagent mass = 0.01 × batch_size_g).
    """
    cation: str
    amount: float = 1.0
    unit: str = UNIT_MOL_PCT

    def to_dict(self) -> dict:
        return {"cation": self.cation, "amount": self.amount, "unit": self.unit}

    @staticmethod
    def from_dict(d: dict) -> "DopantEntry":
        return DopantEntry(
            cation=d["cation"],
            amount=float(d.get("amount", 0.0)),
            unit=d.get("unit", UNIT_MOL_PCT),
        )


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
class DopantRow:
    """Result row for a single additive dopant."""
    cation: str
    amount: float
    unit: str
    atomic_weight: float
    reagent: Optional[Reagent]
    mass: float            # reagent mass before purity correction
    mass_to_weigh: float   # reagent mass after purity correction (if applied)


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
    dopant_rows: List[DopantRow] = field(default_factory=list)
    total_dopant_mass: float = 0.0          # sum of dopant masses before purity
    total_dopant_to_weigh: float = 0.0      # sum of dopant masses with purity
    grand_total_to_weigh: float = 0.0       # base + dopants, what the user actually puts on the scale


def build_composition_coeffs(
    end_members: List[EndMember],
) -> tuple[Dict[str, float], Dict[str, float]]:
    """Combine end-members (weighted by mole fraction) into A/B-site coefficients.

    Dopants are additive (weighed in separately) and do NOT modify these
    coefficients. Returns (A_site_coeffs, B_site_coeffs) keyed by element symbol.
    """
    A: Dict[str, float] = {}
    B: Dict[str, float] = {}

    for em in end_members:
        f = em.fraction or 0.0
        for el, c in em.A.items():
            A[el] = A.get(el, 0.0) + f * c
        for el, c in em.B.items():
            B[el] = B.get(el, 0.0) + f * c

    return A, B


def calculate(
    batch_size_g: float,
    end_members: List[EndMember],
    dopants: List[DopantEntry],
    reagents_by_element: Dict[str, Reagent],
    excess_percent: Dict[str, float],
    apply_purity: bool = True,
) -> CalculationResult:
    """Run the full stoichiometric calculation.

    Parameters
    ----------
    batch_size_g : target mass of the base ceramic product (after firing), in grams.
                   Dopants are added ON TOP of this — they don't reduce the host
                   reagent masses.
    end_members  : list of EndMember objects (their fractions should sum to 1).
    dopants      : list of DopantEntry; each contributes additional reagent mass.
    reagents_by_element : {element_symbol: Reagent} — which reagent supplies each element.
    excess_percent : {element_symbol: percent} — extra moles of that element to
                     compensate volatilization. Applies to base composition only.
    apply_purity : if True, divide reagent mass by (purity/100) to get the mass to weigh.
    """
    a_coeffs, b_coeffs = build_composition_coeffs(end_members)

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

    # Dopants are added on top of the base batch — they do not modify rows above.
    dopant_rows: List[DopantRow] = []
    for dop in dopants or []:
        if not dop.cation or dop.amount <= 0:
            continue
        reagent = reagents_by_element.get(dop.cation)
        if dop.unit == UNIT_WT_PCT:
            # amount wt% of the BASE batch mass, expressed as the dopant reagent.
            mass = batch_size_g * (dop.amount / 100.0)
        else:
            # Default: mol% of base formula-unit moles.
            extra_moles_of_cation = moles * (dop.amount / 100.0)
            if reagent is not None and reagent.atoms > 0:
                mass = extra_moles_of_cation * reagent.mw / reagent.atoms
            else:
                mass = 0.0

        if apply_purity and reagent is not None and reagent.purity > 0:
            mass_to_weigh = mass * (100.0 / reagent.purity)
        else:
            mass_to_weigh = mass

        dopant_rows.append(DopantRow(
            cation=dop.cation,
            amount=dop.amount,
            unit=dop.unit,
            atomic_weight=ATOMIC_WEIGHTS.get(dop.cation, 0.0),
            reagent=reagent,
            mass=mass,
            mass_to_weigh=mass_to_weigh,
        ))

    total_dopant_mass = sum(d.mass for d in dopant_rows)
    total_dopant_to_weigh = sum(d.mass_to_weigh for d in dopant_rows)
    grand_total_to_weigh = total_mass_to_weigh + total_dopant_to_weigh

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
        dopant_rows=dopant_rows,
        total_dopant_mass=total_dopant_mass,
        total_dopant_to_weigh=total_dopant_to_weigh,
        grand_total_to_weigh=grand_total_to_weigh,
    )


def normalize_dict(d: Dict[str, float]) -> Dict[str, float]:
    """Return a copy of d with values rescaled to sum to 1.0. No-op if sum is 0."""
    s = sum(d.values())
    if s <= 0:
        return dict(d)
    return {k: v / s for k, v in d.items()}


def _format_dopant(dop: DopantEntry) -> str:
    unit = "wt%" if dop.unit == UNIT_WT_PCT else "mol%"
    return f"{dop.amount:.4g} {unit} {dop.cation}"


def format_composition(
    end_members: List[EndMember],
    dopants: Optional[List[DopantEntry]] = None,
) -> str:
    """Human-readable composition string. Dopants are listed as additive add-ons."""
    parts = []
    for em in end_members:
        if em.fraction <= 0:
            continue
        a_str = "".join(f"{el}{'' if c == 1 else c}" for el, c in em.A.items())
        b_str = "".join(f"{el}{'' if c == 1 else c}" for el, c in em.B.items())
        prefix = "" if em.fraction == 1.0 else f"{em.fraction:.4g} "
        parts.append(f"{prefix}({a_str})({b_str})O3")
    s = " - ".join(parts)
    active = [d for d in (dopants or []) if d.cation and d.amount > 0]
    if active:
        s += "  + " + ", ".join(_format_dopant(d) for d in active)
    return s


def dopant_entry_from_legacy_dict(d: dict) -> Optional[DopantEntry]:
    """Convert the legacy {enabled, site, cation, level} dopant dict to a
    single mol%-additive DopantEntry. Returns None for disabled / zero-level
    entries. Used when loading recipes saved with the old substitutive model."""
    if not d or not d.get("enabled"):
        return None
    level = float(d.get("level", 0.0) or 0.0)
    cation = d.get("cation", "")
    if level <= 0 or not cation:
        return None
    return DopantEntry(cation=cation, amount=level * 100.0, unit=UNIT_MOL_PCT)
