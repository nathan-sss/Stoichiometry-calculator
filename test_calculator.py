"""Unit tests for the calculator module.

Run with:  pytest test_calculator.py -v
"""

import pytest
from calculator import (
    EndMember, Dopant, Reagent, calculate, build_composition_coeffs,
    format_composition, normalize_dict,
)
from data import ATOMIC_WEIGHTS


# ---------- Test fixtures ----------

def nbt_endmember():
    return EndMember(name="NBT", fraction=1.0,
                     A={"Na": 0.5, "Bi": 0.5}, B={"Ti": 1.0})

def bt_endmember(fraction=1.0):
    return EndMember(name="BT", fraction=fraction,
                     A={"Ba": 1.0}, B={"Ti": 1.0})

def default_reagents():
    return {
        "Na": Reagent("Na2CO3", "Na", 2, 105.99, 99.9, ""),
        "Bi": Reagent("Bi2O3", "Bi", 2, 465.96, 99.9, ""),
        "Ti": Reagent("TiO2", "Ti", 1, 79.87, 99.9, ""),
        "Ba": Reagent("BaCO3", "Ba", 1, 197.34, 99.9, ""),
        "Ca": Reagent("CaCO3", "Ca", 1, 100.09, 99.9, ""),
        "Nb": Reagent("Nb2O5", "Nb", 2, 265.81, 99.9, ""),
    }


# ---------- Composition coefficient tests ----------

def test_nbt_coefficients():
    """Pure NBT should have A: Na=0.5, Bi=0.5; B: Ti=1.0"""
    a, b = build_composition_coeffs([nbt_endmember()], None)
    assert a == {"Na": 0.5, "Bi": 0.5}
    assert b == {"Ti": 1.0}


def test_nbt_bt_solid_solution():
    """0.94 NBT - 0.06 BT solid solution."""
    em = [
        EndMember("NBT", 0.94, {"Na": 0.5, "Bi": 0.5}, {"Ti": 1.0}),
        EndMember("BT", 0.06, {"Ba": 1.0}, {"Ti": 1.0}),
    ]
    a, b = build_composition_coeffs(em, None)
    assert a["Na"] == pytest.approx(0.94 * 0.5)
    assert a["Bi"] == pytest.approx(0.94 * 0.5)
    assert a["Ba"] == pytest.approx(0.06)
    assert b["Ti"] == pytest.approx(0.94 + 0.06)


def test_dopant_a_site():
    """Adding 1.5 mol% Ca dopant on A-site of NBT."""
    dop = Dopant(enabled=True, site="A", cation="Ca", level=0.015)
    a, b = build_composition_coeffs([nbt_endmember()], dop)
    # Host A-site cations scaled by (1 - 0.015)
    assert a["Na"] == pytest.approx(0.5 * 0.985)
    assert a["Bi"] == pytest.approx(0.5 * 0.985)
    # Dopant added at 0.015
    assert a["Ca"] == pytest.approx(0.015)
    # B-site untouched
    assert b == {"Ti": 1.0}
    # Total A-site still sums to 1
    assert sum(a.values()) == pytest.approx(1.0)


def test_dopant_b_site():
    """Adding 0.75 mol% Nb dopant on B-site of NBT."""
    dop = Dopant(enabled=True, site="B", cation="Nb", level=0.0075)
    a, b = build_composition_coeffs([nbt_endmember()], dop)
    assert b["Ti"] == pytest.approx(1.0 * (1 - 0.0075))
    assert b["Nb"] == pytest.approx(0.0075)
    assert sum(b.values()) == pytest.approx(1.0)


def test_dopant_disabled():
    """Disabled dopant should not affect composition."""
    dop = Dopant(enabled=False, site="A", cation="Ca", level=0.015)
    a_no_dop, b_no_dop = build_composition_coeffs([nbt_endmember()], None)
    a_dop, b_dop = build_composition_coeffs([nbt_endmember()], dop)
    assert a_no_dop == a_dop
    assert b_no_dop == b_dop


# ---------- Molar mass and moles ----------

def test_nbt_molar_mass():
    """Compute MW of pure NBT (Na0.5 Bi0.5 Ti O3)."""
    result = calculate(
        batch_size_g=50.0,
        end_members=[nbt_endmember()],
        dopant=None,
        reagents_by_element=default_reagents(),
        excess_percent={},
        apply_purity=False,
    )
    expected_mw = (0.5 * ATOMIC_WEIGHTS["Na"] +
                   0.5 * ATOMIC_WEIGHTS["Bi"] +
                   1.0 * ATOMIC_WEIGHTS["Ti"] +
                   3.0 * ATOMIC_WEIGHTS["O"])
    assert result.total_mw == pytest.approx(expected_mw)
    assert result.moles == pytest.approx(50.0 / expected_mw)


def test_bt_calculation():
    """BaTiO3 sanity check: 1 mol of BaTiO3 = 233.19 g/mol."""
    result = calculate(
        batch_size_g=233.19,
        end_members=[bt_endmember()],
        dopant=None,
        reagents_by_element=default_reagents(),
        excess_percent={},
        apply_purity=False,
    )
    expected_mw = (ATOMIC_WEIGHTS["Ba"] + ATOMIC_WEIGHTS["Ti"] + 3 * ATOMIC_WEIGHTS["O"])
    assert result.total_mw == pytest.approx(expected_mw, rel=1e-4)
    assert result.moles == pytest.approx(1.0, rel=1e-3)

    # 1 mol of BaCO3 needed = 197.34 g
    ba_row = next(r for r in result.rows if r.element == "Ba")
    assert ba_row.mass_target == pytest.approx(197.34, rel=1e-3)
    # 1 mol of TiO2 = 79.87 g
    ti_row = next(r for r in result.rows if r.element == "Ti")
    assert ti_row.mass_target == pytest.approx(79.87, rel=1e-3)


# ---------- Excess additions ----------

def test_excess_addition():
    """A +3% Na excess should give 3% more Na2CO3 mass."""
    result_no_excess = calculate(
        batch_size_g=50.0,
        end_members=[nbt_endmember()],
        dopant=None,
        reagents_by_element=default_reagents(),
        excess_percent={},
        apply_purity=False,
    )
    result_excess = calculate(
        batch_size_g=50.0,
        end_members=[nbt_endmember()],
        dopant=None,
        reagents_by_element=default_reagents(),
        excess_percent={"Na": 3.0},
        apply_purity=False,
    )
    na_no = next(r for r in result_no_excess.rows if r.element == "Na")
    na_with = next(r for r in result_excess.rows if r.element == "Na")
    # mass_target unchanged
    assert na_no.mass_target == pytest.approx(na_with.mass_target)
    # mass_with_excess is 3% higher
    assert na_with.mass_with_excess == pytest.approx(na_no.mass_with_excess * 1.03)


# ---------- Purity correction ----------

def test_purity_correction():
    """With purity = 99.9, mass_to_weigh = mass_with_excess * (100/99.9)."""
    reagents = default_reagents()
    result_with = calculate(
        batch_size_g=50.0,
        end_members=[nbt_endmember()],
        dopant=None,
        reagents_by_element=reagents,
        excess_percent={},
        apply_purity=True,
    )
    result_without = calculate(
        batch_size_g=50.0,
        end_members=[nbt_endmember()],
        dopant=None,
        reagents_by_element=reagents,
        excess_percent={},
        apply_purity=False,
    )
    for r_with, r_without in zip(result_with.rows, result_without.rows):
        # Without purity, mass_to_weigh equals mass_with_excess
        assert r_without.mass_to_weigh == pytest.approx(r_without.mass_with_excess)
        # With purity, mass_to_weigh is slightly higher
        expected = r_with.mass_with_excess * (100.0 / 99.9)
        assert r_with.mass_to_weigh == pytest.approx(expected)


def test_purity_100_percent_no_correction():
    """At 100% purity, purity correction should not change the mass."""
    reagents = {"Na": Reagent("Na2CO3", "Na", 2, 105.99, 100.0, ""),
                "Bi": Reagent("Bi2O3", "Bi", 2, 465.96, 100.0, ""),
                "Ti": Reagent("TiO2", "Ti", 1, 79.87, 100.0, "")}
    result = calculate(
        batch_size_g=50.0,
        end_members=[nbt_endmember()],
        dopant=None,
        reagents_by_element=reagents,
        excess_percent={},
        apply_purity=True,
    )
    for r in result.rows:
        assert r.mass_to_weigh == pytest.approx(r.mass_with_excess)


# ---------- Helpers ----------

def test_normalize_dict():
    d = {"Na": 0.4, "Bi": 0.4}
    n = normalize_dict(d)
    assert n["Na"] == pytest.approx(0.5)
    assert n["Bi"] == pytest.approx(0.5)
    assert sum(n.values()) == pytest.approx(1.0)


def test_normalize_empty():
    assert normalize_dict({}) == {}


def test_format_composition():
    em = [nbt_endmember()]
    s = format_composition(em, None)
    assert "Na" in s and "Bi" in s and "Ti" in s and "O3" in s


def test_format_with_dopant():
    em = [nbt_endmember()]
    dop = Dopant(enabled=True, site="A", cation="Ca", level=0.015)
    s = format_composition(em, dop)
    assert "Ca" in s
    assert "A-site" in s


# ---------- Round-trip serialization ----------

def test_endmember_roundtrip():
    em = nbt_endmember()
    em2 = EndMember.from_dict(em.to_dict())
    assert em2.name == em.name
    assert em2.fraction == em.fraction
    assert em2.A == em.A
    assert em2.B == em.B


def test_dopant_roundtrip():
    d = Dopant(enabled=True, site="B", cation="Nb", level=0.0075)
    d2 = Dopant.from_dict(d.to_dict())
    assert d2 == d


def test_reagent_roundtrip():
    r = Reagent("Na2CO3", "Na", 2, 105.99, 99.9, "Sodium carbonate")
    r2 = Reagent.from_dict(r.to_dict())
    assert r2 == r


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
