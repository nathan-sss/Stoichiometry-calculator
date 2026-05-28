"""Unit tests for the calculator module.

Run with:  pytest test_calculator.py -v
"""

import pytest
from calculator import (
    EndMember, DopantEntry, Reagent, calculate, build_composition_coeffs,
    format_composition, normalize_dict, dopant_entry_from_legacy_dict,
    UNIT_MOL_PCT, UNIT_WT_PCT,
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
        "La": Reagent("La2O3", "La", 2, 325.81, 99.9, ""),
    }


# ---------- Composition coefficient tests ----------

def test_nbt_coefficients():
    """Pure NBT should have A: Na=0.5, Bi=0.5; B: Ti=1.0"""
    a, b = build_composition_coeffs([nbt_endmember()])
    assert a == {"Na": 0.5, "Bi": 0.5}
    assert b == {"Ti": 1.0}


def test_nbt_bt_solid_solution():
    """0.94 NBT - 0.06 BT solid solution."""
    em = [
        EndMember("NBT", 0.94, {"Na": 0.5, "Bi": 0.5}, {"Ti": 1.0}),
        EndMember("BT", 0.06, {"Ba": 1.0}, {"Ti": 1.0}),
    ]
    a, b = build_composition_coeffs(em)
    assert a["Na"] == pytest.approx(0.94 * 0.5)
    assert a["Bi"] == pytest.approx(0.94 * 0.5)
    assert a["Ba"] == pytest.approx(0.06)
    assert b["Ti"] == pytest.approx(0.94 + 0.06)


# ---------- Additive dopant model ----------

def test_dopant_does_not_change_base_masses():
    """A dopant must NOT alter the base reagent masses. The fix the user
    requested: adding 1 mol% Bi excess used to shrink Na2CO3; now it shouldn't."""
    base = calculate(
        batch_size_g=50.0,
        end_members=[nbt_endmember()],
        dopants=[],
        reagents_by_element=default_reagents(),
        excess_percent={},
        apply_purity=False,
    )
    with_dop = calculate(
        batch_size_g=50.0,
        end_members=[nbt_endmember()],
        dopants=[DopantEntry(cation="Bi", amount=1.0, unit=UNIT_MOL_PCT)],
        reagents_by_element=default_reagents(),
        excess_percent={},
        apply_purity=False,
    )
    base_by_el = {r.element: r for r in base.rows}
    dop_by_el = {r.element: r for r in with_dop.rows}
    for el in base_by_el:
        assert dop_by_el[el].mass_target == pytest.approx(base_by_el[el].mass_target), (
            f"Adding a dopant changed the base mass for {el} — that's the bug."
        )
        assert dop_by_el[el].mass_with_excess == pytest.approx(base_by_el[el].mass_with_excess)


def test_mol_pct_dopant_extra_mass():
    """A 1 mol% Bi dopant should add (1% × moles × MW_Bi2O3 / 2) of Bi2O3."""
    result = calculate(
        batch_size_g=50.0,
        end_members=[nbt_endmember()],
        dopants=[DopantEntry(cation="Bi", amount=1.0, unit=UNIT_MOL_PCT)],
        reagents_by_element=default_reagents(),
        excess_percent={},
        apply_purity=False,
    )
    assert len(result.dopant_rows) == 1
    d = result.dopant_rows[0]
    expected = result.moles * 0.01 * 465.96 / 2
    assert d.mass == pytest.approx(expected)
    assert d.mass_to_weigh == pytest.approx(expected)  # purity off


def test_wt_pct_dopant_extra_mass():
    """A 2 wt% La dopant should add 0.02 × batch_size grams of La2O3."""
    result = calculate(
        batch_size_g=50.0,
        end_members=[nbt_endmember()],
        dopants=[DopantEntry(cation="La", amount=2.0, unit=UNIT_WT_PCT)],
        reagents_by_element=default_reagents(),
        excess_percent={},
        apply_purity=False,
    )
    d = result.dopant_rows[0]
    assert d.mass == pytest.approx(50.0 * 0.02)


def test_multiple_dopants():
    """Multiple dopants accumulate, total_dopant_mass sums them."""
    result = calculate(
        batch_size_g=50.0,
        end_members=[nbt_endmember()],
        dopants=[
            DopantEntry(cation="Bi", amount=1.0, unit=UNIT_MOL_PCT),
            DopantEntry(cation="La", amount=2.0, unit=UNIT_WT_PCT),
        ],
        reagents_by_element=default_reagents(),
        excess_percent={},
        apply_purity=False,
    )
    assert len(result.dopant_rows) == 2
    expected = sum(d.mass for d in result.dopant_rows)
    assert result.total_dopant_mass == pytest.approx(expected)
    assert result.grand_total_to_weigh == pytest.approx(
        result.total_mass_to_weigh + result.total_dopant_to_weigh
    )


def test_zero_amount_dopant_skipped():
    """A dopant with amount=0 should produce no row."""
    result = calculate(
        batch_size_g=50.0,
        end_members=[nbt_endmember()],
        dopants=[DopantEntry(cation="Ca", amount=0.0, unit=UNIT_MOL_PCT)],
        reagents_by_element=default_reagents(),
        excess_percent={},
        apply_purity=False,
    )
    assert result.dopant_rows == []


def test_dopant_purity_correction():
    """With purity 99.9, dopant mass_to_weigh = mass × 100/99.9."""
    result = calculate(
        batch_size_g=50.0,
        end_members=[nbt_endmember()],
        dopants=[DopantEntry(cation="La", amount=1.0, unit=UNIT_MOL_PCT)],
        reagents_by_element=default_reagents(),
        excess_percent={},
        apply_purity=True,
    )
    d = result.dopant_rows[0]
    assert d.mass_to_weigh == pytest.approx(d.mass * 100.0 / 99.9)


# ---------- Molar mass and moles ----------

def test_nbt_molar_mass():
    """Compute MW of pure NBT (Na0.5 Bi0.5 Ti O3)."""
    result = calculate(
        batch_size_g=50.0,
        end_members=[nbt_endmember()],
        dopants=[],
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
        dopants=[],
        reagents_by_element=default_reagents(),
        excess_percent={},
        apply_purity=False,
    )
    expected_mw = (ATOMIC_WEIGHTS["Ba"] + ATOMIC_WEIGHTS["Ti"] + 3 * ATOMIC_WEIGHTS["O"])
    assert result.total_mw == pytest.approx(expected_mw, rel=1e-4)
    assert result.moles == pytest.approx(1.0, rel=1e-3)

    ba_row = next(r for r in result.rows if r.element == "Ba")
    assert ba_row.mass_target == pytest.approx(197.34, rel=1e-3)
    ti_row = next(r for r in result.rows if r.element == "Ti")
    assert ti_row.mass_target == pytest.approx(79.87, rel=1e-3)


# ---------- Excess additions ----------

def test_excess_addition():
    """A +3% Na excess should give 3% more Na2CO3 mass."""
    result_no_excess = calculate(
        batch_size_g=50.0,
        end_members=[nbt_endmember()],
        dopants=[],
        reagents_by_element=default_reagents(),
        excess_percent={},
        apply_purity=False,
    )
    result_excess = calculate(
        batch_size_g=50.0,
        end_members=[nbt_endmember()],
        dopants=[],
        reagents_by_element=default_reagents(),
        excess_percent={"Na": 3.0},
        apply_purity=False,
    )
    na_no = next(r for r in result_no_excess.rows if r.element == "Na")
    na_with = next(r for r in result_excess.rows if r.element == "Na")
    assert na_no.mass_target == pytest.approx(na_with.mass_target)
    assert na_with.mass_with_excess == pytest.approx(na_no.mass_with_excess * 1.03)


# ---------- Purity correction ----------

def test_purity_correction():
    """With purity = 99.9, mass_to_weigh = mass_with_excess * (100/99.9)."""
    reagents = default_reagents()
    result_with = calculate(
        batch_size_g=50.0,
        end_members=[nbt_endmember()],
        dopants=[],
        reagents_by_element=reagents,
        excess_percent={},
        apply_purity=True,
    )
    result_without = calculate(
        batch_size_g=50.0,
        end_members=[nbt_endmember()],
        dopants=[],
        reagents_by_element=reagents,
        excess_percent={},
        apply_purity=False,
    )
    for r_with, r_without in zip(result_with.rows, result_without.rows):
        assert r_without.mass_to_weigh == pytest.approx(r_without.mass_with_excess)
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
        dopants=[],
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
    dops = [DopantEntry(cation="Ca", amount=1.5, unit=UNIT_MOL_PCT)]
    s = format_composition(em, dops)
    assert "Ca" in s
    assert "mol%" in s


def test_format_with_wt_pct_dopant():
    em = [nbt_endmember()]
    dops = [DopantEntry(cation="La", amount=2.0, unit=UNIT_WT_PCT)]
    s = format_composition(em, dops)
    assert "La" in s
    assert "wt%" in s


# ---------- Round-trip serialization ----------

def test_endmember_roundtrip():
    em = nbt_endmember()
    em2 = EndMember.from_dict(em.to_dict())
    assert em2.name == em.name
    assert em2.fraction == em.fraction
    assert em2.A == em.A
    assert em2.B == em.B


def test_dopant_entry_roundtrip():
    d = DopantEntry(cation="Nb", amount=0.75, unit=UNIT_MOL_PCT)
    d2 = DopantEntry.from_dict(d.to_dict())
    assert d2 == d


def test_reagent_roundtrip():
    r = Reagent("Na2CO3", "Na", 2, 105.99, 99.9, "Sodium carbonate")
    r2 = Reagent.from_dict(r.to_dict())
    assert r2 == r


# ---------- Legacy migration ----------

def test_legacy_dopant_migration_enabled():
    """The old {enabled, site, cation, level} dict should convert to a
    mol%-additive DopantEntry (legacy.level was a mol fraction, new amount is %)."""
    legacy = {"enabled": True, "site": "A", "cation": "Ca", "level": 0.015}
    entry = dopant_entry_from_legacy_dict(legacy)
    assert entry is not None
    assert entry.cation == "Ca"
    assert entry.amount == pytest.approx(1.5)
    assert entry.unit == UNIT_MOL_PCT


def test_legacy_dopant_migration_disabled():
    """Disabled or zero-level legacy dopants drop to None."""
    assert dopant_entry_from_legacy_dict({"enabled": False, "cation": "Ca", "level": 0.015}) is None
    assert dopant_entry_from_legacy_dict({"enabled": True, "cation": "Ca", "level": 0.0}) is None
    assert dopant_entry_from_legacy_dict({}) is None


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
