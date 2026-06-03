"""Unit tests for the calculator module.

Run with:  pytest test_calculator.py -v
"""

import pytest
from calculator import (
    EndMember, DopantEntry, Reagent, calculate, build_composition_coeffs,
    format_composition, normalize_dict, dopant_entry_from_legacy_dict,
    UNIT_MOL_PCT, UNIT_WT_PCT, MODE_SUBSTITUTIONAL, MODE_ADDITIVE,
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
    """The old {enabled, site, cation, level} dict was substitutional, so it must
    convert to a SUBSTITUTIONAL DopantEntry (level mol fraction → amount in %)."""
    legacy = {"enabled": True, "site": "A", "cation": "Ca", "level": 0.015}
    entry = dopant_entry_from_legacy_dict(legacy)
    assert entry is not None
    assert entry.cation == "Ca"
    assert entry.amount == pytest.approx(1.5)
    assert entry.mode == MODE_SUBSTITUTIONAL
    assert entry.site == "A"


def test_legacy_dopant_migration_disabled():
    """Disabled or zero-level legacy dopants drop to None."""
    assert dopant_entry_from_legacy_dict({"enabled": False, "cation": "Ca", "level": 0.015}) is None
    assert dopant_entry_from_legacy_dict({"enabled": True, "cation": "Ca", "level": 0.0}) is None
    assert dopant_entry_from_legacy_dict({}) is None


# ---------- Rigorous invariant checks on the core calculation ----------
#
# These don't trust any hand-computed number: they assert the *physical*
# relationships the calculation must satisfy, so they catch a wrong formula
# even if every example number happened to look plausible.

def test_base_reagent_delivers_exact_element_moles():
    """THE core invariant. For every base reagent, the mass we tell the user to
    weigh must convert back into exactly `moles × coeff` moles of that element.
    If this holds for an arbitrary composite, the reagent-mass formula is correct."""
    em = [
        EndMember("NBT", 0.94, {"Na": 0.5, "Bi": 0.5}, {"Ti": 1.0}),
        EndMember("BT", 0.06, {"Ba": 1.0}, {"Ti": 1.0}),
    ]
    result = calculate(50.0, em, [], default_reagents(), {}, apply_purity=False)
    for r in result.rows:
        moles_required = result.moles * r.coeff
        moles_delivered = (r.mass_target / r.reagent.mw) * r.reagent.atoms
        assert moles_delivered == pytest.approx(moles_required), r.element


def test_excess_delivers_exact_extra_moles():
    """A +5% excess must deliver exactly 1.05× the element moles, no more no less."""
    result = calculate(50.0, [nbt_endmember()], [], default_reagents(),
                       {"Na": 5.0}, apply_purity=False)
    na = next(r for r in result.rows if r.element == "Na")
    moles_required = result.moles * na.coeff
    moles_delivered = (na.mass_with_excess / na.reagent.mw) * na.reagent.atoms
    assert moles_delivered == pytest.approx(moles_required * 1.05)


def test_purity_delivers_exact_element_after_impurity():
    """After weighing the purity-corrected mass, the *pure* fraction of it must
    still deliver exactly the required element moles."""
    result = calculate(50.0, [nbt_endmember()], [], default_reagents(),
                       {}, apply_purity=True)
    for r in result.rows:
        moles_required = result.moles * r.coeff
        pure_mass = r.mass_to_weigh * (r.reagent.purity / 100.0)
        moles_delivered = (pure_mass / r.reagent.mw) * r.reagent.atoms
        assert moles_delivered == pytest.approx(moles_required), r.element


def test_dopant_mol_pct_delivers_exact_cation_moles():
    """A 1 mol% dopant must add exactly 0.01 × (formula-unit moles) of the cation."""
    result = calculate(50.0, [nbt_endmember()],
                       [DopantEntry("La", 1.0, unit=UNIT_MOL_PCT)],
                       default_reagents(), {}, apply_purity=False)
    d = result.dopant_rows[0]
    moles_delivered = (d.mass / d.reagent.mw) * d.reagent.atoms
    assert moles_delivered == pytest.approx(result.moles * 0.01)


def test_dopant_wt_pct_is_exact_fraction_of_batch():
    """A 2 wt% dopant reagent must weigh exactly 2% of the batch mass."""
    result = calculate(50.0, [nbt_endmember()],
                       [DopantEntry("La", 2.0, unit=UNIT_WT_PCT)],
                       default_reagents(), {}, apply_purity=False)
    assert result.dopant_rows[0].mass == pytest.approx(0.02 * 50.0)


def test_composite_molar_mass_independent_value():
    """Independent hand-computed MW for 0.94 NBT – 0.06 BT (Na0.47 Bi0.47 Ba0.06
    Ti O3): 213.13 g/mol. Pins total_mw against an external number."""
    em = [
        EndMember("NBT", 0.94, {"Na": 0.5, "Bi": 0.5}, {"Ti": 1.0}),
        EndMember("BT", 0.06, {"Ba": 1.0}, {"Ti": 1.0}),
    ]
    result = calculate(100.0, em, [], default_reagents(), {}, apply_purity=False)
    assert result.total_mw == pytest.approx(213.13, abs=0.02)
    assert result.moles == pytest.approx(100.0 / result.total_mw)


def test_target_product_mass_reconstructs_from_moles():
    """moles × total_mw must return the batch size — i.e. the batch really is the
    fired-oxide product mass, oxygen included."""
    result = calculate(37.5, [nbt_endmember()], [], default_reagents(),
                       {}, apply_purity=False)
    assert result.moles * result.total_mw == pytest.approx(37.5)


def test_grand_total_is_base_plus_dopants():
    """Grand total to weigh = base reagents + dopant reagents, exactly."""
    result = calculate(50.0, [nbt_endmember()],
                       [DopantEntry("La", 1.5, unit=UNIT_MOL_PCT),
                        DopantEntry("Mn", 0.5, unit=UNIT_WT_PCT)],
                       default_reagents(), {}, apply_purity=True)
    assert result.grand_total_to_weigh == pytest.approx(
        result.total_mass_to_weigh + result.total_dopant_to_weigh
    )
    assert result.total_dopant_to_weigh == pytest.approx(
        sum(d.mass_to_weigh for d in result.dopant_rows)
    )


def test_element_on_both_sites_sums_coefficients():
    """If a cation appears on A in one end-member and B in another, its reagent
    mass must cover the TOTAL moles across both sites (e.g. Bi in NBT + BiFeO3)."""
    em = [
        EndMember("NBT", 0.5, {"Na": 0.5, "Bi": 0.5}, {"Ti": 1.0}),
        EndMember("BFO", 0.5, {"Bi": 1.0}, {"Fe": 1.0}),
    ]
    reagents = default_reagents()
    reagents["Fe"] = Reagent("Fe2O3", "Fe", 2, 159.69, 99.9, "")
    result = calculate(50.0, em, [], reagents, {}, apply_purity=False)
    bi = next(r for r in result.rows if r.element == "Bi")
    # Bi coeff = 0.5*0.5 (A in NBT) + 0.5*1.0 (A in BFO) = 0.75
    assert bi.coeff == pytest.approx(0.75)
    moles_delivered = (bi.mass_target / bi.reagent.mw) * bi.reagent.atoms
    assert moles_delivered == pytest.approx(result.moles * 0.75)


# ---------- Substitutional dopants (lattice doping) ----------

def _sub(cation, amount, site):
    return DopantEntry(cation=cation, amount=amount,
                       mode=MODE_SUBSTITUTIONAL, site=site)


def test_substitutional_scales_host_and_adds_dopant():
    """2 mol% La on the A-site of (Na0.5Bi0.5)TiO3 → Na,Bi each ×0.98, La=0.02,
    site sum preserved at 1."""
    a, b = build_composition_coeffs([nbt_endmember()], [_sub("La", 2.0, "A")])
    assert a["Na"] == pytest.approx(0.5 * 0.98)
    assert a["Bi"] == pytest.approx(0.5 * 0.98)
    assert a["La"] == pytest.approx(0.02)
    assert sum(a.values()) == pytest.approx(1.0)
    assert b == {"Ti": 1.0}


def test_substitutional_changes_mw_and_reduces_host_mass():
    """A substitutional dopant MUST change host reagent masses and the MW —
    this is the opposite of the additive invariant, and is the correct physics
    for lattice doping."""
    base = calculate(50.0, [nbt_endmember()], [], default_reagents(), {},
                     apply_purity=False)
    doped = calculate(50.0, [nbt_endmember()], [_sub("La", 2.0, "A")],
                      default_reagents(), {}, apply_purity=False)
    na_base = next(r for r in base.rows if r.element == "Na")
    na_doped = next(r for r in doped.rows if r.element == "Na")
    # Na coeff dropped 0.5 → 0.49, so its reagent mass drops too
    assert na_doped.coeff == pytest.approx(0.49)
    assert na_doped.mass_target < na_base.mass_target
    # MW changed (La heavier than the Na/Bi it displaced)
    assert doped.total_mw != pytest.approx(base.total_mw)
    # La is an ordinary element row now, NOT an add-on dopant row
    assert any(r.element == "La" for r in doped.rows)
    assert doped.dopant_rows == []


def test_substitutional_delivers_exact_moles_for_all_elements():
    """The core invariant still holds with substitution: every reagent mass
    (host and dopant alike) converts back to exactly moles × coeff."""
    reagents = default_reagents()
    result = calculate(50.0, [nbt_endmember()], [_sub("La", 3.0, "A")],
                       reagents, {}, apply_purity=False)
    for r in result.rows:
        moles_required = result.moles * r.coeff
        moles_delivered = (r.mass_target / r.reagent.mw) * r.reagent.atoms
        assert moles_delivered == pytest.approx(moles_required), r.element


def test_substitutional_b_site():
    """1 mol% Mn on the B-site (Ti) → Ti=0.99, Mn=0.01."""
    a, b = build_composition_coeffs([nbt_endmember()], [_sub("Mn", 1.0, "B")])
    assert b["Ti"] == pytest.approx(0.99)
    assert b["Mn"] == pytest.approx(0.01)
    assert sum(b.values()) == pytest.approx(1.0)


def test_codoping_same_site_shares_remainder():
    """Two substitutional dopants on the A-site: hosts share 1 − (x1+x2), not
    sequentially double-scaled. 1% La + 1% K on A → Na,Bi ×0.98, La=0.01, K=0.01."""
    dops = [_sub("La", 1.0, "A"), _sub("K", 1.0, "A")]
    a, b = build_composition_coeffs([nbt_endmember()], dops)
    assert a["Na"] == pytest.approx(0.5 * 0.98)
    assert a["Bi"] == pytest.approx(0.5 * 0.98)
    assert a["La"] == pytest.approx(0.01)
    assert a["K"] == pytest.approx(0.01)
    assert sum(a.values()) == pytest.approx(1.0)


def test_codoping_both_sites():
    """Donor/acceptor co-doping: 2% La on A and 1% Mn on B simultaneously."""
    dops = [_sub("La", 2.0, "A"), _sub("Mn", 1.0, "B")]
    a, b = build_composition_coeffs([nbt_endmember()], dops)
    assert a["La"] == pytest.approx(0.02)
    assert sum(a.values()) == pytest.approx(1.0)
    assert b["Ti"] == pytest.approx(0.99)
    assert b["Mn"] == pytest.approx(0.01)
    assert sum(b.values()) == pytest.approx(1.0)


def test_substitutional_and_additive_together():
    """A substitutional dopant (folded into the formula) and an additive dopant
    (weighed on top) coexist without interfering."""
    reagents = default_reagents()
    dops = [_sub("La", 2.0, "A"), DopantEntry("Mn", 0.5, unit=UNIT_WT_PCT)]
    result = calculate(50.0, [nbt_endmember()], dops, reagents, {},
                       apply_purity=False)
    # La folded into base rows, Mn is an add-on
    assert any(r.element == "La" for r in result.rows)
    assert [d.cation for d in result.dopant_rows] == ["Mn"]
    # Mn add-on still 0.5 wt% of batch
    assert result.dopant_rows[0].mass == pytest.approx(0.005 * 50.0)


def test_substitutional_dopant_entry_roundtrip():
    d = DopantEntry(cation="La", amount=2.0, mode=MODE_SUBSTITUTIONAL, site="A")
    d2 = DopantEntry.from_dict(d.to_dict())
    assert d2 == d


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
