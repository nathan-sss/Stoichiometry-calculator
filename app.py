"""Streamlit UI for the Perovskite Recipe Calculator.

Run with:  streamlit run app.py
"""

from __future__ import annotations
import base64
import csv
import json
import re
from copy import deepcopy
from datetime import datetime
from io import StringIO

import streamlit as st

from data import (
    ATOMIC_WEIGHTS, A_SITE_CATIONS, B_SITE_CATIONS,
    PERIODIC_TABLE_LAYOUT, DEFAULT_REAGENTS, PRESETS, NON_CATION_ELEMENTS,
)
from calculator import (
    EndMember, Dopant, Reagent, calculate, build_composition_coeffs,
    format_composition, normalize_dict,
)


# ---------- Chemistry formatting helpers (mirrors main.py) ----------

_SUBSCRIPT_DIGITS = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
_UNSUBSCRIPT_DIGITS = str.maketrans("₀₁₂₃₄₅₆₇₈₉", "0123456789")
_FORMULA_DIGIT_RE = re.compile(r"(?<=[A-Za-z\)])(\d+\.?\d*)")


def prettify_formula(s: str) -> str:
    """Subscript digits that follow a letter or ')' so 'Na2CO3' -> 'Na₂CO₃'
    while leading mole-fraction coefficients (e.g. '0.94 ') stay full-size."""
    return _FORMULA_DIGIT_RE.sub(
        lambda m: m.group(0).translate(_SUBSCRIPT_DIGITS), s
    )


def normalize_formula(s: str) -> str:
    """Reverse of prettify_formula — turn Unicode subscript digits back to ASCII
    so storage and JSON serialization stay canonical."""
    return s.translate(_UNSUBSCRIPT_DIGITS)


def render_html_table(rows: list[dict]) -> str:
    """Build a plain HTML table from a list of row-dicts.

    We must NOT hand a DataFrame to st.dataframe / st.table / st.data_editor
    on the stlite web build: their Apache Arrow IPC pipeline is incompatible
    with stlite's bundled Arrow decoder and throws "Parquet error: Repetition
    level must be defined for a primitive type". Plain HTML sidesteps that
    entire path. The result is non-interactive (no sort, no select) but
    visually equivalent for a read-only result display."""
    import html as _html
    if not rows:
        return ""
    cols = list(rows[0].keys())
    parts = [
        '<div style="overflow-x:auto;margin:8px 0">',
        '<table style="width:100%;border-collapse:collapse;font-size:0.95em">',
        '<thead><tr style="background:#f1f5f9;text-align:left">',
    ]
    for c in cols:
        parts.append(
            f'<th style="padding:8px 10px;border:1px solid #e2e8f0;'
            f'font-weight:600;color:#1e293b">{_html.escape(str(c))}</th>'
        )
    parts.append('</tr></thead><tbody>')
    for row in rows:
        parts.append('<tr>')
        for c in cols:
            v = row.get(c)
            cell = "" if v is None else _html.escape(str(v))
            parts.append(
                f'<td style="padding:8px 10px;border:1px solid #e2e8f0;'
                f'color:#334155">{cell}</td>'
            )
        parts.append('</tr>')
    parts.append('</tbody></table></div>')
    return "".join(parts)


def download_link(label: str, data: str, filename: str, mime: str) -> None:
    """Render a download link that works in both local Streamlit and stlite/GitHub Pages.

    Streamlit's native st.download_button generates a /media/<hash>.<ext> URL
    served by the Streamlit server — that pattern 404s under stlite when
    deployed on GitHub Pages (service-worker scope can't cover a subpath
    deployment). A data: URL is fully browser-handled and works everywhere."""
    b64 = base64.b64encode(data.encode("utf-8")).decode("ascii")
    href = f"data:{mime};base64,{b64}"
    st.markdown(
        f'<a href="{href}" download="{filename}" '
        f'style="display:inline-block;padding:0.4em 0.9em;'
        f'background:#4f46e5;color:white;border-radius:6px;'
        f'text-decoration:none;font-weight:600;font-size:0.9em">'
        f'{label}</a>',
        unsafe_allow_html=True,
    )


# ---------- Page setup ----------

st.set_page_config(
    page_title="Perovskite Recipe Calculator",
    page_icon="🧪",
    layout="wide",
)

# Custom CSS for periodic table and general polish
st.markdown("""
<style>
    .pt-grid {
        display: grid;
        grid-template-columns: repeat(18, minmax(36px, 1fr));
        gap: 2px;
        min-width: 700px;
    }
    .pt-cell {
        aspect-ratio: 1 / 1;
        display: flex; align-items: center; justify-content: center;
        font-family: monospace; font-weight: 600;
        font-size: 12px; border-radius: 3px;
    }
    .pt-a { background: #ddd6fe; color: #5b21b6; }
    .pt-b { background: #bae6fd; color: #075985; }
    .pt-other { background: white; color: #475569; border: 1px solid #e2e8f0; }
    .pt-disabled { background: #f1f5f9; color: #cbd5e1; }
    .pt-selected { background: #4f46e5 !important; color: white !important; outline: 2px solid #312e81; }
    .composition-preview {
        font-family: monospace; background: linear-gradient(to right, #eef2ff, #faf5ff);
        padding: 12px; border-radius: 6px; border: 1px solid #c7d2fe;
    }
    .site-badge-a {
        background: #ede9fe; color: #5b21b6; padding: 2px 6px;
        border-radius: 4px; font-size: 11px; font-weight: 600;
    }
    .site-badge-b {
        background: #e0f2fe; color: #075985; padding: 2px 6px;
        border-radius: 4px; font-size: 11px; font-weight: 600;
    }
    .stat-sum-ok { color: #047857; font-weight: 600; }
    .stat-sum-off { color: #b45309; font-weight: 600; }

    /* --- Periodic-table picker (rendered inside st.dialog) --- */
    /* Streamlit's default button padding is too generous for the 18-column
       grid: each cell ends up ~25 px wide and the element symbol wraps
       vertically ("L" stacked above "i"). Make buttons inside the dialog
       tight, monospace, and forbid wrapping. */
    [data-testid="stDialog"] .stButton > button {
        padding: 0.35em 0.1em !important;
        min-width: 0 !important;
        width: 100% !important;
        font-family: ui-monospace, "SF Mono", Menlo, monospace !important;
    }
    [data-testid="stDialog"] .stButton > button p {
        white-space: nowrap !important;
        font-family: ui-monospace, "SF Mono", Menlo, monospace !important;
        font-size: 0.95rem !important;
        margin: 0 !important;
        line-height: 1.1 !important;
        font-weight: 600 !important;
    }
</style>
""", unsafe_allow_html=True)


# ---------- Session state init ----------

def init_state():
    defaults = {
        "reagents": [Reagent.from_dict(r) for r in DEFAULT_REAGENTS],
        # Fresh start: one empty end-member. User loads a preset or builds custom.
        "end_members": [EndMember(name="EM1", fraction=1.0, A={}, B={})],
        "dopant": Dopant(),
        "batch_size": 50.0,
        "reagent_choice": {},   # element -> reagent name
        "excess_pct": {},       # element -> percent
        "apply_purity": True,
        "saved_recipes": [],
        "active_tab": "Calculator",
        # Cation picker dialog state. target is one of:
        #   None                                  — dialog closed
        #   {"em_idx": int, "site": "A"/"B"}     — adding a cation to an EM site
        #   "dopant"                              — picking the dopant cation
        "cation_picker_target": None,
    }
    for k, v in defaults.items():
        if k not in st.session_state:
            st.session_state[k] = v

init_state()


def sync_dominant_fraction():
    """Set EM[0]'s fraction to 1 − Σ(other fractions), clamped at 0.

    Reads widget state (if rendered last rerun) for the truest current values,
    then writes back to both the dataclass and the widget's session-state key
    so the disabled EM1 spinbox shows the updated number."""
    ems = st.session_state.get("end_members", [])
    if not ems:
        return
    others_sum = 0.0
    for i in range(1, len(ems)):
        key = f"em_frac_{i}"
        if key in st.session_state:
            v = float(st.session_state[key])
            ems[i].fraction = v
            others_sum += v
        else:
            others_sum += ems[i].fraction
    new_em0 = max(0.0, 1.0 - others_sum)
    ems[0].fraction = new_em0
    # Pre-populate the disabled widget's state so display matches.
    st.session_state["em_frac_0"] = new_em0


sync_dominant_fraction()


# ---------- Helpers ----------

def get_category(el: str) -> str:
    if el not in ATOMIC_WEIGHTS:
        return "unknown"
    if el in A_SITE_CATIONS:
        return "a"
    if el in B_SITE_CATIONS:
        return "b"
    return "other"


def reagents_by_element_map() -> dict:
    """Build {element: Reagent} from current selection."""
    out = {}
    for el, name in st.session_state.reagent_choice.items():
        r = next((r for r in st.session_state.reagents if r.name == name), None)
        if r is not None:
            out[el] = r
    return out


def auto_pick_reagents(elements: list[str]):
    """Ensure each element has a reagent chosen; auto-pick first match if missing."""
    for el in elements:
        current = st.session_state.reagent_choice.get(el)
        # Validate current choice
        if current and any(r.name == current and r.element == el for r in st.session_state.reagents):
            continue
        candidate = next((r for r in st.session_state.reagents if r.element == el), None)
        if candidate is not None:
            st.session_state.reagent_choice[el] = candidate.name


def render_periodic_table(
    site: str | None,
    selected: str | None = None,
    excluded: list[str] | None = None,
    key_prefix: str = "pt",
) -> str | None:
    """Render an interactive periodic table. Returns clicked element or None."""
    excluded = excluded or []
    clicked = None

    # Build a 2D layout dict: {(row, col): symbol}
    layout = {(row, col): sym for sym, row, col in PERIODIC_TABLE_LAYOUT}
    max_row = 8

    for row in range(1, max_row + 1):
        if row == 7:
            continue  # we don't use row 7 (actinides skipped)
        cols = st.columns(18, gap="small")
        for col_idx in range(1, 19):
            sym = layout.get((row, col_idx))
            with cols[col_idx - 1]:
                if sym is None:
                    # Place a lanthanide indicator at (6, 3)
                    if row == 6 and col_idx == 3:
                        st.markdown("<div style='text-align:center;color:#94a3b8;font-size:11px;font-style:italic;padding-top:8px'>*</div>", unsafe_allow_html=True)
                    else:
                        st.markdown("<div></div>", unsafe_allow_html=True)
                    continue

                cat = get_category(sym)
                is_disabled = (
                    sym in excluded
                    or sym in NON_CATION_ELEMENTS
                    or sym == "H"
                    or sym not in ATOMIC_WEIGHTS
                )
                is_selected = (selected == sym)

                # Style
                if is_disabled:
                    btn_type = "secondary"
                elif is_selected:
                    btn_type = "primary"
                elif (site == "A" and cat == "a") or (site == "B" and cat == "b"):
                    btn_type = "primary"
                else:
                    btn_type = "secondary"

                # Visual hint via help text
                tw = ATOMIC_WEIGHTS.get(sym)
                help_text = f"{sym}" + (f" · {tw} g/mol" if tw else "")
                if cat == "a":
                    help_text += " · typical A-site"
                elif cat == "b":
                    help_text += " · typical B-site"

                if st.button(
                    sym,
                    key=f"{key_prefix}_{sym}_{row}_{col_idx}",
                    type=btn_type,
                    use_container_width=True,
                    disabled=is_disabled,
                    help=help_text,
                ):
                    clicked = sym
    return clicked


# ---------- Cation picker dialog ----------
#
# Renders the periodic table inside a wide modal so it doesn't get squeezed
# by the nested end-member container. Triggered by setting
# `st.session_state.cation_picker_target` and then re-running.

@st.dialog("Pick cation", width="large")
def cation_picker_dialog():
    target = st.session_state.get("cation_picker_target")
    if target is None:
        return

    if target == "dopant":
        site = st.session_state.dopant.site
        selected = st.session_state.dopant.cation
        excluded = None
        st.markdown(f"Picking **dopant cation** for the **{site}-site**.")
    else:
        em_idx = target["em_idx"]
        site = target["site"]
        em = st.session_state.end_members[em_idx]
        site_dict = em.A if site == "A" else em.B
        excluded = list(site_dict.keys())
        selected = None
        st.markdown(f"Adding to **{em.name}** {site}-site. "
                    f"Already-used cations are disabled.")

    clicked = render_periodic_table(
        site=site, selected=selected, excluded=excluded,
        key_prefix="dlg_pt",
    )
    if clicked:
        if target == "dopant":
            st.session_state.dopant.cation = clicked
        else:
            em_idx = target["em_idx"]
            site = target["site"]
            em = st.session_state.end_members[em_idx]
            site_dict = em.A if site == "A" else em.B
            site_dict[clicked] = 0.0
        st.session_state.cation_picker_target = None
        st.rerun()

    if st.button("Cancel", key="dlg_cancel"):
        st.session_state.cation_picker_target = None
        st.rerun()


# ---------- Sidebar: tab switching ----------

with st.sidebar:
    st.markdown("## 🧪 Perovskite Recipe Calculator")
    st.caption("")
    st.divider()
    tab_choice = st.radio(
        "Navigation",
        ["Calculator", "Raw Materials", "Saved Recipes"],
        key="active_tab",
        label_visibility="collapsed",
    )

# ============================================================
# CALCULATOR TAB
# ============================================================

if tab_choice == "Calculator":
    # Open the cation picker modal if a target was set on a previous rerun.
    if st.session_state.cation_picker_target is not None:
        cation_picker_dialog()

    st.title("Perovskite Recipe Calculator")
    st.caption("Build in Ge Wang's lab. Any problem please contact xiangjie.song@manchester.ac.uk")

    # ---- Quick presets ----
    with st.container(border=True):
        st.markdown("**✨ Quick presets** — click to load a common perovskite system")
        cols = st.columns(min(len(PRESETS), 8))
        for i, p in enumerate(PRESETS):
            with cols[i % len(cols)]:
                if st.button(f"**{p['name']}**\n\n{prettify_formula(p['full'])}",
                             key=f"preset_{p['id']}", use_container_width=True):
                    st.session_state.end_members = [EndMember.from_dict(em) for em in p["end_members"]]
                    st.session_state.dopant = Dopant()
                    st.session_state.excess_pct = {}
                    st.session_state.reagent_choice = {}
                    # Wipe stale EM widget keys so newly-loaded fractions take effect.
                    for k in list(st.session_state.keys()):
                        if isinstance(k, str) and (k.startswith("em_frac_") or k.startswith("em_name_")):
                            del st.session_state[k]
                    st.rerun()

    # ---- Composition builder ----
    with st.container(border=True):
        st.markdown("## Composition (ABO₃ perovskite)")

        em_sum = sum(em.fraction for em in st.session_state.end_members)
        a_coeffs, b_coeffs = build_composition_coeffs(st.session_state.end_members, st.session_state.dopant)
        a_sum = sum(a_coeffs.values())
        b_sum = sum(b_coeffs.values())

        col1, col2, col3 = st.columns(3)
        with col1:
            ok = abs(em_sum - 1.0) < 1e-3
            cls = "stat-sum-ok" if ok else "stat-sum-off"
            st.markdown(f"<span class='{cls}'>Σ end-members = {em_sum:.4f}</span>", unsafe_allow_html=True)
        with col2:
            ok = abs(a_sum - 1.0) < 1e-3
            cls = "stat-sum-ok" if ok else "stat-sum-off"
            st.markdown(f"<span class='{cls}'>Σ A-site = {a_sum:.4f}</span>", unsafe_allow_html=True)
        with col3:
            ok = abs(b_sum - 1.0) < 1e-3
            cls = "stat-sum-ok" if ok else "stat-sum-off"
            st.markdown(f"<span class='{cls}'>Σ B-site = {b_sum:.4f}</span>", unsafe_allow_html=True)

        st.markdown("")

        # End-members
        to_delete = None
        for idx, em in enumerate(st.session_state.end_members):
            with st.container(border=True):
                c1, c2, c3 = st.columns([3, 2, 1])
                with c1:
                    new_name = st.text_input("End-member name", value=em.name, key=f"em_name_{idx}")
                    if new_name != em.name:
                        st.session_state.end_members[idx].name = new_name
                with c2:
                    is_dominant = (idx == 0)
                    if is_dominant:
                        # sync_dominant_fraction() pre-populates
                        # st.session_state["em_frac_0"]. Streamlit forbids
                        # combining a Session State pre-set with a value=
                        # argument, so we deliberately omit value= here.
                        st.number_input(
                            "Mole fraction (auto = 1 − Σ others)",
                            step=0.01, format="%.4f",
                            key=f"em_frac_{idx}", disabled=True,
                        )
                    else:
                        new_frac = st.number_input(
                            "Mole fraction", value=float(em.fraction),
                            step=0.01, format="%.4f",
                            key=f"em_frac_{idx}",
                        )
                        if new_frac != em.fraction:
                            st.session_state.end_members[idx].fraction = new_frac
                with c3:
                    st.markdown("&nbsp;", unsafe_allow_html=True)  # spacer
                    # First end-member is the base — can't be removed.
                    if idx > 0:
                        if st.button("Remove", key=f"em_del_{idx}", use_container_width=True):
                            to_delete = idx

                # A-site and B-site editors
                site_col_a, site_col_b = st.columns(2)
                for site_label, site_key, badge_cls, col in [
                    ("A-site", "A", "site-badge-a", site_col_a),
                    ("B-site", "B", "site-badge-b", site_col_b),
                ]:
                    with col:
                        site_dict = getattr(em, site_key)
                        site_sum = sum(site_dict.values())
                        st.markdown(f"<span class='{badge_cls}'>{site_label}</span> "
                                    f"<span style='color:#64748b;font-size:11px'>Σ = {site_sum:.4f}</span>",
                                    unsafe_allow_html=True)

                        # Existing cations
                        cation_to_remove = None
                        for el, val in list(site_dict.items()):
                            cc1, cc2, cc3 = st.columns([1, 3, 1])
                            with cc1:
                                st.markdown(f"<div style='padding-top:8px;font-family:monospace;font-weight:600'>{el}</div>",
                                            unsafe_allow_html=True)
                            with cc2:
                                new_val = st.number_input(
                                    f"{el} coeff",
                                    value=float(val), step=0.01, format="%.4f",
                                    key=f"em_{idx}_{site_key}_{el}",
                                    label_visibility="collapsed",
                                )
                                if new_val != val:
                                    site_dict[el] = new_val
                            with cc3:
                                if st.button("✕", key=f"em_{idx}_{site_key}_del_{el}", use_container_width=True):
                                    cation_to_remove = el
                        if cation_to_remove:
                            del site_dict[cation_to_remove]
                            st.rerun()

                        # Normalize button
                        if abs(site_sum - 1.0) > 1e-3 and len(site_dict) > 0:
                            if st.button(f"Normalize {site_label} to Σ = 1",
                                         key=f"em_{idx}_{site_key}_norm",
                                         use_container_width=True):
                                normed = normalize_dict(site_dict)
                                site_dict.clear()
                                site_dict.update(normed)
                                st.rerun()

                        # Add cation — opens periodic-table modal dialog.
                        if st.button(f"➕ Add {site_label} cation",
                                     key=f"em_{idx}_{site_key}_add_btn",
                                     use_container_width=True):
                            st.session_state.cation_picker_target = {
                                "em_idx": idx, "site": site_key,
                            }
                            st.rerun()

        if to_delete is not None:
            st.session_state.end_members.pop(to_delete)
            st.rerun()

        # End-member controls (Normalize button removed — EM1 auto-balances now)
        bc1, _ = st.columns([1, 5])
        with bc1:
            if st.button("➕ Add end-member"):
                st.session_state.end_members.append(
                    EndMember(name=f"EM{len(st.session_state.end_members) + 1}",
                              fraction=0.05, A={}, B={})
                )
                st.rerun()

        # Live composition preview (with empty-state friendly message)
        _a, _b = build_composition_coeffs(st.session_state.end_members, st.session_state.dopant)
        _all_elements = [e for e in set(list(_a.keys()) + list(_b.keys()))
                         if (_a.get(e, 0) + _b.get(e, 0)) > 1e-12]
        if not _all_elements:
            st.markdown(
                "<div class='composition-preview'><i style='color:#94a3b8'>"
                "Add A-site and B-site cations to begin, or load a preset above.</i></div>",
                unsafe_allow_html=True,
            )
        else:
            comp_str = prettify_formula(
                format_composition(st.session_state.end_members, st.session_state.dopant)
            )
            st.markdown(f"<div class='composition-preview'>{comp_str}</div>", unsafe_allow_html=True)

    # ---- Dopant ----
    with st.expander(f"⚛️ Dopant — {'enabled' if st.session_state.dopant.enabled else 'disabled'}",
                     expanded=st.session_state.dopant.enabled):
        new_enabled = st.checkbox("Enable dopant",
                                  value=st.session_state.dopant.enabled,
                                  key="dopant_enabled_cb")
        if new_enabled != st.session_state.dopant.enabled:
            st.session_state.dopant.enabled = new_enabled
            st.rerun()

        if st.session_state.dopant.enabled:
            dc1, dc2, dc3 = st.columns(3)
            with dc1:
                new_site = st.radio("Substitute on site", ["A", "B"],
                                    index=0 if st.session_state.dopant.site == "A" else 1,
                                    horizontal=True, key="dopant_site_radio")
                if new_site != st.session_state.dopant.site:
                    st.session_state.dopant.site = new_site
                    # Reset cation to a sensible default for the new site
                    if new_site == "A":
                        st.session_state.dopant.cation = A_SITE_CATIONS[0]
                    else:
                        st.session_state.dopant.cation = B_SITE_CATIONS[0]
                    st.rerun()
            with dc2:
                new_level = st.number_input("Doping level (mol fraction)",
                                            value=float(st.session_state.dopant.level),
                                            step=0.001, format="%.4f",
                                            key="dopant_level_input")
                if new_level != st.session_state.dopant.level:
                    st.session_state.dopant.level = new_level
                st.caption(f"= {st.session_state.dopant.level * 100:.3f} mol%")
            with dc3:
                st.markdown(f"**Current cation:** `{st.session_state.dopant.cation}`")
                if st.session_state.dopant.site == "A" and st.session_state.dopant.cation not in A_SITE_CATIONS:
                    st.warning(f"⚠ {st.session_state.dopant.cation} is unusual on A-site — check ionic radius.")
                if st.session_state.dopant.site == "B" and st.session_state.dopant.cation not in B_SITE_CATIONS:
                    st.warning(f"⚠ {st.session_state.dopant.cation} is unusual on B-site — check ionic radius.")

            if st.button("🔬 Pick from periodic table…", key="dopant_pick_btn"):
                st.session_state.cation_picker_target = "dopant"
                st.rerun()

    # ---- Recompute coefficients for downstream sections ----
    a_coeffs, b_coeffs = build_composition_coeffs(st.session_state.end_members, st.session_state.dopant)
    elements = sorted(set(list(a_coeffs.keys()) + list(b_coeffs.keys())))
    elements = [el for el in elements if (a_coeffs.get(el, 0) + b_coeffs.get(el, 0)) > 1e-12]
    auto_pick_reagents(elements)

    # ---- Batch & reagents & excess ----
    with st.container(border=True):
        st.markdown("### Batch & reagents")

        bc1, bc2, bc3, bc4 = st.columns(4)
        with bc1:
            new_batch = st.number_input("Batch size (g)", value=float(st.session_state.batch_size),
                                        step=1.0, format="%.4f", key="batch_size_input")
            if new_batch != st.session_state.batch_size:
                st.session_state.batch_size = new_batch
        with bc2:
            st.session_state.apply_purity = st.checkbox(
                "Apply purity correction",
                value=st.session_state.apply_purity,
                key="apply_purity_cb",
                help="When on, the final mass is divided by purity/100 to compensate for impurities."
            )

        st.markdown("**Reagents & excess additions**")
        st.caption("Excess = % of element moles to compensate volatilization (e.g. 3 for +3% Na)")

        missing = []
        for el in elements:
            candidates = [r for r in st.session_state.reagents if r.element == el]
            on_a = a_coeffs.get(el, 0) > 0
            on_b = b_coeffs.get(el, 0) > 0
            site_tag = ("A" if on_a else "") + ("B" if on_b else "")

            rc1, rc2, rc3 = st.columns([1, 4, 2])
            with rc1:
                st.markdown(f"<div style='padding-top:8px'><b>{el}</b> "
                            f"<span style='font-size:10px;color:#64748b'>[{site_tag}]</span></div>",
                            unsafe_allow_html=True)
            with rc2:
                if not candidates:
                    st.error(f"⚠ No reagent in database for {el}. Add one in 'Raw Materials'.")
                    missing.append(el)
                else:
                    names = [r.name for r in candidates]
                    current = st.session_state.reagent_choice.get(el, names[0])
                    if current not in names:
                        current = names[0]
                    idx = names.index(current)
                    new_choice = st.selectbox(
                        f"Reagent for {el}",
                        names,
                        index=idx,
                        format_func=lambda n: (
                            f"{prettify_formula(n)} "
                            f"(MW {next(r.mw for r in candidates if r.name == n)})"
                        ),
                        key=f"reagent_choice_{el}",
                        label_visibility="collapsed",
                    )
                    st.session_state.reagent_choice[el] = new_choice
            with rc3:
                current_excess = st.session_state.excess_pct.get(el, 0.0)
                new_excess = st.number_input(
                    f"Excess % for {el}",
                    value=float(current_excess),
                    step=0.5, format="%.4g",
                    key=f"excess_{el}",
                    label_visibility="collapsed",
                )
                st.session_state.excess_pct[el] = new_excess

    # ---- Calculation ----
    result = calculate(
        batch_size_g=st.session_state.batch_size,
        end_members=st.session_state.end_members,
        dopant=st.session_state.dopant,
        reagents_by_element=reagents_by_element_map(),
        excess_percent=st.session_state.excess_pct,
        apply_purity=st.session_state.apply_purity,
    )

    # ---- Results ----
    with st.container(border=True):
        st.markdown("### Mass calculation")

        mc1, mc2, mc3, mc4 = st.columns(4)
        if not result.rows:
            mc1.metric("Total MW (g/mol)", "—")
            mc2.metric("Moles", "—")
            mc3.metric("Total reagent (g)", "—")
            label = "Total to weigh (g)" if not st.session_state.apply_purity else "Total w/ purity (g)"
            mc4.metric(label, "—")
        else:
            mc1.metric("Total MW (g/mol)", f"{result.total_mw:.4f}")
            mc2.metric("Moles", f"{result.moles:.4f}")
            mc3.metric("Total reagent (g)", f"{result.total_mass_with_excess:.3f}")
            label = "Total to weigh (g)" if not st.session_state.apply_purity else "Total w/ purity (g)"
            mc4.metric(label, f"{result.total_mass_to_weigh:.3f}")

        # Build dataframe for display
        rows_data = []
        for r in result.rows:
            site = ("A" if r.on_a_site else "") + ("B" if r.on_b_site else "")
            rows_data.append({
                "Element": r.element,
                "Site": site,
                "At. wt": round(r.atomic_weight, 3),
                "Coefficient": round(r.coeff, 4),
                "Excess %": round(r.excess_pct, 4),
                "Reagent": prettify_formula(r.reagent.name) if r.reagent else "—",
                "Reagent MW": round(r.reagent.mw, 2) if r.reagent else None,
                "Purity %": round(r.reagent.purity, 2) if r.reagent else None,
                "Mass target (g)": round(r.mass_target, 4),
                "Mass + excess (g)": round(r.mass_with_excess, 4),
                ("w/ purity (g)" if st.session_state.apply_purity else "To weigh (g)"): round(r.mass_to_weigh, 4),
            })
        # Render the table as raw HTML, NOT via st.dataframe / st.table /
        # st.data_editor. All three of those serialize through Apache Arrow
        # IPC, which is broken on the stlite web build ("Parquet error:
        # Repetition level must be defined for a primitive type"). See
        # render_html_table() for the rationale.
        if rows_data:
            st.markdown(render_html_table(rows_data), unsafe_allow_html=True)

        st.caption(
            f"Oxygen contributes {result.oxygen_coeff:.3f} atoms × 15.999 = "
            f"{result.oxygen_mass:.4f} g/mol to the formula unit MW (supplied by the reagents)."
        )

    # ---- Save & export ----
    with st.container(border=True):
        st.markdown("### Save & export")
        sc1, sc2, sc3 = st.columns([3, 1, 1])
        with sc1:
            recipe_name = st.text_input("Recipe name", placeholder="e.g. NBT-BT 6mol% with La doping",
                                        key="recipe_name_input", label_visibility="collapsed")
        with sc2:
            if st.button("💾 Save recipe", use_container_width=True, disabled=not recipe_name.strip()):
                recipe = {
                    "id": int(datetime.now().timestamp() * 1000),
                    "name": recipe_name,
                    "timestamp": datetime.now().isoformat(),
                    "batch_size": st.session_state.batch_size,
                    "end_members": [em.to_dict() for em in st.session_state.end_members],
                    "dopant": st.session_state.dopant.to_dict(),
                    "reagent_choice": dict(st.session_state.reagent_choice),
                    "excess_pct": dict(st.session_state.excess_pct),
                    "apply_purity": st.session_state.apply_purity,
                }
                st.session_state.saved_recipes.insert(0, recipe)
                st.success(f"Saved '{recipe_name}'")

        with sc3:
            # Build CSV
            buf = StringIO()
            w = csv.writer(buf)
            w.writerow(["Recipe Export", datetime.now().isoformat()])
            w.writerow(["Composition", format_composition(st.session_state.end_members, st.session_state.dopant)])
            w.writerow(["Batch size (g)", st.session_state.batch_size])
            w.writerow(["Purity correction", "applied" if st.session_state.apply_purity else "not applied"])
            w.writerow(["Total MW", f"{result.total_mw:.4f}"])
            w.writerow(["Moles", f"{result.moles:.6f}"])
            w.writerow([])
            w.writerow(["Element", "Site", "At wt", "Coeff", "Excess %",
                        "Reagent", "Reagent MW", "Purity %",
                        "Mass target (g)", "Mass + excess (g)",
                        "Purity-corrected (g)" if st.session_state.apply_purity else "To weigh (g)"])
            for r in result.rows:
                site = ("A" if r.on_a_site else "") + ("B" if r.on_b_site else "")
                w.writerow([r.element, site, r.atomic_weight, f"{r.coeff:.6f}",
                            r.excess_pct,
                            r.reagent.name if r.reagent else "",
                            r.reagent.mw if r.reagent else "",
                            r.reagent.purity if r.reagent else "",
                            f"{r.mass_target:.4f}", f"{r.mass_with_excess:.4f}",
                            f"{r.mass_to_weigh:.4f}"])
            download_link(
                "⬇️ Export CSV",
                data=buf.getvalue(),
                filename=f"{recipe_name or 'recipe'}_{int(datetime.now().timestamp())}.csv",
                mime="text/csv",
            )


# ============================================================
# RAW MATERIALS TAB
# ============================================================

elif tab_choice == "Raw Materials":
    st.title("📚 Raw materials database")
    st.caption("Add or edit the precursors your lab uses.")

    # Manual edit form. We can't use st.data_editor here: it routes through
    # the same broken Apache Arrow pipeline as st.dataframe under stlite,
    # giving the same "Parquet error" crash. Rendering each cell as a plain
    # widget bypasses Arrow entirely.
    st.caption(
        "Names display with subscripts; type ASCII (e.g. `Na2CO3`) or "
        "Unicode (e.g. `Na₂CO₃`) — both are accepted. "
        "Click **Apply changes** when done."
    )

    _RM_COL_WIDTHS = [3, 1.5, 1, 1.7, 1.4, 4, 0.7]
    _RM_HEADERS = ["Name", "Element", "Atoms", "MW (g/mol)", "Purity %", "Notes", ""]
    hc = st.columns(_RM_COL_WIDTHS)
    for header, col in zip(_RM_HEADERS, hc):
        col.markdown(f"**{header}**" if header else "")

    to_delete = None
    for i, r in enumerate(st.session_state.reagents):
        cc = st.columns(_RM_COL_WIDTHS)
        cc[0].text_input(
            "Name", value=prettify_formula(r.name),
            key=f"rm_name_{i}", label_visibility="collapsed",
        )
        cc[1].text_input(
            "Element", value=r.element,
            key=f"rm_el_{i}", label_visibility="collapsed",
        )
        cc[2].number_input(
            "Atoms", value=int(r.atoms), min_value=1, step=1,
            key=f"rm_atoms_{i}", label_visibility="collapsed",
        )
        cc[3].number_input(
            "MW", value=float(r.mw), min_value=0.0, step=0.01, format="%.4f",
            key=f"rm_mw_{i}", label_visibility="collapsed",
        )
        cc[4].number_input(
            "Purity", value=float(r.purity), min_value=0.0, max_value=100.0,
            step=0.1, format="%.3f",
            key=f"rm_pur_{i}", label_visibility="collapsed",
        )
        cc[5].text_input(
            "Notes", value=r.notes or "",
            key=f"rm_notes_{i}", label_visibility="collapsed",
        )
        if cc[6].button("✕", key=f"rm_del_{i}", help="Delete this row"):
            to_delete = i

    def _clear_rm_widget_keys():
        """Wipe all rm_* widget keys so the next render re-syncs from
        st.session_state.reagents instead of stale per-index widget state."""
        for k in list(st.session_state.keys()):
            if isinstance(k, str) and k.startswith("rm_"):
                del st.session_state[k]

    if to_delete is not None:
        del st.session_state.reagents[to_delete]
        _clear_rm_widget_keys()
        st.rerun()

    st.divider()
    bc = st.columns([1, 1, 4])
    if bc[0].button("➕ Add reagent"):
        st.session_state.reagents.append(Reagent("New", "X", 1, 0.0, 99.9, ""))
        st.rerun()
    if bc[1].button("💾 Apply changes", type="primary"):
        new_reagents = []
        for i in range(len(st.session_state.reagents)):
            try:
                name = (st.session_state.get(f"rm_name_{i}") or "").strip()
                element = (st.session_state.get(f"rm_el_{i}") or "").strip()
                if not name or not element:
                    continue
                new_reagents.append(Reagent(
                    name=normalize_formula(name),
                    element=normalize_formula(element),
                    atoms=int(st.session_state.get(f"rm_atoms_{i}", 1) or 1),
                    mw=float(st.session_state.get(f"rm_mw_{i}", 0.0) or 0.0),
                    purity=float(st.session_state.get(f"rm_pur_{i}", 99.9) or 99.9),
                    notes=str(st.session_state.get(f"rm_notes_{i}", "") or ""),
                ))
            except (ValueError, TypeError):
                continue
        st.session_state.reagents = new_reagents
        _clear_rm_widget_keys()
        st.success(f"Updated database: {len(new_reagents)} reagents")
        st.rerun()

    st.divider()
    st.caption(
        "**Element** = which target element this reagent supplies. "
        "**Atoms** = atoms of that element per reagent formula unit (e.g. Na₂CO₃ → element=Na, atoms=2)."
    )

    # Import/Export JSON
    st.markdown("### Backup / share")
    c1, c2 = st.columns(2)
    with c1:
        backup = json.dumps([r.to_dict() for r in st.session_state.reagents], indent=2)
        download_link(
            "⬇️ Download as JSON",
            data=backup,
            filename=f"reagents_{int(datetime.now().timestamp())}.json",
            mime="application/json",
        )
    with c2:
        uploaded = st.file_uploader("Import reagents from JSON", type="json", key="reagents_upload")
        if uploaded is not None:
            try:
                data = json.load(uploaded)
                st.session_state.reagents = [Reagent.from_dict(r) for r in data]
                st.success(f"Imported {len(data)} reagents")
                st.rerun()
            except Exception as e:
                st.error(f"Failed to import: {e}")


# ============================================================
# SAVED RECIPES TAB
# ============================================================

elif tab_choice == "Saved Recipes":
    st.title("📖 Saved recipes")

    if not st.session_state.saved_recipes:
        st.info("No saved recipes yet. Save one from the Calculator tab.")
    else:
        for recipe in st.session_state.saved_recipes:
            with st.container(border=True):
                rc1, rc2, rc3 = st.columns([4, 1, 1])
                with rc1:
                    st.markdown(f"**{recipe['name']}**")
                    em_list = [EndMember.from_dict(em) for em in recipe["end_members"]]
                    dop = Dopant.from_dict(recipe["dopant"])
                    st.code(prettify_formula(format_composition(em_list, dop)), language=None)
                    st.caption(f"{recipe['batch_size']:.2f} g · {recipe['timestamp']}")
                with rc2:
                    if st.button("📥 Load", key=f"load_{recipe['id']}", use_container_width=True):
                        st.session_state.batch_size = recipe["batch_size"]
                        st.session_state.end_members = [EndMember.from_dict(em) for em in recipe["end_members"]]
                        st.session_state.dopant = Dopant.from_dict(recipe["dopant"])
                        st.session_state.reagent_choice = dict(recipe["reagent_choice"])
                        st.session_state.excess_pct = dict(recipe["excess_pct"])
                        st.session_state.apply_purity = recipe.get("apply_purity", True)
                        # Clear stale EM widget keys so new fractions render correctly.
                        for k in list(st.session_state.keys()):
                            if isinstance(k, str) and (k.startswith("em_frac_") or k.startswith("em_name_")):
                                del st.session_state[k]
                        st.session_state.active_tab = "Calculator"
                        st.rerun()
                with rc3:
                    if st.button("🗑️ Delete", key=f"del_{recipe['id']}", use_container_width=True):
                        st.session_state.saved_recipes = [
                            r for r in st.session_state.saved_recipes if r["id"] != recipe["id"]
                        ]
                        st.rerun()

        st.divider()
        # Export all recipes as JSON
        all_json = json.dumps(st.session_state.saved_recipes, indent=2)
        download_link(
            "⬇️ Export all recipes as JSON",
            data=all_json,
            filename=f"recipes_{int(datetime.now().timestamp())}.json",
            mime="application/json",
        )

        uploaded = st.file_uploader("Import recipes from JSON", type="json", key="recipes_upload")
        if uploaded is not None:
            try:
                data = json.load(uploaded)
                st.session_state.saved_recipes = data + st.session_state.saved_recipes
                st.success(f"Imported {len(data)} recipes")
                st.rerun()
            except Exception as e:
                st.error(f"Failed to import: {e}")
