"""Desktop UI (PySide6) for the Perovskite Recipe Calculator — Phase 1.

Run with:           python main.py
Build .exe (Win):   pyinstaller --onefile --windowed --name PerovskiteCalc main.py
"""
from __future__ import annotations

import csv
import json
import re
import sys
from datetime import datetime
from pathlib import Path

from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import QFont
from PySide6.QtWidgets import (
    QAbstractItemView, QApplication, QCheckBox, QComboBox,
    QDialog, QDoubleSpinBox, QFileDialog, QGridLayout, QGroupBox,
    QHBoxLayout, QHeaderView, QLabel, QLineEdit, QMainWindow, QMessageBox,
    QPushButton, QScrollArea, QTableWidget, QTableWidgetItem,
    QTabWidget, QVBoxLayout, QWidget,
)

from calculator import (
    DopantEntry, EndMember, Reagent, build_composition_coeffs, calculate,
    format_composition, dopant_entry_from_legacy_dict,
    UNIT_MOL_PCT, UNIT_WT_PCT,
)
from data import (
    A_SITE_CATIONS, ATOMIC_WEIGHTS, B_SITE_CATIONS, DEFAULT_REAGENTS,
    NON_CATION_ELEMENTS, PERIODIC_TABLE_LAYOUT, PRESETS,
)


ALL_ELEMENTS = sorted(ATOMIC_WEIGHTS.keys())


class NoWheelDoubleSpinBox(QDoubleSpinBox):
    """Spinbox that ignores scroll-wheel events so values don't drift when
    scrolling the page."""
    def wheelEvent(self, event):
        event.ignore()


class NoWheelComboBox(QComboBox):
    """ComboBox that ignores scroll-wheel events so the selection doesn't
    change when scrolling the page."""
    def wheelEvent(self, event):
        event.ignore()


def make_spin(value: float, step: float = 0.001, decimals: int = 4,
              minimum: float = 0.0, maximum: float = 1e6) -> NoWheelDoubleSpinBox:
    s = NoWheelDoubleSpinBox()
    s.setRange(minimum, maximum)
    s.setDecimals(decimals)
    s.setSingleStep(step)
    s.setValue(value)
    return s


_SUBSCRIPT_DIGITS = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
_UNSUBSCRIPT_DIGITS = str.maketrans("₀₁₂₃₄₅₆₇₈₉", "0123456789")
# Match digits (with optional decimal) that follow a letter or closing paren,
# so leading mole-fraction coefficients (like "0.94 ") stay full-size.
_FORMULA_DIGIT_RE = re.compile(r"(?<=[A-Za-z\)])(\d+\.?\d*)")


def prettify_formula(s: str) -> str:
    """Convert ASCII digits inside a chemical formula to Unicode subscripts.

    Only digits that follow a letter or closing parenthesis are converted,
    so leading coefficients are left readable.

    >>> prettify_formula("0.94 (Bi0.5Na0.5)(Ti)O3 - 0.06 (Ba)(Ti)O3")
    '0.94 (Bi₀.₅Na₀.₅)(Ti)O₃ - 0.06 (Ba)(Ti)O₃'
    >>> prettify_formula("Na2CO3")
    'Na₂CO₃'
    """
    return _FORMULA_DIGIT_RE.sub(
        lambda m: m.group(0).translate(_SUBSCRIPT_DIGITS), s
    )


def normalize_formula(s: str) -> str:
    """Inverse of prettify_formula: turn Unicode subscript digits back to ASCII
    so storage and serialization stay in their canonical form."""
    return s.translate(_UNSUBSCRIPT_DIGITS)


# ---------- Periodic-table picker ----------

class PeriodicTablePicker(QDialog):
    """Modal periodic-table cation picker.

    `site` ("A" / "B" / None) highlights typical cations for that site.
    `excluded` greys out already-used symbols. `selected` outlines the current.
    Use the classmethod `pick(...)` for one-shot use.
    """

    A_BG, A_FG = "#ede9fe", "#6d28d9"
    B_BG, B_FG = "#e0f2fe", "#075985"
    A_HL_BG, A_HL_FG = "#ddd6fe", "#5b21b6"
    B_HL_BG, B_HL_FG = "#bae6fd", "#075985"
    OTHER_BG, OTHER_FG = "#ffffff", "#475569"
    DISABLED_BG, DISABLED_FG = "#f1f5f9", "#cbd5e1"

    def __init__(self, site: str | None = None,
                 excluded: list[str] | None = None,
                 selected: str | None = None,
                 parent=None):
        super().__init__(parent)
        self.setWindowTitle("Pick cation")
        self.setModal(True)
        self.result_symbol: str | None = None
        excluded_set = set(excluded or [])

        outer = QVBoxLayout(self)

        legend = QHBoxLayout()
        for label, bg, fg in [
            ("A-site", self.A_BG, self.A_FG),
            ("B-site", self.B_BG, self.B_FG),
            ("Other cation", self.OTHER_BG, self.OTHER_FG),
            ("Disabled", self.DISABLED_BG, self.DISABLED_FG),
        ]:
            chip = QLabel(label)
            chip.setStyleSheet(
                f"background:{bg}; color:{fg}; padding:3px 8px;"
                f"border:1px solid #cbd5e1; border-radius:3px; font-size:12px;"
            )
            legend.addWidget(chip)
        legend.addStretch(1)
        outer.addLayout(legend)

        grid = QGridLayout()
        grid.setHorizontalSpacing(2)
        grid.setVerticalSpacing(2)

        for sym, row, col in PERIODIC_TABLE_LAYOUT:
            in_a = sym in A_SITE_CATIONS
            in_b = sym in B_SITE_CATIONS
            cat = "a" if in_a else ("b" if in_b else "other")
            disabled = (
                sym in NON_CATION_ELEMENTS
                or sym in excluded_set
                or sym == "H"
                or sym not in ATOMIC_WEIGHTS
            )

            btn = QPushButton(sym)
            btn.setFixedSize(48, 48)
            btn.setEnabled(not disabled)
            btn.setStyleSheet(self._style_for(cat, site, disabled, sym == selected))
            tw = ATOMIC_WEIGHTS.get(sym)
            tip = sym + (f" — {tw} g/mol" if tw else "")
            if in_a:
                tip += " · typical A-site"
            elif in_b:
                tip += " · typical B-site"
            if sym in excluded_set:
                tip += " · already used"
            btn.setToolTip(tip)
            btn.clicked.connect(lambda _checked=False, s=sym: self._pick(s))
            grid.addWidget(btn, row, col)

        # Lanthanide indicator at row 6, col 3
        ind = QLabel("*")
        ind.setStyleSheet("color:#94a3b8; font-style:italic;")
        ind.setAlignment(Qt.AlignCenter)
        grid.addWidget(ind, 6, 3)

        # Visual gap between row 6 and the lanthanide row
        grid.setRowMinimumHeight(7, 10)

        outer.addLayout(grid)

        footer = QHBoxLayout()
        footer.addStretch(1)
        cancel = QPushButton("Cancel")
        cancel.clicked.connect(self.reject)
        footer.addWidget(cancel)
        outer.addLayout(footer)

    def _style_for(self, cat: str, site: str | None,
                   disabled: bool, is_selected: bool) -> str:
        if disabled:
            return (
                f"QPushButton {{ background:{self.DISABLED_BG}; color:{self.DISABLED_FG};"
                f"border:1px solid #e2e8f0; border-radius:3px;"
                f"font-family:monospace; font-weight:600; }}"
            )
        # Highlight if matches the requested site
        highlight = (site == "A" and cat == "a") or (site == "B" and cat == "b")
        if highlight:
            bg = self.A_HL_BG if cat == "a" else self.B_HL_BG
            fg = self.A_HL_FG if cat == "a" else self.B_HL_FG
        elif cat == "a":
            bg, fg = self.A_BG, self.A_FG
        elif cat == "b":
            bg, fg = self.B_BG, self.B_FG
        else:
            bg, fg = self.OTHER_BG, self.OTHER_FG

        border = "2px solid #4f46e5" if is_selected else "1px solid #e2e8f0"
        return (
            f"QPushButton {{ background:{bg}; color:{fg};"
            f"border:{border}; border-radius:3px;"
            f"font-family:monospace; font-weight:600; }}"
            f"QPushButton:hover {{ background:#f1f5f9; }}"
        )

    def _pick(self, symbol: str):
        self.result_symbol = symbol
        self.accept()

    @classmethod
    def pick(cls, parent, site: str | None = None,
             excluded: list[str] | None = None,
             selected: str | None = None) -> str | None:
        dlg = cls(site=site, excluded=excluded, selected=selected, parent=parent)
        if dlg.exec() == QDialog.Accepted:
            return dlg.result_symbol
        return None


# ---------- Cation row ----------

class CationRow(QWidget):
    changed = Signal()
    remove_requested = Signal(object)

    def __init__(self, cation: str, coeff: float, parent=None):
        super().__init__(parent)
        h = QHBoxLayout(self)
        h.setContentsMargins(0, 0, 0, 0)

        self.cation_combo = NoWheelComboBox()
        self.cation_combo.setEditable(True)
        self.cation_combo.addItems(ALL_ELEMENTS)
        self.cation_combo.setCurrentText(cation)
        self.cation_combo.setMinimumWidth(70)
        self.cation_combo.setMaximumWidth(90)

        self.coeff_spin = make_spin(coeff, step=0.01, decimals=4, maximum=10.0)
        self.coeff_spin.setMinimumWidth(90)

        remove_btn = QPushButton("✕")
        remove_btn.setFixedWidth(28)

        h.addWidget(self.cation_combo)
        h.addWidget(self.coeff_spin, stretch=1)
        h.addWidget(remove_btn)

        self.cation_combo.activated.connect(self.changed)
        self.cation_combo.lineEdit().editingFinished.connect(self.changed)
        self.coeff_spin.valueChanged.connect(self.changed)
        remove_btn.clicked.connect(lambda: self.remove_requested.emit(self))

    def cation(self) -> str:
        return self.cation_combo.currentText().strip()

    def coeff(self) -> float:
        return self.coeff_spin.value()


# ---------- End-member card ----------

class EndMemberCard(QGroupBox):
    changed = Signal()
    remove_requested = Signal(object)

    def __init__(self, em: EndMember, is_dominant: bool = False, parent=None):
        super().__init__(parent)
        self.is_dominant = is_dominant
        v = QVBoxLayout(self)

        head = QHBoxLayout()
        head.addWidget(QLabel("Name:"))
        self.name_edit = QLineEdit(em.name)
        self.name_edit.setMaximumWidth(180)
        head.addWidget(self.name_edit)
        head.addSpacing(12)
        head.addWidget(QLabel("Mole fraction:"))
        self.frac_spin = make_spin(em.fraction, step=0.01, decimals=4, maximum=10.0)
        self.frac_spin.setMaximumWidth(110)
        head.addWidget(self.frac_spin)
        if is_dominant:
            self.frac_spin.setEnabled(False)
            self.frac_spin.setToolTip("Auto-calculated as 1 − Σ(other end-members)")
            auto_hint = QLabel(
                "<span style='color:#94a3b8;font-size:12px;font-style:italic'>"
                "auto = 1 − Σ others</span>"
            )
            head.addWidget(auto_hint)
        head.addStretch(1)
        remove_btn = QPushButton("🗑  Remove")
        if is_dominant:
            remove_btn.setEnabled(False)
            remove_btn.setToolTip("The base end-member can't be removed")
        head.addWidget(remove_btn)
        v.addLayout(head)

        sites = QHBoxLayout()
        self.a_box, self.a_rows_layout, self.a_sum_label = self._site_box("A-site")
        self.b_box, self.b_rows_layout, self.b_sum_label = self._site_box("B-site")
        sites.addWidget(self.a_box)
        sites.addWidget(self.b_box)
        v.addLayout(sites)

        self.a_rows: list[CationRow] = []
        self.b_rows: list[CationRow] = []
        for el, c in em.A.items():
            self._add_cation_row("A", el, c)
        for el, c in em.B.items():
            self._add_cation_row("B", el, c)
        self._update_sum_labels()

        self.name_edit.editingFinished.connect(self.changed)
        self.frac_spin.valueChanged.connect(self._emit_changed)
        remove_btn.clicked.connect(lambda: self.remove_requested.emit(self))

    def _site_box(self, label: str):
        box = QGroupBox(label)
        layout = QVBoxLayout(box)
        sum_label = QLabel("Σ = 0.000")
        layout.addWidget(sum_label)
        rows_layout = QVBoxLayout()
        layout.addLayout(rows_layout)
        btns = QHBoxLayout()
        add_btn = QPushButton(f"+ Add cation")
        norm_btn = QPushButton("Normalize Σ→1")
        btns.addWidget(add_btn)
        btns.addWidget(norm_btn)
        layout.addLayout(btns)
        add_btn.clicked.connect(lambda: self._pick_and_add_cation(label[0]))
        norm_btn.clicked.connect(lambda: self._normalize_site(label[0]))
        return box, rows_layout, sum_label

    def _pick_and_add_cation(self, site: str):
        existing = [r.cation() for r in (self.a_rows if site == "A" else self.b_rows)]
        picked = PeriodicTablePicker.pick(self, site=site, excluded=existing)
        if picked:
            self._add_cation_row(site, cation=picked, coeff=0.0)

    def _add_cation_row(self, site: str, cation: str | None = None, coeff: float = 0.0):
        defaults = A_SITE_CATIONS if site == "A" else B_SITE_CATIONS
        existing = {r.cation() for r in (self.a_rows if site == "A" else self.b_rows)}
        if cation is None:
            cation = next((c for c in defaults if c not in existing), defaults[0])
        row = CationRow(cation, coeff)
        row.changed.connect(self._emit_changed)
        row.remove_requested.connect(lambda r, s=site: self._remove_cation_row(s, r))
        rows = self.a_rows if site == "A" else self.b_rows
        layout = self.a_rows_layout if site == "A" else self.b_rows_layout
        rows.append(row)
        layout.addWidget(row)
        self._emit_changed()

    def _remove_cation_row(self, site: str, row: CationRow):
        rows = self.a_rows if site == "A" else self.b_rows
        if row in rows:
            rows.remove(row)
            row.setParent(None)
            row.deleteLater()
        self._emit_changed()

    def _normalize_site(self, site: str):
        rows = self.a_rows if site == "A" else self.b_rows
        s = sum(r.coeff() for r in rows)
        if s <= 0:
            return
        for r in rows:
            r.coeff_spin.blockSignals(True)
            r.coeff_spin.setValue(r.coeff() / s)
            r.coeff_spin.blockSignals(False)
        self._emit_changed()

    def _update_sum_labels(self):
        a_sum = sum(r.coeff() for r in self.a_rows)
        b_sum = sum(r.coeff() for r in self.b_rows)
        for lbl, val in ((self.a_sum_label, a_sum), (self.b_sum_label, b_sum)):
            ok = abs(val - 1.0) < 1e-3
            color = "#047857" if ok else "#b45309"
            lbl.setText(f"Σ = {val:.4f}")
            lbl.setStyleSheet(f"color: {color}; font-size: 12px; font-weight: 600;")

    def _emit_changed(self):
        self._update_sum_labels()
        self.changed.emit()

    def to_endmember(self) -> EndMember:
        A: dict[str, float] = {}
        for r in self.a_rows:
            c = r.cation()
            if c:
                A[c] = A.get(c, 0.0) + r.coeff()
        B: dict[str, float] = {}
        for r in self.b_rows:
            c = r.cation()
            if c:
                B[c] = B.get(c, 0.0) + r.coeff()
        return EndMember(name=self.name_edit.text() or "EM",
                         fraction=self.frac_spin.value(), A=A, B=B)


# ---------- Main window ----------

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Perovskite Recipe Calculator")
        self.resize(1100, 850)

        self.reagents: list[Reagent] = [Reagent.from_dict(r) for r in DEFAULT_REAGENTS]
        self.reagent_choice: dict[str, str] = {}
        self.excess_pct: dict[str, float] = {}
        self.end_member_cards: list[EndMemberCard] = []
        self.dopants: list[DopantEntry] = []
        self.saved_recipes: list[dict] = []
        self._suppress = False
        self._current_reagent_sig: tuple | None = None

        self.tabs = QTabWidget()
        self.setCentralWidget(self.tabs)
        self.tabs.addTab(self._build_calculator_tab(), "Calculator")
        self.tabs.addTab(self._build_raw_materials_tab(), "Raw materials")
        self.tabs.addTab(self._build_saved_recipes_tab(), "Saved recipes")

        # Start with one empty dominant end-member; user can load a preset
        # or build a custom composition from scratch.
        self._add_end_member_card(EndMember(name="EM1", fraction=1.0, A={}, B={}),
                                  is_dominant=True)
        self._recalculate()

    def _build_calculator_tab(self) -> QWidget:
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        body = QWidget()
        scroll.setWidget(body)
        self.body_layout = QVBoxLayout(body)
        self.body_layout.setContentsMargins(16, 16, 16, 16)
        self.body_layout.setSpacing(12)

        self._build_header()
        self._build_presets()
        self._build_composition()
        self._build_dopant()
        self._build_batch_reagents()
        self._build_results()
        self._build_export()
        self.body_layout.addStretch(1)
        return scroll

    # ----- builders -----

    def _build_header(self):
        title = QLabel("Calculator")
        f = QFont(); f.setPointSize(18); f.setBold(True)
        title.setFont(f)
        self.body_layout.addWidget(title)
        sub = QLabel("Build for Ge Wang's lab. Any problem please contact xiangjie.song@manchester.ac.uk")
        sub.setStyleSheet("color: #a032a8;")
        self.body_layout.addWidget(sub)

    def _build_presets(self):
        box = QGroupBox("Quick presets")
        h = QHBoxLayout(box)
        self.preset_combo = NoWheelComboBox()
        for p in PRESETS:
            self.preset_combo.addItem(f"{p['name']} — {prettify_formula(p['full'])}", userData=p)
        load_btn = QPushButton("Load preset")
        load_btn.clicked.connect(lambda: self._load_preset(self.preset_combo.currentData()))
        h.addWidget(self.preset_combo, stretch=1)
        h.addWidget(load_btn)
        self.body_layout.addWidget(box)

    def _build_composition(self):
        box = QGroupBox("Composition (ABO₃ perovskite)")
        v = QVBoxLayout(box)

        sums = QHBoxLayout()
        self.em_sum_label = QLabel()
        self.a_sum_label = QLabel()
        self.b_sum_label = QLabel()
        for lbl in (self.em_sum_label, self.a_sum_label, self.b_sum_label):
            sums.addWidget(lbl)
        sums.addStretch(1)
        v.addLayout(sums)

        self.cards_layout = QVBoxLayout()
        v.addLayout(self.cards_layout)

        btns = QHBoxLayout()
        add_btn = QPushButton("+ Add end-member")
        add_btn.clicked.connect(self._add_blank_end_member)
        btns.addWidget(add_btn)
        btns.addStretch(1)
        v.addLayout(btns)

        self.composition_preview = QLabel()
        self.composition_preview.setStyleSheet(
            "font-family: monospace; background: #eef2ff;"
            "padding: 10px; border-radius: 4px; border: 1px solid #c7d2fe;"
        )
        self.composition_preview.setWordWrap(True)
        v.addWidget(self.composition_preview)

        self.body_layout.addWidget(box)

    def _build_dopant(self):
        box = QGroupBox("Dopants (added on top of the batch — not substituted)")
        v = QVBoxLayout(box)

        caption = QLabel(
            "Each dopant is weighed in addition to the base batch. Use "
            "<b>mol %</b> for excess of the cation relative to the formula-unit moles, "
            "or <b>wt %</b> for percent of the base batch mass."
        )
        caption.setWordWrap(True)
        caption.setStyleSheet("color: #64748b; font-size: 12px;")
        v.addWidget(caption)

        # Header row
        head = QHBoxLayout()
        for label, w in [("Cation", 70), ("Amount", 110), ("Unit", 110), ("", 32)]:
            lbl = QLabel(f"<b>{label}</b>" if label else "")
            lbl.setMinimumWidth(w)
            head.addWidget(lbl)
        head.addStretch(1)
        v.addLayout(head)

        self.dopant_rows_layout = QVBoxLayout()
        v.addLayout(self.dopant_rows_layout)

        btn_row = QHBoxLayout()
        add_btn = QPushButton("+ Add dopant")
        add_btn.clicked.connect(self._add_dopant_via_picker)
        btn_row.addWidget(add_btn)
        btn_row.addStretch(1)
        v.addLayout(btn_row)

        self._dopant_row_widgets: list[QWidget] = []
        self.body_layout.addWidget(box)

    def _build_batch_reagents(self):
        box = QGroupBox("Batch & reagents")
        v = QVBoxLayout(box)

        top = QHBoxLayout()
        top.addWidget(QLabel("Batch size (g):"))
        self.batch_spin = make_spin(50.0, step=1.0, decimals=4, maximum=100000.0)
        self.batch_spin.setMinimumWidth(130)
        top.addWidget(self.batch_spin)
        top.addSpacing(20)
        self.purity_cb = QCheckBox("Apply purity correction")
        self.purity_cb.setChecked(True)
        top.addWidget(self.purity_cb)
        top.addStretch(1)
        v.addLayout(top)

        v.addWidget(QLabel("Reagent & excess % per element:"))
        self.reagents_layout = QVBoxLayout()
        v.addLayout(self.reagents_layout)

        self.batch_spin.valueChanged.connect(self._recalculate)
        self.purity_cb.toggled.connect(self._recalculate)

        self.body_layout.addWidget(box)

    def _build_results(self):
        box = QGroupBox("Mass calculation")
        v = QVBoxLayout(box)

        metrics = QHBoxLayout()
        self.mw_label = QLabel("Total MW: —")
        self.moles_label = QLabel("Moles: —")
        self.mass_excess_label = QLabel("Total mass: —")
        self.mass_weigh_label = QLabel("To weigh: —")
        for lbl in (self.mw_label, self.moles_label, self.mass_excess_label, self.mass_weigh_label):
            lbl.setStyleSheet("font-weight: 600; padding: 4px 10px; background: #f8fafc; border: 1px solid #e2e8f0; border-radius: 4px;")
            metrics.addWidget(lbl)
        metrics.addStretch(1)
        v.addLayout(metrics)

        self.results_table = QTableWidget(0, 9)
        self.results_table.setHorizontalHeaderLabels([
            "Element", "Site", "At. wt", "Coeff", "Excess (mol %)",
            "Reagent", "Reagent MW", "Purity %", "To weigh (g)"
        ])
        self.results_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.results_table.setEditTriggers(QTableWidget.NoEditTriggers)
        self.results_table.setMinimumHeight(180)
        v.addWidget(self.results_table)

        # Add-on dopants table (hidden until any dopants are present)
        self.dopants_section_label = QLabel(
            "<b>Add-on dopants</b> (weighed in addition to the base batch)"
        )
        self.dopants_section_label.setStyleSheet("padding-top: 6px;")
        self.dopants_section_label.setVisible(False)
        v.addWidget(self.dopants_section_label)

        self.dopants_table = QTableWidget(0, 7)
        self.dopants_table.setHorizontalHeaderLabels([
            "Cation", "Amount", "At. wt",
            "Reagent", "Reagent MW", "Purity %",
            "To weigh (g)",
        ])
        self.dopants_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.dopants_table.setEditTriggers(QTableWidget.NoEditTriggers)
        self.dopants_table.setMinimumHeight(80)
        self.dopants_table.setVisible(False)
        v.addWidget(self.dopants_table)

        self.oxygen_label = QLabel()
        self.oxygen_label.setStyleSheet("color: #64748b; font-size: 12px;")
        v.addWidget(self.oxygen_label)

        self.body_layout.addWidget(box)

    def _build_export(self):
        box = QGroupBox("Save & export")
        h = QHBoxLayout(box)
        h.addWidget(QLabel("Recipe name:"))
        self.recipe_name = QLineEdit()
        self.recipe_name.setPlaceholderText("e.g. NBT-BT 6mol% with La doping")
        h.addWidget(self.recipe_name, stretch=1)
        save_btn = QPushButton("💾 Save recipe")
        save_btn.clicked.connect(self._save_current_recipe)
        h.addWidget(save_btn)
        export_btn = QPushButton("⬇ Export CSV…")
        export_btn.clicked.connect(self._export_csv)
        h.addWidget(export_btn)
        self.body_layout.addWidget(box)

    # ----- preset / end-member ops -----

    def _load_preset(self, preset: dict | None):
        if preset is None:
            return
        self._suppress = True
        for card in list(self.end_member_cards):
            card.setParent(None)
            card.deleteLater()
        self.end_member_cards.clear()
        for i, em_dict in enumerate(preset["end_members"]):
            self._add_end_member_card(EndMember.from_dict(em_dict), is_dominant=(i == 0))
        self.dopants.clear()
        self._rebuild_dopant_rows()
        self.excess_pct.clear()
        self.reagent_choice.clear()
        self._current_reagent_sig = None
        self._suppress = False
        self._on_card_changed()

    def _add_end_member_card(self, em: EndMember, is_dominant: bool = False):
        card = EndMemberCard(em, is_dominant=is_dominant)
        card.changed.connect(self._on_card_changed)
        card.remove_requested.connect(self._remove_end_member_card)
        self.end_member_cards.append(card)
        self.cards_layout.addWidget(card)

    def _add_blank_end_member(self):
        em = EndMember(name=f"EM{len(self.end_member_cards) + 1}",
                       fraction=0.05, A={}, B={})
        self._add_end_member_card(em, is_dominant=False)
        self._on_card_changed()

    def _remove_end_member_card(self, card: EndMemberCard):
        if getattr(card, "is_dominant", False):
            return
        if len(self.end_member_cards) <= 1:
            return
        if card in self.end_member_cards:
            self.end_member_cards.remove(card)
            card.setParent(None)
            card.deleteLater()
        self._on_card_changed()

    def _on_card_changed(self):
        """Card-level changes trigger dominant-fraction recompute then full recalc."""
        self._update_dominant_fraction()
        self._recalculate()

    def _update_dominant_fraction(self):
        """If the first card is marked dominant, set its fraction to 1 − Σ(others), clamped at 0."""
        if not self.end_member_cards:
            return
        first = self.end_member_cards[0]
        if not getattr(first, "is_dominant", False):
            return
        others_sum = sum(c.frac_spin.value() for c in self.end_member_cards[1:])
        new_value = max(0.0, 1.0 - others_sum)
        if abs(first.frac_spin.value() - new_value) > 1e-9:
            first.frac_spin.blockSignals(True)
            first.frac_spin.setValue(new_value)
            first.frac_spin.blockSignals(False)
        self._recalculate()

    def _base_cations_in_use(self) -> list[str]:
        """All cations currently used in any end-member's A or B site."""
        seen: set[str] = set()
        for card in self.end_member_cards:
            em = card.to_endmember()
            seen.update(em.A.keys())
            seen.update(em.B.keys())
        return sorted(seen)

    def _add_dopant_via_picker(self):
        # Exclude base composition cations and any already-added dopant.
        excluded = self._base_cations_in_use() + [d.cation for d in self.dopants]
        picked = PeriodicTablePicker.pick(self, site=None, excluded=excluded)
        if picked:
            self.dopants.append(DopantEntry(cation=picked, amount=1.0, unit=UNIT_MOL_PCT))
            self._rebuild_dopant_rows()
            self._recalculate()

    def _rebuild_dopant_rows(self):
        """Clear and re-render the dopant rows from self.dopants."""
        # Tear down previous widgets
        while self.dopant_rows_layout.count():
            item = self.dopant_rows_layout.takeAt(0)
            w = item.widget()
            if w is not None:
                w.setParent(None)
                w.deleteLater()
        self._dopant_row_widgets.clear()

        for i, dop in enumerate(self.dopants):
            row = QWidget()
            h = QHBoxLayout(row)
            h.setContentsMargins(0, 0, 0, 0)

            cation_lbl = QLabel(f"<b>{dop.cation}</b>")
            cation_lbl.setMinimumWidth(70)
            h.addWidget(cation_lbl)

            amount_spin = make_spin(dop.amount, step=0.1, decimals=4,
                                    minimum=0.0, maximum=100.0)
            amount_spin.setMinimumWidth(110)
            amount_spin.setMaximumWidth(110)
            h.addWidget(amount_spin)

            unit_combo = NoWheelComboBox()
            unit_combo.addItem("mol %", userData=UNIT_MOL_PCT)
            unit_combo.addItem("wt %", userData=UNIT_WT_PCT)
            unit_combo.setCurrentIndex(0 if dop.unit == UNIT_MOL_PCT else 1)
            unit_combo.setMinimumWidth(110)
            unit_combo.setMaximumWidth(110)
            h.addWidget(unit_combo)

            remove_btn = QPushButton("✕")
            remove_btn.setFixedWidth(32)
            h.addWidget(remove_btn)
            h.addStretch(1)

            amount_spin.valueChanged.connect(
                lambda v, idx=i: self._on_dopant_amount_changed(idx, v))
            unit_combo.currentIndexChanged.connect(
                lambda _i, idx=i, c=unit_combo: self._on_dopant_unit_changed(idx, c))
            remove_btn.clicked.connect(lambda _c=False, idx=i: self._remove_dopant(idx))

            self.dopant_rows_layout.addWidget(row)
            self._dopant_row_widgets.append(row)

    def _on_dopant_amount_changed(self, idx: int, value: float):
        if 0 <= idx < len(self.dopants):
            self.dopants[idx].amount = value
            self._recalculate()

    def _on_dopant_unit_changed(self, idx: int, combo: QComboBox):
        if 0 <= idx < len(self.dopants):
            unit = combo.currentData() or UNIT_MOL_PCT
            self.dopants[idx].unit = unit
            self._recalculate()

    def _remove_dopant(self, idx: int):
        if 0 <= idx < len(self.dopants):
            del self.dopants[idx]
            self._rebuild_dopant_rows()
            self._recalculate()

    # ----- state -----

    def _current_endmembers(self) -> list[EndMember]:
        return [c.to_endmember() for c in self.end_member_cards]

    def _current_dopants(self) -> list[DopantEntry]:
        return [DopantEntry(cation=d.cation, amount=d.amount, unit=d.unit)
                for d in self.dopants]

    def _reagents_by_element(self) -> dict[str, Reagent]:
        out: dict[str, Reagent] = {}
        for el, name in self.reagent_choice.items():
            r = next((r for r in self.reagents if r.name == name), None)
            if r is not None:
                out[el] = r
        return out

    # ----- recalc -----

    def _recalculate(self):
        if self._suppress:
            return
        end_members = self._current_endmembers()
        dopants = self._current_dopants()

        em_sum = sum(em.fraction for em in end_members)
        a_coeffs, b_coeffs = build_composition_coeffs(end_members)
        a_sum = sum(a_coeffs.values())
        b_sum = sum(b_coeffs.values())

        def fmt(label, val):
            ok = abs(val - 1.0) < 1e-3
            color = "#047857" if ok else "#b45309"
            return f"<span style='color:{color};font-weight:600'>Σ {label} = {val:.4f}</span>"
        self.em_sum_label.setText(fmt("end-members", em_sum))
        self.a_sum_label.setText(fmt("A-site", a_sum))
        self.b_sum_label.setText(fmt("B-site", b_sum))

        elements = sorted([el for el in set(list(a_coeffs.keys()) + list(b_coeffs.keys()))
                           if (a_coeffs.get(el, 0) + b_coeffs.get(el, 0)) > 1e-12])

        if not elements:
            self.composition_preview.setText(
                "<i style='color:#94a3b8'>Add A-site and B-site cations to begin, "
                "or load a preset above.</i>"
            )
        else:
            self.composition_preview.setText(prettify_formula(format_composition(end_members, dopants)))

        # Dopant cations also need a reagent picked so their masses can be computed.
        dopant_elements = sorted({d.cation for d in dopants if d.cation})
        for el in elements + dopant_elements:
            current = self.reagent_choice.get(el)
            valid = any(r.name == current and r.element == el for r in self.reagents)
            if not valid:
                cand = next((r for r in self.reagents if r.element == el), None)
                self.reagent_choice[el] = cand.name if cand else ""

        sig = tuple((el, a_coeffs.get(el, 0) > 0, b_coeffs.get(el, 0) > 0) for el in elements)
        if sig != self._current_reagent_sig:
            self._rebuild_reagent_rows(elements, a_coeffs, b_coeffs)
            self._current_reagent_sig = sig

        result = calculate(
            batch_size_g=self.batch_spin.value(),
            end_members=end_members,
            dopants=dopants,
            reagents_by_element=self._reagents_by_element(),
            excess_percent=self.excess_pct,
            apply_purity=self.purity_cb.isChecked(),
        )

        if not elements:
            self.mw_label.setText("Total MW: —")
            self.moles_label.setText("Moles: —")
            self.mass_excess_label.setText("Base reagent: —")
            self.mass_weigh_label.setText("Grand total: —")
            self.results_table.setRowCount(0)
            self.dopants_table.setRowCount(0)
            self.dopants_table.setVisible(False)
            self.dopants_section_label.setVisible(False)
            self.oxygen_label.setText("")
        else:
            base_label = "Base w/ purity" if self.purity_cb.isChecked() else "Base reagent"
            self.mw_label.setText(f"Total MW: {result.total_mw:.4f} g/mol")
            self.moles_label.setText(f"Moles: {result.moles:.4f}")
            self.mass_excess_label.setText(f"{base_label}: {result.total_mass_to_weigh:.3f} g")
            if result.dopant_rows:
                self.mass_weigh_label.setText(
                    f"Grand total: {result.grand_total_to_weigh:.3f} g "
                    f"(+{result.total_dopant_to_weigh:.3f} g dopant)"
                )
            else:
                self.mass_weigh_label.setText(f"Grand total: {result.total_mass_to_weigh:.3f} g")

            self.results_table.setRowCount(len(result.rows))
            for i, r in enumerate(result.rows):
                site = ("A" if r.on_a_site else "") + ("B" if r.on_b_site else "")
                vals = [
                    r.element, site, f"{r.atomic_weight:.3f}",
                    f"{r.coeff:.4f}", f"{r.excess_pct:.4g}",
                    prettify_formula(r.reagent.name) if r.reagent else "—",
                    f"{r.reagent.mw:.2f}" if r.reagent else "",
                    f"{r.reagent.purity:.2f}" if r.reagent else "",
                    f"{r.mass_to_weigh:.4f}",
                ]
                for j, v in enumerate(vals):
                    self.results_table.setItem(i, j, QTableWidgetItem(str(v)))

            # Add-on dopants table
            if result.dopant_rows:
                self.dopants_table.setVisible(True)
                self.dopants_section_label.setVisible(True)
                self.dopants_table.setRowCount(len(result.dopant_rows))
                for i, d in enumerate(result.dopant_rows):
                    unit_label = "wt %" if d.unit == UNIT_WT_PCT else "mol %"
                    vals = [
                        d.cation,
                        f"{d.amount:.4g} {unit_label}",
                        f"{d.atomic_weight:.3f}",
                        prettify_formula(d.reagent.name) if d.reagent else "—",
                        f"{d.reagent.mw:.2f}" if d.reagent else "",
                        f"{d.reagent.purity:.2f}" if d.reagent else "",
                        f"{d.mass_to_weigh:.4f}",
                    ]
                    for j, v in enumerate(vals):
                        self.dopants_table.setItem(i, j, QTableWidgetItem(str(v)))
            else:
                self.dopants_table.setRowCount(0)
                self.dopants_table.setVisible(False)
                self.dopants_section_label.setVisible(False)

            self.oxygen_label.setText(
                f"Oxygen: {result.oxygen_coeff:.3f} atoms × 15.999 = "
                f"{result.oxygen_mass:.4f} g/mol (supplied by reagents)."
            )

        self._last_result = result
        self._last_end_members = end_members
        self._last_dopants = dopants

    def _rebuild_reagent_rows(self, elements: list[str],
                              a_coeffs: dict[str, float], b_coeffs: dict[str, float]):
        while self.reagents_layout.count():
            item = self.reagents_layout.takeAt(0)
            w = item.widget()
            if w is not None:
                w.setParent(None)
                w.deleteLater()

        for el in elements:
            row = QWidget()
            h = QHBoxLayout(row)
            h.setContentsMargins(0, 0, 0, 0)
            on_a = a_coeffs.get(el, 0) > 0
            on_b = b_coeffs.get(el, 0) > 0
            site = ("A" if on_a else "") + ("B" if on_b else "")

            lbl = QLabel(f"<b>{el}</b> <span style='color:#64748b;font-size:12px'>[{site}]</span>")
            lbl.setMinimumWidth(80)
            h.addWidget(lbl)

            combo = NoWheelComboBox()
            candidates = [r for r in self.reagents if r.element == el]
            if candidates:
                for r in candidates:
                    combo.addItem(f"{prettify_formula(r.name)} — MW {r.mw:.2f}, {r.purity}%", userData=r.name)
                current = self.reagent_choice.get(el, candidates[0].name)
                idx = next((i for i in range(combo.count())
                            if combo.itemData(i) == current), 0)
                combo.setCurrentIndex(idx)
            else:
                combo.addItem(f"⚠ No reagent for {el}")
                combo.setEnabled(False)
            combo.setMinimumWidth(280)
            h.addWidget(combo)

            h.addWidget(QLabel("Excess %:"))
            spin = make_spin(self.excess_pct.get(el, 0.0), step=0.5, decimals=4,
                             minimum=-100.0, maximum=1000.0)
            spin.setMinimumWidth(100)
            h.addWidget(spin)
            h.addStretch(1)

            combo.currentIndexChanged.connect(
                lambda _i, e=el, c=combo: self._on_reagent_changed(e, c))
            spin.valueChanged.connect(
                lambda v, e=el: self._on_excess_changed(e, v))

            self.reagents_layout.addWidget(row)

    def _on_reagent_changed(self, el: str, combo: QComboBox):
        name = combo.currentData()
        if name:
            self.reagent_choice[el] = name
            self._recalculate()

    def _on_excess_changed(self, el: str, value: float):
        self.excess_pct[el] = value
        self._recalculate()

    # ----- export -----

    def _export_csv(self):
        if not hasattr(self, "_last_result"):
            return
        name = self.recipe_name.text().strip() or "recipe"
        ts = int(datetime.now().timestamp())
        default_path = str(Path.home() / f"{name}_{ts}.csv")
        path, _ = QFileDialog.getSaveFileName(self, "Export CSV", default_path, "CSV (*.csv)")
        if not path:
            return
        try:
            with open(path, "w", newline="", encoding="utf-8") as f:
                w = csv.writer(f)
                r = self._last_result
                w.writerow(["Recipe Export", datetime.now().isoformat()])
                w.writerow(["Name", name])
                w.writerow(["Composition", format_composition(self._last_end_members, self._last_dopants)])
                w.writerow(["Batch size (g)", self.batch_spin.value()])
                w.writerow(["Purity correction", "applied" if self.purity_cb.isChecked() else "not applied"])
                w.writerow(["Total MW", f"{r.total_mw:.4f}"])
                w.writerow(["Moles", f"{r.moles:.6f}"])
                w.writerow([])
                last_col = "Purity-corrected (g)" if self.purity_cb.isChecked() else "To weigh (g)"
                w.writerow(["Element", "Site", "At wt", "Coeff", "Excess (mol %)",
                            "Reagent", "Reagent MW", "Purity %",
                            "Mass target (g)", "Mass + excess (g)", last_col])
                for row in r.rows:
                    site = ("A" if row.on_a_site else "") + ("B" if row.on_b_site else "")
                    w.writerow([row.element, site, row.atomic_weight, f"{row.coeff:.6f}",
                                row.excess_pct,
                                row.reagent.name if row.reagent else "",
                                row.reagent.mw if row.reagent else "",
                                row.reagent.purity if row.reagent else "",
                                f"{row.mass_target:.4f}", f"{row.mass_with_excess:.4f}",
                                f"{row.mass_to_weigh:.4f}"])
                if r.dopant_rows:
                    w.writerow([])
                    w.writerow(["Add-on dopants (in addition to the base batch)"])
                    w.writerow(["Cation", "Amount", "Unit", "At wt",
                                "Reagent", "Reagent MW", "Purity %",
                                "Mass (g)", last_col])
                    for d in r.dopant_rows:
                        unit_label = "wt %" if d.unit == UNIT_WT_PCT else "mol %"
                        w.writerow([d.cation, f"{d.amount:.6f}", unit_label,
                                    d.atomic_weight,
                                    d.reagent.name if d.reagent else "",
                                    d.reagent.mw if d.reagent else "",
                                    d.reagent.purity if d.reagent else "",
                                    f"{d.mass:.4f}", f"{d.mass_to_weigh:.4f}"])
                    w.writerow([])
                    w.writerow(["Base total to weigh (g)", f"{r.total_mass_to_weigh:.4f}"])
                    w.writerow(["Dopant total to weigh (g)", f"{r.total_dopant_to_weigh:.4f}"])
                    w.writerow(["Grand total to weigh (g)", f"{r.grand_total_to_weigh:.4f}"])
            QMessageBox.information(self, "Exported", f"Saved to:\n{path}")
        except Exception as e:
            QMessageBox.critical(self, "Export failed", str(e))

    # ----- Raw materials tab -----

    def _build_raw_materials_tab(self) -> QWidget:
        tab = QWidget()
        v = QVBoxLayout(tab)
        v.setContentsMargins(16, 16, 16, 16)
        v.setSpacing(10)

        title = QLabel("Raw materials database")
        f = QFont(); f.setPointSize(16); f.setBold(True); title.setFont(f)
        v.addWidget(title)
        sub = QLabel("Add or edit the precursors your lab uses. Changes update the Calculator tab live.")
        sub.setStyleSheet("color: #64748b;")
        v.addWidget(sub)

        toolbar = QHBoxLayout()
        add_btn = QPushButton("+ Add reagent")
        remove_btn = QPushButton("− Remove selected row(s)")
        reset_btn = QPushButton("↺ Reset to defaults")
        add_btn.clicked.connect(self._add_reagent_row)
        remove_btn.clicked.connect(self._remove_selected_reagents)
        reset_btn.clicked.connect(self._reset_reagents_to_defaults)
        toolbar.addWidget(add_btn)
        toolbar.addWidget(remove_btn)
        toolbar.addWidget(reset_btn)
        toolbar.addStretch(1)
        v.addLayout(toolbar)

        self.reagents_table = QTableWidget(0, 6)
        self.reagents_table.setHorizontalHeaderLabels(
            ["Name", "Element", "Atoms", "MW (g/mol)", "Purity %", "Notes"])
        self.reagents_table.horizontalHeader().setSectionResizeMode(QHeaderView.Interactive)
        self.reagents_table.horizontalHeader().setStretchLastSection(True)
        self.reagents_table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.reagents_table.verticalHeader().setVisible(False)
        self._populate_reagents_table()
        self.reagents_table.cellChanged.connect(self._on_reagent_cell_changed)
        v.addWidget(self.reagents_table, stretch=1)

        hint = QLabel(
            "<b>Element</b> = which target element this reagent supplies. "
            "<b>Atoms</b> = atoms of that element per reagent formula unit "
            "(e.g. Na₂CO₃ → element=Na, atoms=2)."
        )
        hint.setStyleSheet("color: #64748b; font-size: 12px;")
        hint.setWordWrap(True)
        v.addWidget(hint)

        io = QHBoxLayout()
        io.addWidget(QLabel("Backup / share:"))
        import_btn = QPushButton("⬆ Import JSON…")
        export_btn = QPushButton("⬇ Export JSON…")
        import_btn.clicked.connect(self._import_reagents_json)
        export_btn.clicked.connect(self._export_reagents_json)
        io.addWidget(import_btn)
        io.addWidget(export_btn)
        io.addStretch(1)
        v.addLayout(io)

        return tab

    def _populate_reagents_table(self):
        self.reagents_table.blockSignals(True)
        self.reagents_table.setRowCount(len(self.reagents))
        for i, r in enumerate(self.reagents):
            self.reagents_table.setItem(i, 0, QTableWidgetItem(prettify_formula(r.name)))
            self.reagents_table.setItem(i, 1, QTableWidgetItem(r.element))
            self.reagents_table.setItem(i, 2, QTableWidgetItem(str(r.atoms)))
            self.reagents_table.setItem(i, 3, QTableWidgetItem(f"{r.mw:.4f}"))
            self.reagents_table.setItem(i, 4, QTableWidgetItem(f"{r.purity:.3f}"))
            self.reagents_table.setItem(i, 5, QTableWidgetItem(r.notes))
        self.reagents_table.blockSignals(False)

    def _add_reagent_row(self):
        self.reagents.append(Reagent("New", "X", 1, 0.0, 99.9, ""))
        self._populate_reagents_table()
        self._invalidate_reagent_cache_and_recalc()

    def _remove_selected_reagents(self):
        selected = sorted({i.row() for i in self.reagents_table.selectedIndexes()}, reverse=True)
        if not selected:
            return
        for row in selected:
            if 0 <= row < len(self.reagents):
                del self.reagents[row]
        self._populate_reagents_table()
        self._invalidate_reagent_cache_and_recalc()

    def _reset_reagents_to_defaults(self):
        confirm = QMessageBox.question(
            self, "Reset reagents",
            "Replace the current reagent database with the built-in defaults? "
            "Any custom entries will be lost.",
            QMessageBox.Yes | QMessageBox.No,
        )
        if confirm != QMessageBox.Yes:
            return
        self.reagents = [Reagent.from_dict(r) for r in DEFAULT_REAGENTS]
        self._populate_reagents_table()
        self._invalidate_reagent_cache_and_recalc()

    def _on_reagent_cell_changed(self, row: int, col: int):
        if row >= len(self.reagents):
            return
        r = self.reagents[row]
        item = self.reagents_table.item(row, col)
        val = item.text() if item else ""
        try:
            if col == 0:
                r.name = normalize_formula(val).strip()
            elif col == 1:
                r.element = normalize_formula(val).strip()
            elif col == 2:
                r.atoms = max(1, int(float(val)))
            elif col == 3:
                r.mw = float(val)
            elif col == 4:
                r.purity = float(val)
            elif col == 5:
                r.notes = val
        except (ValueError, TypeError):
            self._populate_reagents_table()
            return
        # Re-render the edited cell with prettified subscripts so the table
        # stays visually consistent regardless of what the user typed.
        if col == 0 and item is not None:
            self.reagents_table.blockSignals(True)
            item.setText(prettify_formula(r.name))
            self.reagents_table.blockSignals(False)
        self._invalidate_reagent_cache_and_recalc()

    def _invalidate_reagent_cache_and_recalc(self):
        self._current_reagent_sig = None
        self._recalculate()

    def _export_reagents_json(self):
        ts = int(datetime.now().timestamp())
        default = str(Path.home() / f"reagents_{ts}.json")
        path, _ = QFileDialog.getSaveFileName(self, "Export reagents JSON", default, "JSON (*.json)")
        if not path:
            return
        try:
            data = [r.to_dict() for r in self.reagents]
            with open(path, "w", encoding="utf-8") as f:
                json.dump(data, f, indent=2)
            QMessageBox.information(self, "Exported", f"Saved {len(data)} reagents to:\n{path}")
        except Exception as e:
            QMessageBox.critical(self, "Export failed", str(e))

    def _import_reagents_json(self):
        path, _ = QFileDialog.getOpenFileName(self, "Import reagents JSON",
                                              str(Path.home()), "JSON (*.json)")
        if not path:
            return
        try:
            with open(path, "r", encoding="utf-8") as f:
                data = json.load(f)
            self.reagents = [Reagent.from_dict(r) for r in data]
            self._populate_reagents_table()
            self._invalidate_reagent_cache_and_recalc()
            QMessageBox.information(self, "Imported", f"Loaded {len(data)} reagents.")
        except Exception as e:
            QMessageBox.critical(self, "Import failed", str(e))

    # ----- Saved recipes tab -----

    def _build_saved_recipes_tab(self) -> QWidget:
        tab = QWidget()
        v = QVBoxLayout(tab)
        v.setContentsMargins(16, 16, 16, 16)
        v.setSpacing(10)

        title = QLabel("Saved recipes")
        f = QFont(); f.setPointSize(16); f.setBold(True); title.setFont(f)
        v.addWidget(title)
        sub = QLabel("Recipes saved from the Calculator tab. Use the JSON tools below to share recipes across the lab.")
        sub.setStyleSheet("color: #64748b;")
        sub.setWordWrap(True)
        v.addWidget(sub)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        inner = QWidget()
        scroll.setWidget(inner)
        self.recipes_layout = QVBoxLayout(inner)
        self.recipes_layout.setSpacing(8)
        self.recipes_layout.setContentsMargins(0, 0, 0, 0)
        v.addWidget(scroll, stretch=1)

        io = QHBoxLayout()
        io.addWidget(QLabel("Backup / share:"))
        import_btn = QPushButton("⬆ Import JSON…")
        export_btn = QPushButton("⬇ Export all JSON…")
        import_btn.clicked.connect(self._import_recipes_json)
        export_btn.clicked.connect(self._export_recipes_json)
        io.addWidget(import_btn)
        io.addWidget(export_btn)
        io.addStretch(1)
        v.addLayout(io)

        self._refresh_saved_recipes_list()
        return tab

    def _refresh_saved_recipes_list(self):
        while self.recipes_layout.count():
            item = self.recipes_layout.takeAt(0)
            w = item.widget()
            if w is not None:
                w.setParent(None)
                w.deleteLater()

        if not self.saved_recipes:
            empty = QLabel("No saved recipes yet. Save one from the Calculator tab.")
            empty.setStyleSheet("color:#94a3b8; padding:20px;")
            empty.setAlignment(Qt.AlignCenter)
            self.recipes_layout.addWidget(empty)
        else:
            for recipe in self.saved_recipes:
                self.recipes_layout.addWidget(self._make_recipe_card(recipe))
        self.recipes_layout.addStretch(1)

    @staticmethod
    def _dopants_from_recipe(recipe: dict) -> list[DopantEntry]:
        """Read dopants from either the new list format or the legacy single-dopant dict."""
        if "dopants" in recipe and isinstance(recipe["dopants"], list):
            return [DopantEntry.from_dict(d) for d in recipe["dopants"]]
        legacy = dopant_entry_from_legacy_dict(recipe.get("dopant") or {})
        return [legacy] if legacy is not None else []

    def _make_recipe_card(self, recipe: dict) -> QGroupBox:
        card = QGroupBox()
        h = QHBoxLayout(card)
        info = QVBoxLayout()
        name_lbl = QLabel(f"<b>{recipe['name']}</b>")
        em_list = [EndMember.from_dict(em) for em in recipe["end_members"]]
        dops = self._dopants_from_recipe(recipe)
        comp = QLabel(prettify_formula(format_composition(em_list, dops)))
        comp.setStyleSheet("font-family:monospace; color:#475569;")
        comp.setWordWrap(True)
        meta = QLabel(f"{recipe['batch_size']:.2f} g · {recipe['timestamp']}")
        meta.setStyleSheet("color:#94a3b8; font-size:12px;")
        info.addWidget(name_lbl)
        info.addWidget(comp)
        info.addWidget(meta)
        h.addLayout(info, stretch=1)

        load_btn = QPushButton("📥 Load")
        delete_btn = QPushButton("🗑 Delete")
        load_btn.clicked.connect(lambda: self._load_saved_recipe(recipe))
        delete_btn.clicked.connect(lambda: self._delete_saved_recipe(recipe["id"]))
        h.addWidget(load_btn)
        h.addWidget(delete_btn)
        return card

    def _save_current_recipe(self):
        name = self.recipe_name.text().strip()
        if not name:
            QMessageBox.warning(self, "Name required", "Please enter a recipe name first.")
            return
        recipe = {
            "id": int(datetime.now().timestamp() * 1000),
            "name": name,
            "timestamp": datetime.now().isoformat(timespec="seconds"),
            "batch_size": self.batch_spin.value(),
            "end_members": [c.to_endmember().to_dict() for c in self.end_member_cards],
            "dopants": [d.to_dict() for d in self._current_dopants()],
            "reagent_choice": dict(self.reagent_choice),
            "excess_pct": dict(self.excess_pct),
            "apply_purity": self.purity_cb.isChecked(),
        }
        self.saved_recipes.insert(0, recipe)
        self._refresh_saved_recipes_list()
        QMessageBox.information(self, "Saved",
                                f"'{name}' saved. See the Saved Recipes tab.")

    def _delete_saved_recipe(self, recipe_id: int):
        self.saved_recipes = [r for r in self.saved_recipes if r["id"] != recipe_id]
        self._refresh_saved_recipes_list()

    def _load_saved_recipe(self, recipe: dict):
        self._suppress = True
        for card in list(self.end_member_cards):
            card.setParent(None)
            card.deleteLater()
        self.end_member_cards.clear()
        for i, em_dict in enumerate(recipe["end_members"]):
            self._add_end_member_card(EndMember.from_dict(em_dict), is_dominant=(i == 0))

        self.dopants = self._dopants_from_recipe(recipe)
        self._rebuild_dopant_rows()

        self.batch_spin.setValue(float(recipe["batch_size"]))
        self.purity_cb.setChecked(bool(recipe.get("apply_purity", True)))
        self.reagent_choice = dict(recipe.get("reagent_choice", {}))
        self.excess_pct = {k: float(v) for k, v in recipe.get("excess_pct", {}).items()}
        self._current_reagent_sig = None
        self._suppress = False
        self._on_card_changed()
        self.tabs.setCurrentIndex(0)

    def _export_recipes_json(self):
        if not self.saved_recipes:
            QMessageBox.information(self, "Nothing to export", "No saved recipes yet.")
            return
        ts = int(datetime.now().timestamp())
        default = str(Path.home() / f"recipes_{ts}.json")
        path, _ = QFileDialog.getSaveFileName(self, "Export recipes JSON", default, "JSON (*.json)")
        if not path:
            return
        try:
            with open(path, "w", encoding="utf-8") as f:
                json.dump(self.saved_recipes, f, indent=2)
            QMessageBox.information(self, "Exported",
                                    f"Saved {len(self.saved_recipes)} recipes to:\n{path}")
        except Exception as e:
            QMessageBox.critical(self, "Export failed", str(e))

    def _import_recipes_json(self):
        path, _ = QFileDialog.getOpenFileName(self, "Import recipes JSON",
                                              str(Path.home()), "JSON (*.json)")
        if not path:
            return
        try:
            with open(path, "r", encoding="utf-8") as f:
                data = json.load(f)
            self.saved_recipes = list(data) + self.saved_recipes
            self._refresh_saved_recipes_list()
            QMessageBox.information(self, "Imported", f"Loaded {len(data)} recipes.")
        except Exception as e:
            QMessageBox.critical(self, "Import failed", str(e))


def main():
    app = QApplication(sys.argv)
    base_font = app.font()
    base_font.setPointSize(14)
    app.setFont(base_font)
    win = MainWindow()
    win.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
