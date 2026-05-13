# Perovskite Recipe Calculator

A Streamlit app for calculating stoichiometric raw material additions for ABO₃ perovskite ceramics — NBT, KNN, BFBT, PZT, BiFeO₃, BaTiO₃, SrTiO₃, and any user-defined composition.

Supports:
- **Solid solutions** (multi-end-member compositions like 0.7 BiFeO₃ – 0.3 BaTiO₃)
- **A-site and B-site doping** with periodic-table cation picker
- **Excess additions** to compensate volatilization (e.g. +3% Na for NBT)
- **Optional purity correction**
- **Custom reagent database** (editable, with JSON import/export)
- **Saved recipes** (JSON import/export for sharing across the lab)
- **CSV export** of any calculation for lab records

---

## Project structure

```
.
├── app.py              # Streamlit UI (run this)
├── calculator.py       # Pure calculation engine (no Streamlit dependency)
├── data.py             # Atomic weights, default reagents, presets, periodic table layout
├── test_calculator.py  # Unit tests
├── requirements.txt
└── README.md
```

The calculation logic lives entirely in `calculator.py` and is independent of the UI. You can import and use it from a Jupyter notebook or other Python scripts.

---

## Setup

### 1. Clone or download the project

Drop these files into a folder (e.g. `perovskite-calc/`).

### 2. Create a virtual environment (recommended)

```bash
python -m venv .venv
source .venv/bin/activate     # macOS/Linux
.venv\Scripts\activate         # Windows
```

### 3. Install dependencies

```bash
pip install -r requirements.txt
```

### 4. Run the app

```bash
streamlit run app.py
```

This will open the app in your browser at `http://localhost:8501`.

---

## Running tests

```bash
pytest test_calculator.py -v
```

You should see all tests pass. The tests verify NBT and BaTiO₃ stoichiometry, dopant substitution, excess additions, purity correction, and serialization round-trips.

---

## Using the calculation engine from Python

```python
from calculator import EndMember, Dopant, Reagent, calculate

# Define the composition: 0.94 NBT - 0.06 BT
end_members = [
    EndMember("NBT", 0.94, A={"Na": 0.5, "Bi": 0.5}, B={"Ti": 1.0}),
    EndMember("BT",  0.06, A={"Ba": 1.0},            B={"Ti": 1.0}),
]

# Optional A-site dopant: 1.5 mol% Ca
dopant = Dopant(enabled=True, site="A", cation="Ca", level=0.015)

# Which reagent supplies each element
reagents = {
    "Na": Reagent("Na2CO3", "Na", 2, 105.99, 99.9, ""),
    "Bi": Reagent("Bi2O3",  "Bi", 2, 465.96, 99.9, ""),
    "Ba": Reagent("BaCO3",  "Ba", 1, 197.34, 99.9, ""),
    "Ti": Reagent("TiO2",   "Ti", 1,  79.87, 99.9, ""),
    "Ca": Reagent("CaCO3",  "Ca", 1, 100.09, 99.9, ""),
}

# Run the calculation
result = calculate(
    batch_size_g=50.0,
    end_members=end_members,
    dopant=dopant,
    reagents_by_element=reagents,
    excess_percent={"Na": 3.0, "Bi": 3.0},   # +3% excess for volatile cations
    apply_purity=True,
)

print(f"Total MW: {result.total_mw:.4f} g/mol")
print(f"Moles in batch: {result.moles:.6f}")
for row in result.rows:
    print(f"  {row.element}: weigh {row.mass_to_weigh:.4f} g of {row.reagent.name}")
```

---

## Extending the app

- **Add reagents permanently**: edit `DEFAULT_REAGENTS` in `data.py`.
- **Add presets**: append to the `PRESETS` list in `data.py`.
- **Change A/B-site cation categories**: edit `A_SITE_CATIONS` and `B_SITE_CATIONS` in `data.py` — these only affect periodic-table color-coding, not the math.
- **Add new calculation features**: modify `calculator.py` and add corresponding tests in `test_calculator.py`.

---

## Notes on the calculation

- **Batch size** = mass of the *final fired ceramic product*. The mass of starting precursors you weigh out is usually larger because carbonates lose CO₂ during calcination (loss on ignition).
- **Excess additions** are applied per element as a percentage of that element's moles in the formula unit.
- **Purity correction**, when enabled, divides the reagent mass by `purity / 100`. With 99.9% purity reagents this adds only ~0.1% — negligible compared to LOI but useful for lower-purity precursors.
- **Oxygen** is assumed to come from the reagents (oxides, carbonates) and contributes 3 atoms per ABO₃ formula unit, scaled by the sum of end-member fractions.