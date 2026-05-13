# Building PerovskiteCalc.exe on Windows

This is the procedure for producing the standalone `PerovskiteCalc.exe`.
You run these steps **inside your Windows VM** (or any Windows PC) — building
the `.exe` on macOS produces a Mac binary, not a Windows one.

## Prerequisites

- Windows 10 or 11
- Python **3.10 or newer** ([download from python.org](https://www.python.org/downloads/windows/))
  - During install, tick **"Add python.exe to PATH"**.

## One-time setup

Open PowerShell or Command Prompt in the project folder, then:

```powershell
python -m venv .venv
.\.venv\Scripts\activate
pip install -r requirements.txt
```

This installs PySide6, PyInstaller, and pytest into a private virtual
environment (the `.venv` folder). The Streamlit web-version dependencies
are kept in `requirements-web.txt` and are **not** needed on Windows.

**If PowerShell blocks `activate`** with an "execution policy" error, run
this once to allow local scripts for your user (safe — still blocks
downloaded unsigned scripts):

```powershell
Set-ExecutionPolicy -ExecutionPolicy RemoteSigned -Scope CurrentUser
```

## (Optional) Add a custom icon

Place a Windows icon file named **`icon.ico`** in the project folder
(next to `main.py`). The `.spec` file picks it up automatically.

No icon? Skip this step — PyInstaller falls back to a default icon.

## Build the .exe

With the venv active:

```powershell
pyinstaller PerovskiteCalc.spec --noconfirm --clean
```

The build takes 1–3 minutes. When it finishes, the executable is at:

```
dist\PerovskiteCalc.exe
```

That single file is what you distribute. Double-click to launch.

## Smoke-test before distributing

Always run the produced `.exe` once on your build machine before sending
it to colleagues:

```powershell
.\dist\PerovskiteCalc.exe
```

Try:
- Load the **NBT-BT** preset → confirm Total MW ≈ 213 g/mol.
- Click **+ Add cation** on an A-site → confirm the periodic-table picker opens.
- Save a recipe → confirm it appears on the **Saved recipes** tab.
- Export CSV → confirm the file is written.

## Distributing

The `.exe` is fully self-contained — no Python required on the recipient's
machine. Approximate size: ~70 MB on Windows.

**Caveats:**
- First launch takes 3–5 seconds (PyInstaller extracts to a temp folder).
  Subsequent launches are faster once Windows caches the file.
- Some antivirus tools flag freshly-built PyInstaller binaries as suspicious
  (false positive due to the self-extraction pattern). If a recipient sees a
  warning, they can either whitelist the file or you can sign the binary
  with a code-signing certificate.

## Troubleshooting

**`PySide6` not found during build**
- The active Python isn't the venv's. Re-run `.\.venv\Scripts\activate`.

**`Failed to load the Qt platform plugin "windows"` at runtime**
- Usually means antivirus quarantined a Qt DLL inside the extracted bundle.
  Whitelist `dist\PerovskiteCalc.exe` and rebuild.

**Bundle is huge (>200 MB)**
- Check that `excludes` in `PerovskiteCalc.spec` still lists `streamlit`,
  `pandas`, `numpy`, etc. — they're large and not needed for `main.py`.

**Need to debug a crash on launch**
- Edit `PerovskiteCalc.spec`, set `console=True`, rebuild. A terminal window
  will open alongside the GUI and print any Python tracebacks. Set back to
  `False` before distributing.

## Rebuilding after code changes

Just re-run the build command — `--clean` ensures stale cache is removed:

```powershell
pyinstaller PerovskiteCalc.spec --noconfirm --clean
```

No need to recreate the venv unless you bump dependencies in `requirements.txt`.
