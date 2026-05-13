# PerovskiteCalc.spec — PyInstaller build configuration
#
# Build with:
#   pyinstaller PerovskiteCalc.spec
#
# Produces a single-file, windowed executable in dist/:
#   - Windows: dist/PerovskiteCalc.exe
#   - macOS:   dist/PerovskiteCalc (binary) + dist/PerovskiteCalc.app (bundle)
#
# Drop an `icon.ico` (Windows) or `icon.icns` (macOS) into this directory
# to embed a custom application icon. The spec auto-detects whichever
# format matches the build platform.

import os
import sys

block_cipher = None


def _pick_icon():
    """Return the path to a platform-appropriate icon, or None."""
    if sys.platform == "darwin" and os.path.exists("icon.icns"):
        return "icon.icns"
    if os.path.exists("icon.ico"):
        return "icon.ico"
    return None


a = Analysis(
    ['main.py'],
    pathex=[],
    binaries=[],
    datas=[],
    hiddenimports=[],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    # Streamlit / pandas / pytest are only used by app.py and tests.
    # Excluding them shrinks the bundle by ~80 MB.
    excludes=[
        'streamlit', 'pandas', 'pytest', 'numpy', 'altair', 'pyarrow',
        'matplotlib', 'IPython', 'jupyter',
    ],
    win_no_prefer_redirects=False,
    win_private_assemblies=False,
    cipher=block_cipher,
    noarchive=False,
)
pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.zipfiles,
    a.datas,
    [],
    name='PerovskiteCalc',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    upx_exclude=[],
    runtime_tmpdir=None,
    console=False,            # --windowed: no terminal window on launch
    disable_windowed_traceback=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
    icon=_pick_icon(),
)
