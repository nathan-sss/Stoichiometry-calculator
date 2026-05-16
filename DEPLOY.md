# Publishing as a free website via GitHub Pages + stlite

stlite compiles Streamlit to WebAssembly so the whole app runs inside the
visitor's browser. No Python server, no rent — GitHub Pages just serves
the static files, the visitor's browser does the rest.

Result: a URL like `https://<your-username>.github.io/<repo-name>/`
that anyone in the world can open.

## One-time setup

1. **Push the project to GitHub** (public repo). The repo needs to contain
   at minimum: `app.py`, `calculator.py`, `data.py`, `index.html`,
   `.nojekyll`.

2. On the GitHub repo page: **Settings → Pages**.
   - Under "Build and deployment", set **Source = Deploy from a branch**.
   - Set **Branch = main**, folder = **/ (root)**.
   - Click **Save**.

3. Wait ~1 minute. GitHub shows the live URL at the top of the Pages
   settings page once the build finishes.

## Updating the live site

Every `git push` to the deployment branch triggers a redeploy automatically:

```bash
git add -A
git commit -m "your description of changes"
git push
```

GitHub rebuilds within ~30 seconds. Visitors get the new version on their
next page load (no caching shenanigans — stlite caches Python, but your
`app.py` is fetched fresh each load).

## How it works (and what to expect)

The first time a visitor opens the site:

- Their browser downloads `index.html` (tiny, ~2 KB).
- Then downloads stlite (~10 MB) and pandas (~10 MB) from a CDN.
- Pyodide initialises a Python interpreter in WebAssembly inside the page.
- Finally, your `app.py` runs and the UI appears.

Total first-load time: **20–40 seconds** depending on network. Subsequent
visits are near-instant — the browser caches everything.

## Caveats

- **No file uploads to the server**, because there is no server. The JSON
  import buttons (`st.file_uploader`) read files locally in the browser
  via `FileReader`, so they work, but the user has to pick a file from
  *their* device.
- **No backend persistence**. Anything a user types is in their own
  browser session only — nothing is saved between visits or shared with
  other visitors.
- **Mobile**: the page works, but the periodic-table picker is still
  cramped on phones (this is a Streamlit layout limit, not a stlite limit).
- **stlite version**: `index.html` pins `@stlite/browser@0.83.0`. If that
  CDN URL ever 404s (a yanked release), bump to the latest version listed
  at <https://www.npmjs.com/package/@stlite/browser>.

## Troubleshooting

**Page shows the loading spinner forever**

- Open browser devtools → Console. Look for red errors.
- Most common cause: a typo in `index.html` or one of the `app-file` URLs
  is wrong.
- Second-most common: stlite version mismatch. Try
  `https://cdn.jsdelivr.net/npm/@stlite/browser/build/stlite.js` (no
  version pin = latest).

**ModuleNotFoundError for some package**

- Add the package to `<app-requirements>` in `index.html`.
- Not every PyPI package works in WebAssembly. Pure-Python or
  Pyodide-pre-built packages are safe; anything with native C extensions
  may not work.

**Periodic-table picker doesn't open**

- `@st.dialog` requires Streamlit ≥1.31. stlite 0.83.0 bundles a recent
  Streamlit, so this should work. If you bump to an older stlite version,
  the dialog won't render.
