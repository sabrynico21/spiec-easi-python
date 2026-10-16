# Packaging

Notes on shipping `spieceasi` so it can be driven by the Electron application.

## Wheel / source distribution

```bash
python -m build           # produces dist/spieceasi-*.whl and *.tar.gz
```

`pyproject.toml` declares `numpy` as the only runtime dependency, so the wheel
is platform-independent (`py3-none-any`) as long as the port stays pure NumPy.

## Standalone executable (PyInstaller)

The Electron backend can bundle the analysis process as a frozen executable.
Because the port is pure NumPy (no compiled extension in this repository), the
usual PyInstaller pitfalls are limited to:

- ensure NumPy's binaries are collected (PyInstaller handles this via the
  bundled hooks);
- if StARS parallelism is implemented with `multiprocessing`, call
  `multiprocessing.freeze_support()` at the entry point and prefer a
  **thread** backend in frozen builds (process pools re-execute the binary);
- install `setuptools` in the build environment so NumPy hooks resolve.

Example:

```bash
pyinstaller --name spieceasi-cli -F \
  --collect-all numpy \
  -m spieceasi
```

## CLI contract for Electron

- stdin/file in, stdout out, single JSON object, no extra prints on stdout.
- Warnings and progress on stderr only.
- Non-zero exit code on error; the error message is written to stderr.
- Sparse matrices serialised as edge lists with 0-based indices.

See `README.md` for the request/response schema.
