# Publishing rockmoc to PyPI

Package is configured for modern **PEP 517** installs via hatchling.

- **Author:** shanthoshkv (`shanthoshkv@gmail.com`)
- **Repo:** https://github.com/shanthoshkv/rockmoc
- **Version source:** `src/rockmoc/_version.py` (read by hatch)

## Local install (development)

```bash
cd rockmoc
python -m venv .venv
# Windows: .venv\Scripts\activate
# Unix:    source .venv/bin/activate

pip install -U pip
pip install -e ".[plot,cli,dev]"
# If editable metadata fails on older pip, use:
#   pip install -e ".[plot,cli,dev]" --config-settings editable_mode=compat
python -c "import rockmoc; print(rockmoc.__version__)"
rockmoc --version
pytest -q
```

## Build distributions

```bash
pip install -e ".[dev]"   # includes build + twine
python -m build
# → dist/rockmoc-<version>.tar.gz
# → dist/rockmoc-<version>-py3-none-any.whl
```

Smoke-test the wheel in a **clean** venv:

```bash
python -m venv /tmp/rockmoc-smoke
# activate it
pip install dist/rockmoc-*-py3-none-any.whl
python -c "import rockmoc; from rockmoc import design_mln, DesignSpec, PerfectGas; print(rockmoc.__version__)"
pip install "rockmoc[cli]@file://$(pwd)/dist/..."   # optional: reinstall with extras from path
```

Or simply:

```bash
pip install dist/rockmoc-*-py3-none-any.whl
pip install typer matplotlib   # if you want CLI/plots without extras reinstall
```

Check metadata:

```bash
twine check dist/*
```

## TestPyPI (recommended first)

1. Create account + API token: https://test.pypi.org/manage/account/  
2. Upload:

```bash
twine upload --repository testpypi dist/*
```

3. Install from TestPyPI:

```bash
pip install -i https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple/ rockmoc
```

(`extra-index-url` supplies numpy/scipy from real PyPI.)

## PyPI (production)

1. Create account + API token: https://pypi.org/manage/account/  
2. Confirm the name `rockmoc` is free (or change `name` in `pyproject.toml`).  
3. Upload:

```bash
twine upload dist/*
```

4. Verify:

```bash
pip install rockmoc
pip install "rockmoc[plot,cli]"
python -c "import rockmoc; print(rockmoc.__version__)"
rockmoc --version
```

## Bump version

Keep these two in sync:

1. `pyproject.toml` → `[project] version = "…"`
2. `src/rockmoc/_version.py` → `__version__ = "…"`

```python
__version__ = "0.6.1"
```

Then rebuild.

## GitHub release workflow

1. Commit and push.  
2. Tag:

```bash
git tag v0.6.0
git push origin v0.6.0
```

3. Optional: store `PYPI_API_TOKEN` in GitHub repo secrets and use  
   `.github/workflows/publish.yml` (uploads on `v*` tags).

## What goes in the wheel vs sdist

| Artifact | Contents |
|----------|----------|
| **Wheel** | `rockmoc` package under `site-packages` only |
| **Sdist** | source + tests, examples, docs (for packagers / `pip install` from sdist) |

Not published: `website/`, `refs/*.pdf` (large), generated graphs (gitignored).

## Extras users can install

```bash
pip install rockmoc                # core
pip install "rockmoc[plot]"        # matplotlib
pip install "rockmoc[cli]"         # typer CLI
pip install "rockmoc[cea]"         # RocketCEA
pip install "rockmoc[all]"         # everything
```
