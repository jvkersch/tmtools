# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

TM-Tools is a Python package providing bindings for the TM-align algorithm for protein structure comparison. The package consists of:

- **Python interface** (`tmtools/`): High-level Python API with NumPy integration
- **C++ extension** (`src/`): pybind11 bindings around the original TM-align C++ code
- **Core algorithm** (`src/extern/TMalign-modified.cpp`): Modified version of the original TM-align implementation

## Architecture

- `tmtools/__init__.py`: Main entry point exposing `tm_align`, `print_version`, and `TMResult`
- `tmtools/io.py`: BioPython-based utilities for reading PDB/mmCIF files (optional dependency)
- `tmtools/helpers.py`: Utility functions for coordinate manipulation
- `tmtools/testing.py`: Test data utilities for accessing bundled PDB files
- `tmtools/tests/`: Unit test suite using Python's unittest framework
- `src/_bindings.cpp`: pybind11 interface between Python and C++
- `src/_wrapper.cpp/_wrapper.h`: C++ wrapper around TM-align functionality

## Development Commands

**Installation for development:**
```bash
# Create virtual environment (recommended)
python3 -m venv tmtools-dev
source tmtools-dev/bin/activate

# Install in development mode
pip install -e . -v

# Install optional dependencies
pip install biopython
```

**Run tests:**
```bash
python -m unittest discover -v .
```

**Run single test file:**
```bash
python -m unittest tmtools.tests.test_io -v
```

**Linting and formatting:**
```bash
ruff check .
```

**Build requirements:**
- C++ compiler with C++14 support
- pybind11 ~= 2.10.4
- NumPy

## Testing

- Tests use the standard unittest framework with given/when/then pattern
- Test data includes bundled PDB files accessible via `tmtools.testing.get_pdb_path()`
- BioPython is required for IO-related tests but not core functionality
- Tests are run from a separate directory to ensure proper package installation

## Key Dependencies

- **Required**: numpy
- **Optional**: biopython (for PDB/mmCIF file reading)
- **Build**: setuptools, wheel, pybind11, C++ compiler

## Code Style

- Uses ruff for linting and formatting (configuration in pyproject.toml)
- Follows E, F, W rule selections from ruff
- CI automatically checks formatting via GitHub Actions

## Release Process

**1. Update changelog:**
- Add release notes to `CHANGES.md` following the existing format
- Include version number as heading
- List changes with PR references (e.g., `[GH #XX]`)

**2. Update version numbers:**
- `setup.py` - Update the `version` parameter (currently line 23)
- `tmtools/__init__.py` - Update the `__version__` variable (currently line 3)

**3. Commit and tag:**
- Commit the version number changes
- Create an annotated git tag with format `vX.Y.Z` (e.g., `v0.3.0`)
- The tag should point to the commit that changes the version numbers
- Example: `git tag -a v0.3.0 -m "Release v0.3.0" && git push origin v0.3.0`

**4. Create GitHub release:**
- Publish a new GitHub release using the tag created above
- This automatically triggers the build and upload workflow

**5. Automated build:**
- CI automatically builds wheels for Python 3.10-3.14 on Linux, macOS (x86_64/arm64), and Windows
- Creates source distribution
- Uploads to PyPI using the `PYPI_PASSWORD` secret
- All configured in `.github/workflows/build-and-upload.yml`