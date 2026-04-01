# Linux Offline Python Setup

This guide is for the offline Linux PC when you want to:
- run the full-year clonality scripts
- build the app from source on Linux
- avoid errors like `No module named yaml`

The `yaml` import comes from `PyYAML`, which is included in `requirements.txt`.

## 1. Download Linux Wheels On A Machine With Internet

From the repo root:

Preferred on macOS:

```bash
cd /path/to/OUS
bash packaging/download_linux_wheels_docker.sh
```

Fallback without Docker:

```bash
cd /path/to/OUS
bash packaging/download_linux_wheels.sh
```

This creates:

```bash
packaging/linux_offline_deps/
```

That folder contains the Python wheels for:
- runtime dependencies from `requirements.txt`
- build dependencies from `packaging/build-requirements.txt`

Why the Docker version is better:
- it downloads the wheels from inside a real Linux environment
- that avoids macOS host-resolution issues for packages like `kiwisolver`
- it is the recommended path when the host-side downloader says:
  `No matching distribution found`

## 2. Copy These To The Offline Linux PC

Copy to the Linux machine:
- the whole `OUS/` repo folder
- the folder `packaging/linux_offline_deps/`

## 3. Create A Python Environment On Linux

Use Python `3.10` if possible.

```bash
cd /path/to/OUS
python3.10 -m venv .venv
source .venv/bin/activate
python --version
```

If `python3.10` is not available but `python3` is version `3.10`, use:

```bash
python3 -m venv .venv
source .venv/bin/activate
python --version
```

## 4. Install Python Dependencies Offline

```bash
cd /path/to/OUS
source .venv/bin/activate
pip install --no-index --find-links=packaging/linux_offline_deps -r requirements.txt
```

That step installs `PyYAML`, which fixes:

```bash
ModuleNotFoundError: No module named yaml
```

## 5. Check That `yaml` Works

```bash
source .venv/bin/activate
python -c "import yaml; print(yaml.__version__)"
```

If this prints a version number, the `yaml` error is fixed.

## 6. Run The Full-Year Scripts On Linux

Example:

```bash
source .venv/bin/activate
python scripts/run_clonality_2025_monthly.py \
  --input-root /path/to/data/Klonalitet/2025_data \
  --output-root /path/to/output/full_year_runs \
  --run-name full_2025_run \
  --max-workers 1 \
  --folder-workers 1
```

Combine the result into one workbook:

```bash
source .venv/bin/activate
python scripts/combine_clonality_2025_overview.py \
  --run-root /path/to/output/full_year_runs/full_2025_run \
  --year-label 2025
```

## 7. Build The Linux App Offline

If you want to build the desktop app from source on Linux too:

1. install Linux system packages
2. install offline Python wheels
3. build the app

System packages:

```bash
sudo bash packaging/linux_system_deps.sh
```

Then build:

```bash
source .venv/bin/activate
python build_qt.py
```

## 8. Minimal Fix If You Only Need The `yaml` Error Gone

If the environment already exists and only `yaml` is missing:

```bash
source .venv/bin/activate
pip install --no-index --find-links=packaging/linux_offline_deps PyYAML
```

But the better option is still:

```bash
pip install --no-index --find-links=packaging/linux_offline_deps -r requirements.txt
```

because that keeps the whole environment complete.
