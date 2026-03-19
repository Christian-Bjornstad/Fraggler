# Fraggler Diagnostics — Linux Source Build Guide

Use this guide if the packaged Linux app does not run on the target machine and you want to rebuild it yourself.

## 1. Unpack the source bundle

```bash
unzip Fraggler_Source_with_Guide.zip -d ~/fraggler-source
cd ~/fraggler-source/OUS
```

## 2. Quickest rebuild path: Docker

This is the easiest way to rebuild the offline Linux app on a Linux machine with Docker installed.

```bash
./packaging/build_linux.sh
```

Expected output:

```bash
dist/Fraggler_Linux
dist/releases/Fraggler_Linux_offline.zip
```

## 3. Run directly from source

If you just want to launch Fraggler without packaging it first:

```bash
python3.10 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r packaging/build-requirements.txt
python -m pip install -r requirements.txt
python qt_app.py
```

## 4. Build without Docker

If Docker is not available, you can still build the Linux app directly on a Linux machine:

```bash
python3.10 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r packaging/build-requirements.txt
python -m pip install -r requirements.txt
python build_qt.py
```

Expected output:

```bash
dist/Fraggler_Linux
dist/releases/Fraggler_Linux_offline.zip
```

## 5. Notes

- Target platform: Linux `x86_64`
- Recommended runtime target: Fedora 35 or another Linux machine with `glibc >= 2.31`
- The packaged Linux app uses `QT_QPA_PLATFORM=xcb`
- Keep the full extracted bundle together, including the `_internal` folder
