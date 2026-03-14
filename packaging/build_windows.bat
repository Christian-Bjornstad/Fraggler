@echo off
REM ============================================================
REM Build Fraggler Diagnostics — Windows Desktop Bundle
REM ============================================================
REM Prerequisites:
REM   - Python 3.10+ installed
REM   - Run from the OUS\ project root
REM
REM Usage:
REM   packaging\build_windows.bat
REM
REM Output:
REM   dist\Fraggler_Windows and dist\releases\Fraggler_Windows.zip
REM ============================================================

echo ============================================================
echo   Building Fraggler Diagnostics for Windows
echo ============================================================

cd /d "%~dp0\.."

REM Create and activate venv if needed
if not exist "fraggler-win-venv" (
    echo Creating virtual environment...
    python -m venv fraggler-win-venv
)

call fraggler-win-venv\Scripts\activate.bat

REM Install dependencies
pip install -r requirements.txt
pip install pyinstaller

REM Clean previous builds
if exist build rmdir /s /q build
if exist dist rmdir /s /q dist

REM Build
echo.
echo Running unified desktop build...
echo.

python build_qt.py

echo.
echo ============================================================
echo   Build complete!
echo   Folder: dist\Fraggler_Windows
echo   Zip   : dist\releases\Fraggler_Windows.zip
echo.
echo   To run:
echo     dist\Fraggler_Windows\Fraggler.exe
echo ============================================================

pause
