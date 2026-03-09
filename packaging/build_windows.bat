@echo off
REM ============================================================
REM Build Fraggler Diagnostics — Windows Executable
REM ============================================================
REM Prerequisites:
REM   - Python 3.10+ installed
REM   - Run from the OUS\ project root
REM
REM Usage:
REM   packaging\build_windows.bat
REM
REM Output:
REM   packaging\dist\fraggler-diagnostics\
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
if exist packaging\build rmdir /s /q packaging\build
if exist packaging\dist rmdir /s /q packaging\dist

REM Build
echo.
echo Running PyInstaller...
echo.

pyinstaller packaging\fraggler_diagnostics.spec ^
    --distpath packaging\dist ^
    --workpath packaging\build ^
    --clean ^
    --noconfirm

echo.
echo ============================================================
echo   Build complete!
echo   Executable: packaging\dist\fraggler-diagnostics\
echo.
echo   To run:
echo     packaging\dist\fraggler-diagnostics\fraggler-diagnostics.exe
echo ============================================================

pause
