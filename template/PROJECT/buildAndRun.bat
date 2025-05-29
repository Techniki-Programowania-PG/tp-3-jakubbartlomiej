@echo off
setlocal

REM Change to the directory where this script is located
cd /d "%~dp0"

echo [INFO] Installing the package locally from: %cd%
pip install --no-cache-dir --force-reinstall .

REM Change to the tests folder
cd tests

echo [INFO] Running test_basic.py with pytest...
python -m pytest test_basic.py -v -s

endlocal
