@echo off
cd /d "%~dp0"
if not exist .venv\Scripts\python.exe (
  echo Execute instalar.bat primeiro.
  pause
  exit /b 1
)
.venv\Scripts\python.exe -c "import keyring" >nul 2>&1
if errorlevel 1 .venv\Scripts\python.exe -m pip install -r requirements.txt
.venv\Scripts\python.exe app.py
if errorlevel 1 pause
