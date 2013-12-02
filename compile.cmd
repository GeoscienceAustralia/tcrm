@echo off
echo.
echo Compiling C code for TCRM
echo -------------------------
@copy installer\setup.py . >nul
@copy installer\setup.cfg . >nul
python setup.py build_ext -i
@del setup.py
@del setup.cfg
echo.
pause
