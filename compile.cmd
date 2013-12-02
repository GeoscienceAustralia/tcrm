@echo off
echo.
echo Compiling C code for TCRM
echo -------------------------
@copy installer\setup.py .
@copy installer\setup.cfg .
python setup.py build_ext -i
@del setup.py
@del setup.cfg
echo.
pause
