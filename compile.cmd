@echo off
echo.
echo Compiling C code for TCRM
echo -------------------------
python installer\setup.py build_ext -i
echo.
pause
