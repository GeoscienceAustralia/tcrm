@echo off

rem Set title of command prompt
title Tropical Cyclone Risk Model
goto RunSimulation

:RunSimulation
cls
echo.
echo    Tropical Cyclone Risk Model (TCRM) - Version 1.0 (beta release)
echo    Copyright (C) 2011 Commonwealth of Australia (Geoscience Australia)
echo.
echo    This program is free software: you can redistribute it and/or modify
echo    it under the terms of the GNU General Public License as published by
echo    the Free Software Foundation, either version 3 of the License, or
echo    (at your option) any later version.
echo.
echo    This program is distributed in the hope that it will be useful,
echo    but WITHOUT ANY WARRANTY; without even the implied warranty of
echo    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
echo    GNU General Public License for more details.
echo.
echo    You should have received a copy of the GNU General Public License
echo    along with this program.  If not, see http://www.gnu.org/licenses/.  
echo.
echo.
set %configfile=''
echo Please enter a configuration file (or press 'Q' to exit) and
set /P configfile= then press enter: %=%
if '%configfile%' == 'q' goto EndProgram
if '%configfile%' == 'Q' goto EndProgram
if '%configfile%' == 'quit' goto EndProgram
if '%configfile%' == 'exit' goto EndProgram
if %configfile% == '' goto RunSimulation
if not exist %configfile% goto FILENOTFOUND
goto StartMainpy

:StartMainpy
echo.
echo.
echo ---------  Starting Hazard Simulation  ---------
python -Wignore tcrm.py -c %configfile%
echo ------------------------------------------------
echo.
echo TCRM has finished running.
echo Please check the log file to ensure that the simulation completed successfully.
echo.
pause
goto EndProgram

:FILENOTFOUND
echo.
echo Configuration file '%configfile%' not found
pause
goto RunSimulation

:EndProgram
echo.
goto EOF

:EOF
