set PYTHON=C:\Python27
set DEST=dist

@rem copy the important files to the root, so we don't have to hardcode paths
@rem all over the place

@if exist "installer\build.bat" goto continue
@echo You must run build.bat from within the root directory
:continue

@cd installer

@copy setup.py ..
@copy setup.cfg ..

@cd ..

@del /F /S /Q build
@del /F /S /Q %DEST%

%PYTHON%\python.exe setup.py py2exe
@if errorlevel 1 goto error

@rem %PYTHON%\python.exe installer\winprepnsi.py windows_installer\installer.nsi installer.temp.nsi
@rem if errorlevel 1 goto error
@rem "C:\Program Files\NSIS\makensis.exe" installer.temp.nsi
@rem if errorlevel 1 goto error
@rem del installer.temp.nsi
@rem if errorlevel 1 goto error

@rem cleanup
@del setup.py
@del setup.cfg
@del /F /S /Q build

@rem comment out the next line if you do not want example data included.
@rem xcopy input\*.* %DEST%\input\

@goto done

:error
@echo ------------------------------------------------------------------------
@echo Build failed.

:done
