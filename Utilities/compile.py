"""
    Tropical Cyclone Risk Model (TCRM) - Version 1.0 (beta release)
    Copyright (C) 2011 Commonwealth of Australia (Geoscience Australia)

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""
import os, sys
import numpy
import platform

operatingSystem = platform.system()
pythonPath = sys.prefix
numpyPath = os.path.join(numpy.get_include(), 'numpy')
pythonVersion = sys.version_info

# Move to directory of current file
os.chdir(os.path.dirname(os.path.realpath(__file__)))

if operatingSystem == 'Windows':
 
    pythonIncludePath = '"' + os.path.join(pythonPath, 'include') + '"'
    pythonLibPath = '"' + os.path.join(pythonPath, 'libs') + '"'
    pythonLibFile = 'python' + str(pythonVersion[0]) + str(pythonVersion[1])
   
    #build Cstats
    cmd1 = 'gcc -c Cstats.c -I' + pythonIncludePath + ' -o Cstats.o -Wall -O -std=c99'
    cmd2 = 'gcc -shared Cstats.o -L' + pythonLibPath + ' -l ' + pythonLibFile + ' -o Cstats.pyd'
    print cmd1
    os.system(cmd1)
    print cmd2
    os.system(cmd2)
    
    #build Cmap    
    cmd1 = 'gcc -c Cmap.c -I' + pythonIncludePath + ' -o Cmap.o -Wall -O'
    cmd2 = 'gcc -shared Cmap.o -L' + pythonLibPath + ' -l ' + pythonLibFile + ' -o Cmap.pyd'
    print cmd1
    os.system(cmd1)
    print cmd2
    os.system(cmd2)

    #build KPDF
    cmd1 = 'gcc -c KPDF.c -I' + pythonIncludePath +  ' -I"' + numpyPath + '" -o KPDF.o -Wall -O -std=c99'
    cmd2 = 'gcc -shared KPDF.o -L' + pythonLibPath + ' -l ' + pythonLibFile + ' -o KPDF.pyd'
    print cmd1
    os.system(cmd1)
    print cmd2
    os.system(cmd2)

elif operatingSystem == 'Linux':
    
    pythonIncludePath = os.path.join(pythonPath, 'include', 'python' + str(pythonVersion[0]) + '.' + str(pythonVersion[1]))
    pythonLibPath = os.path.join(pythonPath, 'lib', 'python' + str(pythonVersion[0]) + '.' + str(pythonVersion[1]))
    
    #build Cstats
    cmd1 = 'gcc -c Cstats.c -I' + pythonIncludePath + ' -I' + pythonLibPath + ' -o Cstats.o -Wall -O -fPIC -std=c99'
    cmd2 = 'gcc -shared Cstats.o -o Cstats.so'
    print cmd1
    os.system(cmd1)
    print cmd2
    os.system(cmd2)

    #build Cmap
    cmd1 = 'gcc -c Cmap.c -I' + pythonIncludePath + ' -I' + pythonLibPath + ' -o Cmap.o -Wall -O -fPIC'
    cmd2 = 'gcc -shared Cmap.o -o Cmap.so'
    print cmd1
    os.system(cmd1)
    print cmd2
    os.system(cmd2)
    
    #build KPDF
    cmd1 = 'gcc -c KPDF.c -I' + pythonIncludePath + ' -I' + pythonLibPath + ' -I' + numpyPath + ' -o KPDF.o -Wall -O -fPIC -std=c99'
    cmd2 = 'gcc -shared KPDF.o -o KPDF.so'
    print cmd1
    os.system(cmd1)
    print cmd2
    os.system(cmd2)
