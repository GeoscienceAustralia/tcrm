/*
    Tropical Cyclone Risk Model (TCRM) - Version 1.0 (beta release)
    Copyright (C) 2011  Geoscience Australia

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

Title: Cstats.c - helper functions for statistical methods. This is the C version of the
       Python file stat.py
Author: Nariman Habili, nariman.habili@ga.gov.au
CreationDate: 2008-01-15
Description: Miscellaneous tools required for statistics-related classes.
SeeAlso:
Constraints:
To compile:

    (Windows)
    gcc -c Cstats.c -I"C:\Python25" -o Cstats.o -Wall -O -std=c99
    gcc -shared Cstats.o -L"C:\Python25\libs" -l Python25 -o Cstats.pyd

    (Unix)
    gcc -c Cstats.c -I/usr/local/include/python2.5 -I/usr/local/lib/python2.5 -o Cstats.o -Wall -O -fPIC -std=c99
    gcc -shared Cstats.o -o Cstats.so

Version: 68
ModifiedBy: Craig Arthur, craig.arthur@ga.gov.au
ModifiedDate: 2009-01-21
Modification: Added getCellLonLat()

Version: $Rev: 512 $

$Id: Cstats.c 512 2011-10-31 07:20:38Z nsummons $
*/

#include "Python.h"
#include <math.h>

int _getCellNum(double oLon, double oLat, int xMin, int xMax, int yMin, int yMax, int x, int y);
void _getCellLonLat(int cellNum, int xMin, int xMax, int yMin, int yMax, int x, int y, int nx, int ny, double *oLon, double *oLat);

int _getCellNum(double oLon, double oLat, int xMin, int xMax, int yMin, int yMax, int x, int y)
{
    int lon;
    int lat;
    double j;
    double i;

    lon = (int)(floor(oLon));
    lat = (int)(ceil(oLat));

    if (lon < xMin || lon >= xMax || lat <= yMin || lat > yMax)
	{
		return -1;
	}

	j = abs(abs((double)(lon)) - abs((double)(xMin)))/abs((double)(x));
	i = abs(abs((double)(lat)) - abs((double)(yMax)))/abs((double)(y));

	return (int)(i*abs((xMax - xMin)/x) + j);
}

void _getCellLonLat(int cellNum, int xMin, int xMax, int yMin, int yMax, int x, int y, int nx, int ny, double *oLon, double *oLat)
{
    double lon[nx];
    double lat[ny];
    int ix;
    int jy;

    for (int i=0; i<nx; i++)
    {
        lon[i] = xMin + (double)i*x;
    }

    for (int j=0; j<ny; j++)
    {
        lat[j] = yMax - (double)j*y;
    }
    jy = cellNum/nx;
    ix = cellNum%nx;
    *oLon = lon[ix];
    *oLat = lat[jy];
    return;
}

PyObject* getCellNum(PyObject* self, PyObject* args)
{
	double oLon;
	double oLat;
	PyObject* gridLimit;
	PyObject* gridSpace;
	int xMin;
	int xMax;
	int yMin;
	int yMax;
	int x;
	int y;
	int cell;

	// Convert Python arguments to C
	if (!PyArg_ParseTuple(args, "ddOO", &oLon, &oLat, &gridLimit, &gridSpace))
	{
        PyErr_SetString(PyExc_RuntimeError, "Cstat_tool could not parse input");
		return NULL;
	}

	xMin = PyInt_AsLong(PyDict_GetItemString(gridLimit, "xMin"));
	xMax = PyInt_AsLong(PyDict_GetItemString(gridLimit, "xMax"));
	yMin = PyInt_AsLong(PyDict_GetItemString(gridLimit, "yMin"));
	yMax = PyInt_AsLong(PyDict_GetItemString(gridLimit, "yMax"));

	x = PyInt_AsLong(PyDict_GetItemString(gridSpace, "x"));
	y = PyInt_AsLong(PyDict_GetItemString(gridSpace, "y"));


	//Call underlying routine
	cell = _getCellNum(oLon, oLat, xMin, xMax, yMin, yMax, x, y);
	if (cell < 0)
	{
        PyErr_SetString(PyExc_ValueError, "Invalid input on cellNum: cell number is out of range");
		return NULL;
	}

	return Py_BuildValue("i", cell);
}

PyObject* getCellLonLat(PyObject* self, PyObject* args)
{
        int cellNum;
        PyObject* gridLimit;
        PyObject* gridSpace;
        int xMin;
        int xMax;
        int yMin;
        int yMax;
        int x;
        int y;
        int nx;
        int ny;
        double oLon;
        double oLat;

        if (!PyArg_ParseTuple(args, "iOO", &cellNum, &gridLimit, &gridSpace))
        {
        PyErr_SetString(PyExc_RuntimeError, "getCellLonLat could not parse input");
            return NULL;
        }
	xMin = PyInt_AsLong(PyDict_GetItemString(gridLimit, "xMin"));
	xMax = PyInt_AsLong(PyDict_GetItemString(gridLimit, "xMax"));
	yMin = PyInt_AsLong(PyDict_GetItemString(gridLimit, "yMin"));
	yMax = PyInt_AsLong(PyDict_GetItemString(gridLimit, "yMax"));

	x = PyInt_AsLong(PyDict_GetItemString(gridSpace, "x"));
	y = PyInt_AsLong(PyDict_GetItemString(gridSpace, "y"));
        if (cellNum < 0)
	{
        PyErr_SetString(PyExc_ValueError, "Invalid input on cellNum: cell number is out of range");
		return NULL;
	}
        nx = (xMax - xMin)/x;
        ny = (yMax - yMin)/y;
        _getCellLonLat(cellNum, xMin, xMax, yMin, yMax, x, y, nx, ny, &oLon, &oLat);

        return Py_BuildValue("dd", oLon, oLat );
}
// Method table for python module
static struct PyMethodDef MethodTable[] =
{
	{"getCellNum", getCellNum, METH_VARARGS, "Print out"},
        {"getCellLonLat", getCellLonLat, METH_VARARGS, "Print out"},
	{NULL, NULL, 0, NULL}
};

// Module initialisation
void initCstats(void)
{
        Py_InitModule("Cstats", MethodTable);
}
