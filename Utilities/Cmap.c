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

Title: Cmap.c. C version of maputils.py

Author: Nariman Habili
Email: nariman.habili@ga.gov.au
CreationDate: 2008-01-15

Description: Contains mapping functions that operate on arrays
Supercedes some of the functions lying around in some unusual places (like
DataProcess)
SeeAlso:
Constraints:

To compile:

    (Windows)
	gcc -c Cmap.c -I"C:\Python25" -o Cmap.o -Wall -O
	gcc -shared Cmap.o -L"C:\Python25\libs" -l Python25 -o Cmap.pyd
    
    (Unix)
    gcc -c Cmap.c -I/usr/local/include/python2.5 -I/usr/local/lib/python2.5 -o Cmap.o -Wall -O -fPIC
    gcc -shared Cmap.o -o Cmap.so

Version: $Rev: 512 $

$Id: Cmap.c 512 2011-10-31 07:20:38Z nsummons $
*/

#include "Python.h"
#include <math.h>

int _bear2LatLon(double bearing, double distance, double oLon, double oLat, double* nLon, double* nLat)
{
	double radius = 6367.0; //Earth radius (km)
    double toRads = 0.0174532925;
	double toDegs = 57.2957795130;

	oLon = oLon*toRads;
    oLat = oLat*toRads;
    bearing = bearing*toRads;

	double sin0 = sin(oLat);
	double sin1 = sin(distance/radius);

	double cos0 = cos(oLat);
	double cos1 = cos(distance/radius);

    *nLat = asin(sin0*cos1 + cos0*sin1*cos(bearing));

	double temp0 = sin(bearing)*sin1*cos0;
    double temp1 = cos1 - sin0*sin(*nLat);

    *nLon = oLon + atan2(temp0, temp1);

    *nLon = *nLon * toDegs;
	*nLat = *nLat * toDegs;

	return 0;
}

PyObject* bear2LatLon(PyObject* self, PyObject* args)
{
	double bearing;
	double distance;
	double oLon;
	double oLat;
	double nLon;
	double nLat;

	// Convert Python arguments to C
	if (!PyArg_ParseTuple(args, "dddd", &bearing, &distance, &oLon, &oLat))
	{
        PyErr_SetString(PyExc_RuntimeError, "bear2LatLon could not parse input");
		return NULL;
	}

	//Call underlying routine
	_bear2LatLon(bearing, distance, oLon, oLat, &nLon, &nLat);

	//NOTE: return number of points inside..
	return Py_BuildValue("(dd)", nLon, nLat);
}


// Method table for python module
static struct PyMethodDef MethodTable[] =
{
		/* The cast of the function is necessary since PyCFunction values
		* only take two PyObject* parameters, and rotate() takes
		* three.
		*/
		{"bear2LatLon", bear2LatLon, METH_VARARGS, "Print out"},
		{NULL, NULL, 0, NULL}   /* sentinel */
};

// Module initialisation
void initCmap(void)
{
		Py_InitModule("Cmap", MethodTable);
}
