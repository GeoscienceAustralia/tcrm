/****************
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

		 KPDF.c

    Kernel-based probability density function estimation.

    Mathematical details: B. W. Silverman, 1986, "Density Estimation for
    Statistics and Data Analysis", 1st ed., Chapman and Hall, London, 175 pp.

    Copyright (C) Jon Saenz, Jesus Fernandez, Juan Zubillaga, September, 2000.

    C-recoding of a previous Python extension using pure Python.
    It was extremely slow for two-dimensional estimations.


    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

    Jon Saenz
    Departamento de Fisica Aplicada II
    Facultad de Ciencias
    Universidad del Pais Vasco
    Apdo. 644
    48080 - Bilbao
    Spain
    jsaenz@wm.lc.ehu.es

    2009-08-03
    Craig Arthur, craig.arthur@ga.gov.au
    Added Gaussian kernel to UPDF
	Version: 286

    2010-09-16
	Craig Arthur, craig.arthur@ga.gov.au
	Updated to use PyArray_SimpleNew, replacing PyArray_FromDims
	Version: $Rev: 512 $

 To compile:

    (Windows)
    gcc -c KPDF.c -I"C:\Python25" -o KPDF.o -Wall -O -fPIC -std=c99
    gcc -shared KPDF.o -L"C:\Python25\libs" -l Python25 -o KPDF.pyd

    (Unix)
    gcc -c KPDF.c -I/usr/local/include/python2.5 -I/usr/local/lib/python2.5 -o KPDF.o -Wall -O -fPIC -std=c99
    gcc -shared KPDF.o -o KPDF.so

    $Id: KPDF.c 512 2011-10-31 07:20:38Z nsummons $

***************/

#include "Python.h"
/* #include "numpy/arrayobject.h" */
#include "arrayobject.h" 

#include <math.h>
#include <stdlib.h>


/* Some declarations which will be needed soon */

/******
These constants, of internal use, allow the selection
of different kernels
 ****/
typedef double (*KernelFunctionPtr)( double x );
typedef enum { updfepanechnikov, updfbiweight, updftriangular, updfgaussian } UPDFMode;
typedef double (*MKernelFunctionPtr)( double x );
typedef enum { mpdfepanechnikov , mpdfmvgaussian } MPDFMode;

/* Initialize some internal constants when the module is first imported */
static void initKPDFconstants( void );

/* Forward declaration of the prototypes needed in the method
initialization table */
static PyObject *UPDFEpanechnikov( PyObject *self , PyObject *args );
static PyObject *UPDFBiweight( PyObject *self , PyObject *args );
static PyObject *UPDFTriangular( PyObject *self , PyObject *args );
static PyObject *UPDFGaussian( PyObject *self , PyObject *args );
static PyObject *UPDFOptimumBandwidth( PyObject *self , PyObject *args );

static PyObject *MPDFEpanechnikov( PyObject *self , PyObject *args );
static PyObject *MPDFGaussian( PyObject *self , PyObject *args );
static PyObject *MPDFOptimumBandwidth( PyObject *self, PyObject *args );

static PyObject *MPDF2DGrid2Array( PyObject *self, PyObject *args );
static PyObject *MPDF3DGrid2Array( PyObject *self, PyObject *args );
/*static PyObject *set_callback( PyObject *self, PyObject *args );*/

/* Declare this function because it is accessed before it is defined */
static PyObject *MPDFGrid2Array( PyArrayObject **numpys, int fastest );

/* Define function for passing python callback function
used for updating progress bar. */
static PyObject *my_callback = NULL;
static PyObject *
set_callback(dummy, arg)
    PyObject *dummy, *arg;
{
    Py_XDECREF(my_callback); /* Dispose of previous callback */
    Py_XINCREF(arg); /* Add a reference to new callback */
    my_callback = arg; /* Remember new callback */
    /* Boilerplate to return "None" */
    Py_INCREF(Py_None);
    return Py_None;
}

/***
Initialize the module and allow access from Python to some functions
 **/
static PyMethodDef KPDFMethods[]={
  {"UPDFEpanechnikov",UPDFEpanechnikov,METH_VARARGS},
  {"UPDFBiweight",UPDFBiweight,METH_VARARGS},
  {"UPDFTriangular",UPDFTriangular,METH_VARARGS},
  {"UPDFGaussian",UPDFGaussian,METH_VARARGS},
  {"UPDFOptimumBandwidth",UPDFOptimumBandwidth,METH_VARARGS},
  {"MPDFEpanechnikov",MPDFEpanechnikov,METH_VARARGS},
  {"MPDFGaussian",MPDFGaussian,METH_VARARGS},
  {"MPDFOptimumBandwidth",MPDFOptimumBandwidth,METH_VARARGS},
  {"MPDF2DGrid2Array",MPDF2DGrid2Array,METH_VARARGS},
  {"MPDF3DGrid2Array",MPDF3DGrid2Array,METH_VARARGS},
  {"set_callback",set_callback},
  {NULL,NULL} /* Last entry */
};

void initKPDF( void )
{
  PyObject *m, *d;

  /* Initialize the Python Module */
  m=Py_InitModule("KPDF",KPDFMethods);
  /* Give access to Numeric Arrays */
  import_array();
  /* Intialize the dictionary */
  d=PyModule_GetDict(m);
  initKPDFconstants();
}

/*****
Initialize some internal constants for all the kernels.
Just once at the start.
******/
static double __sqrt5=0.0;
static double __minussqrt5=0.0;
static double __sqrt2pi=0.0;
#define KPDFMAXDIMS 3
static double __Cd[KPDFMAXDIMS];

static void initKPDFconstants( void )
{
  __sqrt5=sqrt(5.);
  __minussqrt5=-1*__sqrt5;
  __sqrt2pi=sqrt(2.*3.141592658);
  __Cd[0]=2.;
  __Cd[1]=acos(-1.);
  __Cd[2]=4*acos(-1.)/3.;
}

static PyObject *UPDFOptimumBandwidth( PyObject *self , PyObject *args )
{
  PyArrayObject *data;
  PyArrayObject *cdata;
  int N,i;
  double *xptr, mean, variance, std, hopt;

  if (!PyArg_ParseTuple(args,"O",&data))
    return NULL;
  cdata=(PyArrayObject*)PyArray_ContiguousFromObject((PyObject*)data,PyArray_DOUBLE,1,1);
  if (cdata==NULL)
    return NULL;
  if (cdata->nd!=1){
    PyErr_SetString(PyExc_TypeError,"Univariate Optimum bandwidth estimation:1-D arrays required");
    return NULL;
  }
  N=cdata->dimensions[0];
  mean=0.;
  for(i=0;i<N;i++){
    xptr=(double*)(cdata->data+i*cdata->strides[0]);
    mean+=*xptr;
  }
  mean=mean/N;
  variance=0.0;
  for(i=0;i<N;i++){
    xptr=(double*)(cdata->data+i*cdata->strides[0]);
    variance+=(*xptr-mean)*(*xptr-mean);
  }
  variance=variance/((double)N);
  std=sqrt(variance);
  hopt=pow(4./3.,0.2)*std/pow((double)N,0.2);
  Py_DECREF(cdata);
  return Py_BuildValue("d",hopt);
}

/*****
General calling order to these functions:
1st argument: Dataset (1-D NumpyArray).
2nd argument: Points to return the PDF (1-D NumpyArray)
3rd argument: Bandwidth
 ***/
static int parseUPDFarguments( PyObject *self, PyObject *args,
			   PyArrayObject **Xi, PyArrayObject **x,
			   double *bandwidth )
{
  if (!PyArg_ParseTuple(args,"OOd",Xi,x,bandwidth))
    return 0;
  else
    return 1;
}

/************
The next functions are internal, too.
Just different kernels to be called from the C functions
I want to avoid iterations in Python, so, I simply avoid allowing
access to these functions and iterate over the Numpy arrays using C
cycles. The alternative using Python is quite time consuming,
mainly for two-dimensional (and up-dimensional) spaces.
************/
static double updfkekernel( double x )
{
  double dummy=0.0;
  double temp;
  if ((x >= __minussqrt5) && (x <= __sqrt5)){
    temp=x/__sqrt5;
    /** Avoid catastrophic cancellation **/
    dummy=(1.-temp)*(1.+temp);
  }
  return dummy;
}

static double updfbwkernel( double x )
{
  double dummy=0.0;
  if ((x>=-1.)&&(x<=1.))
    dummy=pow((1.-x)*(1.+x),2.);
  return dummy;
}

static double updftkernel( double x )
{
  double dummy=0.0;
  if ((x>=-1.)&&(x<=1.))
    dummy=1.-fabs(x);
  return dummy;
}

static double updfgkernel( double x )
{
  double dummy=0.0;
  dummy=exp(-0.5*fabs(x));
  return dummy;
}

static KernelFunctionPtr updfselectkf( UPDFMode mode, double *k )
{
  KernelFunctionPtr f;

  switch (mode){
  case updfepanechnikov:
    f=updfkekernel;
    *k=3./4./__sqrt5;
    break;
  case updfbiweight:
    f=updfbwkernel;
    *k=15./16.;
    break;
  case updftriangular:
    f=updftkernel;
    *k=1.;
    break;
  case updfgaussian:
    f=updfgkernel;
    *k=1./__sqrt2pi;
    break;
  default:
    /*** Return null here ****/
    PyErr_SetString(PyExc_ValueError,"Only Epanechnikov, Biweight, Gaussian and Triangular kernels allowed");
    f=NULL;
    *k=1.;
    break;
  }
  return f;
}

/********
This function iterates over the input Python collection 'data'
(Numeric array) and, for each value in the input collection
(Numeric array) 'xpoints', evaluates the PDF kernel based estimation
and returns the PDF at each point using the
requested bandwidth. It returns, therefore, a Numeric array of the same
shape as 'xpoints'
*******/
static PyObject *UPDFEstimator( PyArrayObject *data , PyArrayObject *xpoints ,
				double bandwidth , UPDFMode mode )
{
  KernelFunctionPtr f;
  double theconstant;
  PyArrayObject *rarray;
  PyArrayObject *cdata;
  PyArrayObject *cxpoints;
  int i,N,Nxpoints,j;
  npy_intp dimensions[1]={0};
  double sum, *xpos, *dpos, *rpos;
  double x, Xj, u;

  /** Test whether the input arrays are linear ***/
  if (data->nd!=1 || xpoints->nd!=1){
    PyErr_SetString(PyExc_TypeError,"Univariate PDF estimator: 1-D Numeric arrays required");
    return NULL;
  }
  /** Select the kernel and normalizing constant ***/
  f=updfselectkf(mode,&theconstant);
  if (f==NULL)
    return NULL;
  /** Allocate internal contiguous arrays of type DOUBLE */
  cdata=(PyArrayObject*)PyArray_ContiguousFromObject(
						     (PyObject*)data,
						     PyArray_DOUBLE,
						     1,1);
  if(cdata==NULL)
    return NULL;
  cxpoints=(PyArrayObject*)PyArray_ContiguousFromObject(
						     (PyObject*)xpoints,
						     PyArray_DOUBLE,
						     1,1);
  if(cxpoints==NULL){
    Py_DECREF(cdata);
    return NULL;
  }
  /*** These numbers are the items in each of the input arrays **/
  N=cdata->dimensions[0];
  Nxpoints=cxpoints->dimensions[0];
  dimensions[0]=Nxpoints;
  /* Create the returned array */
  rarray=(PyArrayObject*)PyArray_SimpleNew(1,dimensions,PyArray_DOUBLE);
  if (rarray==NULL)
    return NULL;
  /**
      For each point in the target array, iterate through the
      input data array
  **/
  for(i=0;i<Nxpoints;i++){
    sum=0.0;
    xpos=(double*)(cxpoints->data+i*cxpoints->strides[0]);
    rpos=(double*)(rarray->data+i*rarray->strides[0]);
    x=*xpos;
    for(j=0;j<N;j++){
      dpos=(double*)(cdata->data+j*cdata->strides[0]);
      Xj=*dpos;
      u=(x-Xj)/bandwidth;
      sum+=(*f)(u);
    }
    *rpos=sum*theconstant/N/bandwidth;
  }
  Py_DECREF(cdata);
  Py_DECREF(cxpoints);
  return (PyObject*)PyArray_Return((PyArrayObject*)rarray);
}


static PyObject *UPDFEpanechnikov( PyObject *self , PyObject *args )
{
  PyArrayObject *Xi;
  PyArrayObject *x;
  double bandwidth;
  if (!parseUPDFarguments(self,args,&Xi,&x,&bandwidth))
    return NULL;
  else
    return UPDFEstimator(Xi,x,bandwidth,updfepanechnikov);
}

static PyObject *UPDFBiweight( PyObject *self , PyObject *args )
{
  PyArrayObject *Xi;
  PyArrayObject *x;
  double bandwidth;
  if (!parseUPDFarguments(self,args,&Xi,&x,&bandwidth))
    return NULL;
  else
    return UPDFEstimator(Xi,x,bandwidth,updfbiweight);
}

static PyObject *UPDFTriangular( PyObject *self , PyObject *args )
{
  PyArrayObject *Xi;
  PyArrayObject *x;
  double bandwidth;
  if (!parseUPDFarguments(self,args,&Xi,&x,&bandwidth))
    return NULL;
  else
    return UPDFEstimator(Xi,x,bandwidth,updftriangular);
}

static PyObject *UPDFGaussian( PyObject *self , PyObject *args )
{
  PyArrayObject *Xi;
  PyArrayObject *x;
  double bandwidth;
  if (!parseUPDFarguments(self,args,&Xi,&x,&bandwidth))
    return NULL;
  else
    return UPDFEstimator(Xi,x,bandwidth,updfgaussian);
}
/*****
General calling order to these functions:
1st argument: Dataset (n-D NumpyArray).
2nd argument: Points to return the PDF (n-D NumpyArray)
3rd argument: Bandwidth (double)
4th (optional) argument: Inverse of the covariance matrix
5th (optional) argument: Squareroot of the determinant of S
 ***/
static int parseMPDFarguments( PyObject *self, PyObject *args,
			       PyArrayObject **Xi, PyArrayObject **x,
			       double *bandwidth,
			       PyArrayObject **Sm1,
			       double *sqrtdetS )
{
  /** Assign NULL to default arguments */
  *sqrtdetS=0.0;
  *Sm1=NULL;

  if (!PyArg_ParseTuple(args,"OOd|Od",Xi,x,bandwidth,Sm1,sqrtdetS))
    return 0;
  else
    return 1;
}


static double mekernel( double x )
{
  if (x<1.)
    return (1-x);
  else
    return 0;
}

static double mgkernel(double x )
{
  return exp(-0.5*x);
}

/* Select the kernel function and the dimension-independent
   normalizing constant associated to each kernel function */
static MKernelFunctionPtr mpdfselectkf( MPDFMode mode, double *k )
{
  MKernelFunctionPtr f;

  switch (mode){
  case mpdfepanechnikov:
    f=mekernel;
    *k=0.5;
    break;
  case mpdfmvgaussian:
    f=mgkernel;
    *k=1.;
    break;
  default:
    PyErr_SetString(PyExc_ValueError,"Only Epanechnikov or Gaussian kernels allowed");
    f=NULL;
    *k=0.0;
    break;
  }
  return f;
}

static int adequatedimensions(PyArrayObject *data, PyArrayObject *xpoints,
			      PyArrayObject *Sm1 )
{
  int OK=1;
  /* 1st. Are the input datasets 2-dimensional arrays? */
  if (data->nd>2 || xpoints->nd>2 ){
    PyErr_SetString(PyExc_ValueError,
		    "For Multivariate PDF Estimation, data arrays must be 2-D");
    OK=0;
  }
  /* 2nd. Is the second axis shorter than the maximum allowed? (KPDFMAXDIMS)*/
  if (data->dimensions[1]>KPDFMAXDIMS || xpoints->dimensions[1]>KPDFMAXDIMS){
    PyErr_SetString(PyExc_ValueError,
		    "Multivariate PDF: Dimensionality too high.");
    OK=0;
  }
  /* Test whether both datasets have the same length in last dimension */
  if (data->dimensions[1]!=xpoints->dimensions[1]){
    PyErr_SetString(PyExc_ValueError,
		    "Multivariate PDF: Dimension of last axis must \
be equal for both data sets.");
    OK=0;
  }
  if (Sm1!=NULL){
    /* Test if the covariance matrix and the data sets are compatible */
    if (Sm1->nd!=2){
      PyErr_SetString(PyExc_ValueError,"MPDF Fukunaga: Covariance matrix must be two-dimensional!!");
      OK=0;
    }
    if (Sm1->dimensions[0]!=Sm1->dimensions[1]){
      PyErr_SetString(PyExc_ValueError,"MPDF Fukunaga: Covariance matrix must be square");
      OK=0;
    }
    if (Sm1->dimensions[0]!=data->dimensions[1]){
      PyErr_SetString(PyExc_ValueError,"MPDF Fukunaga: Covariance matrix must be compatible for multiplication into data sets");
      OK=0;
    }
  }
  return OK;
}

/** Some people like BLAS and LAPACK, others only use LAPACK_LITE
 To avoid a mess during compilation, I write these very small functions
 They will not be so perfect, but allow a simpler Makefile and
 provide all the functionality actually needed in this module
*/
static void vdifference( double *a, double *b , double *c , int D ,
			 double bandwidth)
{
  register int i;
  for(i=0;i<D;i++)
    c[i]=(a[i]-b[i])/bandwidth;
}

static double dotprod( double *v , double *w , int D )
{
  register int i;
  double sum=0.0;
  for(i=0;i<D;i++)
    sum+=v[i]*w[i];
  return sum;
}

static void matrixtimesvector(PyArrayObject *Sm1,
			      double *v, double *w,
			      int D)
{
  register int i,j;
  register double sum, *dpos;

  for(i=0;i<D;i++){
    sum=0.0;
    for(j=0;j<D;j++){
      dpos=(double*)(Sm1->data+i*Sm1->strides[0]+j*Sm1->strides[1]);
      sum+=*dpos*v[j];
    }
    w[i]=sum;
  }
}

static double getDconstant( int D , MPDFMode mode )
{
  double K=0.0;
  if (D==1){
    PyErr_SetString(PyExc_ValueError,"For univariate data, you should use UPDF<kernel> estimators");
    return 0;
  }
  switch(mode){
  case mpdfepanechnikov:
    K=(D+2.)/__Cd[D-1];
    break;
  case mpdfmvgaussian:
    K=pow(2.*acos(-1.),-(((double)D)/2.));
    break;
  default:
    PyErr_SetString(PyExc_ValueError,"MPDF: Only Epanechnikov and Gaussian kernels allowed");
    break;
  }
  return K;
}

static PyObject *MPDFEstimator( PyArrayObject *data,
				PyArrayObject *xpoints,
				double bandwidth,
				PyArrayObject *Sm1,
				double sqrtdetS, MPDFMode mode)
{
  MKernelFunctionPtr f;
  double theconstant, Dconstant, K;

  PyArrayObject *cdata;
  PyArrayObject *cxpoints;
  PyArrayObject *rarray=NULL;
  int i, Nd, Nx,j, D;
  npy_intp dims[1]={0};
  double sum, *dpos, *xpos, *rpos;
  /* double *x, *X, u; */
  double vdif[KPDFMAXDIMS];
  double wdif[KPDFMAXDIMS];
  double normdiff;
  int isfukunaga;
  PyObject *arglist;
  PyObject *result;
  
  isfukunaga=(Sm1==NULL)?0:1;
  f=mpdfselectkf(mode,&theconstant);
  /* If no kernel is selected, pass the exception set by mpdfselectkf() */
  if (f==NULL)
    return NULL;
  /* Now, test whether the dimensionality of the problem exceeds
     current capabilities, and pass the exception */
  if (!adequatedimensions(data,xpoints,Sm1))
    return NULL;
  /* OK, everything should be correct now */
  Nd=data->dimensions[0];
  Nx=xpoints->dimensions[0];
  D=data->dimensions[1];
  Dconstant=getDconstant(D,mode);
  if (isfukunaga)
    K=Dconstant*theconstant/sqrtdetS/((double)Nd)/pow(bandwidth,((double)D));
  else
    K=Dconstant*theconstant/((double)Nd)/pow(bandwidth,((double)D));

  dims[0]=Nx;
  rarray=(PyArrayObject*)PyArray_SimpleNew(1,dims,PyArray_DOUBLE);
  if (rarray==NULL)
     return NULL;
  /* We need contiguous arrays for algebra operations */
  cdata=(PyArrayObject*)PyArray_ContiguousFromObject((PyObject*)data,
						     PyArray_DOUBLE,2,2);
  if (cdata==NULL){
    Py_DECREF(rarray);
    return NULL;
  }
  cxpoints=(PyArrayObject*)PyArray_ContiguousFromObject((PyObject*)xpoints,
							PyArray_DOUBLE,2,2);
  if (cxpoints==NULL){
    Py_DECREF(rarray);
    Py_DECREF(cdata);
    return NULL;
  }
  if (isfukunaga){
    /* For each point in the target array, iterate
       through the whole input data set */
    for(i=0;i<Nx;i++){
      sum=0.0;
      xpos=(double*)(cxpoints->data+i*cxpoints->strides[0]);
      rpos=(double*)(rarray->data+i*rarray->strides[0]);
      for(j=0;j<Nd;j++){
        dpos=(double*)(cdata->data+j*cdata->strides[0]);
        vdifference(xpos,dpos,vdif,D,bandwidth);
        matrixtimesvector(Sm1,vdif,wdif,D);
        normdiff=dotprod(wdif,vdif,D);
        sum+=(*f)(normdiff);
      }
      *rpos=sum*K;
    }
  } else {
    /* For each point in the target array, iterate
       through the whole input data set */
    for(i=0;i<Nx;i++){
      /* Python callback to update progress bar */
      if ( i % 10000 == 0 )
      {
          if(my_callback)
          {
            arglist = Py_BuildValue("(ii)", i, Nx);
            result = PyEval_CallObject(my_callback, arglist);
            Py_DECREF(arglist);
          }
      }    
      sum=0.0;
      xpos=(double*)(cxpoints->data+i*cxpoints->strides[0]);
      rpos=(double*)(rarray->data+i*rarray->strides[0]);
      for(j=0;j<Nd;j++){
        dpos=(double*)(cdata->data+j*cdata->strides[0]);
        vdifference(xpos,dpos,vdif,D,bandwidth);
        normdiff=dotprod(vdif,vdif,D);
        sum+=(*f)(normdiff);
      }
      *rpos=sum*K;
    }
  }
  /** Free local objects */
  Py_DECREF(cdata);
  Py_DECREF(cxpoints);
  return PyArray_Return(rarray);
}

static PyObject *MPDFEpanechnikov( PyObject *self , PyObject *args )
{
  PyArrayObject *Xi;
  PyArrayObject *x;
  double bandwidth;
  PyArrayObject *Sm1;
  double sqrtdetS;

  /*
       If the argument corresponding to the inverse of the covariance matrix
       is None, then, use a simple kernel-based estimator, otherwise,
       apply a Fukunaga-type estimator. The selection is done
       inside MPDFEstimator()
  */
  if (!parseMPDFarguments(self,args,&Xi,&x,&bandwidth,&Sm1,&sqrtdetS))
    return NULL;
  else
    return MPDFEstimator(Xi,x,bandwidth,Sm1,sqrtdetS,mpdfepanechnikov);
}


static PyObject *MPDFGaussian( PyObject *self , PyObject *args )
{
  PyArrayObject *Xi;
  PyArrayObject *x;
  double bandwidth;
  PyArrayObject *Sm1;
  double sqrtdetS;

  /*
       If the argument corresponding to the inverse of the covariance matrix
       is None, then, use a simple kernel-based estimator, otherwise,
       apply a Fukunaga-type estimator. Selected inside MPDFEstimator()
  */
  if (!parseMPDFarguments(self,args,&Xi,&x,&bandwidth,&Sm1,&sqrtdetS))
    return NULL;
  else
    return MPDFEstimator(Xi,x,bandwidth,Sm1,sqrtdetS,mpdfmvgaussian);
}


static PyObject *MPDFOptimumBandwidth( PyObject *self , PyObject *args )
{
  PyArrayObject *data;
  PyArrayObject *cdata;
  int N,D;
  double hopt, Cd, Ak;

  if (!PyArg_ParseTuple(args,"O",&data))
    return NULL;
  if (data->nd==1 || data->nd > KPDFMAXDIMS ){
    PyErr_SetString(PyExc_ValueError,"Multivariate Optimum bandwidth estimation:Unimplemented dimensionality");
    return NULL;
  }
  cdata=(PyArrayObject*)PyArray_ContiguousFromObject((PyObject*)data,
						     PyArray_DOUBLE,2,2);
  if (cdata==NULL)
    return NULL;
  D=cdata->nd;
  N=cdata->dimensions[0];
  Cd=__Cd[D-1];
  Ak=8.*(D+4.)*pow(2*sqrt(acos(-1.)),(double)D)/Cd;
  Ak=pow(Ak,(1./(D+4)));
  hopt=Ak*pow((double)N,(-1./(D+4.)));
  Py_DECREF(cdata);
  return Py_BuildValue("d",hopt);
}

static PyObject *MPDF2DGrid2Array( PyObject *self, PyObject *args )
{
  PyObject *arrays[KPDFMAXDIMS];
  PyArrayObject *numpys[KPDFMAXDIMS];
  int i;
  PyObject *retarray;
  int fastest=0;

  /* Initialize the array objects to NULL */
  for(i=0;i<KPDFMAXDIMS;i++){
    arrays[i]=NULL;
    numpys[i]=NULL;
  }

  /* Parse the arguments */
  if (!PyArg_ParseTuple(args,"OO|i",&(arrays[0]),&(arrays[1]),&fastest))
    return NULL;
  /* Get numeric arrays */
  numpys[0]=(PyArrayObject*)PyArray_ContiguousFromObject(arrays[0],
							   PyArray_DOUBLE,
							   1,1);
  if(numpys[0]==NULL)
    return NULL;
  numpys[1]=(PyArrayObject*)PyArray_ContiguousFromObject(arrays[1],
							   PyArray_DOUBLE,
							   1,1);
  if(numpys[1]==NULL){
    Py_DECREF(numpys[0]);
    return NULL;
  }
  retarray=(PyObject*) MPDFGrid2Array(numpys,fastest);
  Py_DECREF(numpys[0]);
  Py_DECREF(numpys[1]);
  return retarray;
}

static PyObject *MPDF3DGrid2Array( PyObject *self, PyObject *args )
{
  PyObject *arrays[KPDFMAXDIMS];
  PyArrayObject *numpys[KPDFMAXDIMS];
  int i;
  PyObject *retarray;
  int fastest=0;

  /* Initialize the array objects to NULL */
  for(i=0;i<KPDFMAXDIMS;i++){
    arrays[i]=NULL;
    numpys[i]=NULL;
  }

  /* Parse the arguments */
  if (!PyArg_ParseTuple(args,"OOO|i",&(arrays[0]),&(arrays[1]),&(arrays[2]),
			&fastest))
    return NULL;
  /* Get numeric arrays */
  numpys[0]=(PyArrayObject*)PyArray_ContiguousFromObject(arrays[0],
							   PyArray_DOUBLE,
							   1,1);
  if(numpys[0]==NULL)
    return NULL;
  numpys[1]=(PyArrayObject*)PyArray_ContiguousFromObject(arrays[1],
							   PyArray_DOUBLE,
							   1,1);
  if(numpys[1]==NULL){
    Py_DECREF(numpys[0]);
    return NULL;
  }
  numpys[2]=(PyArrayObject*)PyArray_ContiguousFromObject(arrays[2],
							   PyArray_DOUBLE,
							   1,1);
  if(numpys[2]==NULL){
    Py_DECREF(numpys[0]);
    Py_DECREF(numpys[1]);
    return NULL;
  }
  retarray=(PyObject*) MPDFGrid2Array(numpys,fastest);
  Py_DECREF(numpys[0]);
  Py_DECREF(numpys[1]);
  Py_DECREF(numpys[2]);
  return retarray;
}

/****
This function iterates and recurses to allow a flexible implementation.
ijk -> Current indexes in the iteration/recursion cycles
IJK -> Limits to the indices for each dimension
nofparams -> Number of existing NumPy arrays (parameters to the input function)
retarray -> NumPy array which will hold the output grid in a linear format
istart -> Dimension to start the recursion from. (0|last)
istep -> Step in dimensions (1-> Lowest varies fastest, 0-> Highest varies
fastest)
*ipos -> Position of the corresponding cell in the output array
 ***/
static void iterateandrecurse(int *ijk, int *IJK,
			      int nofparams, PyArrayObject *retarray,
			      PyArrayObject **numpys, int istart,
			      int istep , int *ipos )
{
  int i;
  int dim;
  double *target;
  double *origin;

  /* Exit conditions */
  if ( (istart==0 && istep==-1) || (istep==1 && istart==nofparams-1) ){
    /* Iterate through the last dimension before exiting, building the
     cells corresponding to this last dimension */
    for(i=0;i<IJK[istart];i++){
      ijk[istart]=i;
      for(dim=0;dim<nofparams;dim++){
	target=(double*) (retarray->data+(*ipos)*retarray->strides[0]+
			  dim*retarray->strides[1]);
	origin=(double*) ((numpys[dim])->data +
			  (numpys[dim])->strides[0]*ijk[dim]);
	*target=*origin;
      }
      *ipos=*ipos+1;
    }
    return;
  }
  for(i=0;i<IJK[istart];i++){
    /*
       Setup the values of the current index corresponding to this dimension
       and iterate and recurse to get the one(s) corresponding to the
       next dimension(s)
    */
    ijk[istart]=i;
    iterateandrecurse(ijk,IJK,nofparams,retarray,numpys,
		      istart+istep,istep,ipos);
  }
}

static PyObject *MPDFGrid2Array( PyArrayObject **numpys , int fastest )
{
  int nofparams=0;
  int i;
  PyArrayObject *retarray;
  int numofcells=1;
  int cell=0;
  npy_intp dimensions[2];
  int ijk[KPDFMAXDIMS];
  int IJK[KPDFMAXDIMS];

  /* Initialize the counters */
  for(i=0;i<KPDFMAXDIMS;i++){
    ijk[i]=0;
    IJK[i]=0;
  }
  /* Get the number of valid arrays */
  for(i=0;i<KPDFMAXDIMS;i++){
    if (numpys[i]==NULL)
      break;
    nofparams=i+1;
    numofcells=numofcells*numpys[i]->dimensions[0];
    IJK[i]=numpys[i]->dimensions[0];
  }
  dimensions[0]=numofcells;
  dimensions[1]=nofparams;
  retarray=(PyArrayObject*)PyArray_SimpleNew(2,dimensions,PyArray_DOUBLE);
  if(retarray==NULL)
    return NULL;

  if (fastest==0){
    /* Last dimension varies fastest */
    iterateandrecurse(ijk,IJK,nofparams,retarray,numpys,0,1,&cell);
  } else {
    /* First dimension varies fastest */
    iterateandrecurse(ijk,IJK,nofparams,retarray,numpys,nofparams-1,-1,&cell);
  }
  return (PyObject*)PyArray_Return(retarray);
}

