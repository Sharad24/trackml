/*
########################################################################
# ======================  TrackML CHALLENGE MODEL  =====================
########################################################################
# Author: Isabelle Guyon, Victor Estrade
# Date: Apr 10, 2018

# ALL INFORMATION, SOFTWARE, DOCUMENTATION, AND DATA ARE PROVIDED "AS-IS".
# PARIS-SUD UNIVERSITY, THE ORGANIZERS OR CODE AUTHORS DISCLAIM
# ANY EXPRESSED OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR ANY PARTICULAR PURPOSE, AND THE
# WARRANTY OF NON-INFRIGEMENT OF ANY THIRD PARTY'S INTELLECTUAL PROPERTY RIGHTS.
# IN NO EVENT SHALL PARIS-SUD UNIVERSITY AND/OR OTHER ORGANIZERS BE LIABLE FOR ANY SPECIAL,
# INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER ARISING OUT OF OR IN
# CONNECTION WITH THE USE OR PERFORMANCE OF SOFTWARE, DOCUMENTS, MATERIALS,
# PUBLICATIONS, OR INFORMATION MADE AVAILABLE FOR THE CHALLENGE.
*/
   
// READ THE DOC !
// https://docs.scipy.org/doc/numpy/user/c-info.how-to-extend.html
// https://docs.python.org/3.5/c-api/arg.html

#include <math.h>
#include "Python.h"
#include "numpy/arrayobject.h"


void SGTrackerProcessEvent( float *x, float *y, float *z, int *id, int *vol, int *layer, int *labels, int nHits);

void SGTrackerInit( const char * );

 
static PyObject *
trackerInterfaceInit(PyObject *dummy, PyObject *args)
{  
  char *path;

   if (!PyArg_ParseTuple(args, "s",  &path)) {
     return NULL;
   }
   SGTrackerInit( path );
   Py_INCREF(Py_None);
   return Py_None;
}


static PyObject *
trackerInterfaceExec(PyObject *dummy, PyObject *args)
{
  //printf("SG: run trackerInterface().. \n");
   // trackerInterface() expects 3 arguments
  
  PyObject *arg[7] = {0,0,0,0, 0,0,0};

  if (!PyArg_ParseTuple(args, "O!O!O!O!O!O!O!",
                        &PyArray_Type, &arg[0], // Extract arg1 and check that it is an array-like object
                        &PyArray_Type, &arg[1], 
			&PyArray_Type, &arg[2], 
                        &PyArray_Type, &arg[3], 
			&PyArray_Type, &arg[4], 
			&PyArray_Type, &arg[5], 
			&PyArray_Type, &arg[6]
			)      
      ) // Extract arg3 and check that it is an array-like object
    return NULL;

  int nHits = 0;

  PyArrayObject *arr[7] = {0,0,0,0, 0,0,0};
  void *array[7] = {0,0,0,0, 0,0,0};

  for( int i=0; i<7; i++ ){
    // Convert arguments to PyArrayObject

#if NPY_API_VERSION >= 0x0000000c
    if( i<3 ){
      arr[i] = (PyArrayObject*) PyArray_FROM_OTF(arg[i], NPY_FLOAT32, NPY_ARRAY_INOUT_ARRAY2);
    } else {
      arr[i] = (PyArrayObject*) PyArray_FROM_OTF(arg[i], NPY_INT32,   NPY_ARRAY_INOUT_ARRAY2);
    }
#else
    if( i<3 ){
      arr[i] = (PyArrayObject*) PyArray_FROM_OTF(arg[i], NPY_FLOAT32, NPY_ARRAY_INOUT_ARRAY);
    } else {
      arr[i] = (PyArrayObject*) PyArray_FROM_OTF(arg[i], NPY_INT32,   NPY_ARRAY_INOUT_ARRAY);
    }    
#endif
    if (arr[i] == NULL) goto fail;
    array[i] = (void*) PyArray_DATA(arr[i]);

    int nd =  PyArray_NDIM(arr[i]);  //-- number of dimensions
    if (nd > 1) goto fail;
    npy_intp * arr_shape =  PyArray_SHAPE(arr[i]); 
    nHits = arr_shape[0];
  }

     
  // Access to the data of the numpy array :
  float *x = (float *) array[0];
  float *y = (float *) array[1];
  float *z = (float *) array[2];

  int *id     = (int*) array[3];
  int *vol    = (int*) array[4];
  int *layer  = (int*) array[5];
  int *labels = (int*) array[6];

  SGTrackerProcessEvent( x, y, z, id, vol, layer, labels, nHits );

// Resolve copy and reference counting :

#if NPY_API_VERSION >= 0x0000000c
  for( int i=0; i<7; i++ ){
    PyArray_ResolveWritebackIfCopy(arr[i]);
  }
#endif
  for( int i=0; i<7; i++ ){
    Py_DECREF(arr[i]);
  }

  Py_INCREF(Py_None);
  return Py_None;
  
// If something failed this section carefully clean the arrays
 fail:
  
#if NPY_API_VERSION >= 0x0000000c
   for( int i=0; i<7; i++ ){
     PyArray_DiscardWritebackIfCopy(arr[i]);
   }
#endif
   for( int i=0; i<7; i++ ){
     Py_XDECREF(arr[i]);
   }
  
   return NULL;
}

// Define the Python/C methods
static PyMethodDef PythonInterfaceMethods[] = {
    { "trackerInterfaceInit", trackerInterfaceInit, METH_VARARGS, "Doc string"},
    { "trackerInterfaceExec", trackerInterfaceExec, METH_VARARGS, "Doc string"},
    {NULL, NULL, 0, NULL} /* Sentinel */
};


// Define the module
static struct PyModuleDef trackerInterfaceModule = {
   PyModuleDef_HEAD_INIT,
   "trackerInterface",   /* name of module */
   NULL, /* module documentation, may be NULL */
   -1,       /* size of per-interpreter state of the module,
                or -1 if the module keeps state in global variables. */
   PythonInterfaceMethods
};

// Define the initialization
PyMODINIT_FUNC
PyInit_trackerInterface(void)
{
    import_array();
    return PyModule_Create(&trackerInterfaceModule);
}
