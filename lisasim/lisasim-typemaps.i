%include typemaps.i

// The following from Scott Ransom's numpy.i (ransom@cfa.harvard.edu)
// needs arrayobject.h from the Numeric distribution

%{
#include "arrayobject.h"

#define ISCONTIGUOUS(m) ((m)->flags & CONTIGUOUS)

#define PyArray_CONTIGUOUS(m) (ISCONTIGUOUS(m) ? Py_INCREF(m), m : \
(PyArrayObject *)(PyArray_ContiguousFromObject((PyObject *)(m), (m)->descr->type_num, 0,0)))
%}

// Map a Numeric array into an array of doubles
// Modified to work with numarray, as well (?)

// I suppose we could go to making this generic for double*, assuming that's different from double[]

%typemap(python, in) double* IN_1D_DOUBLE {
  PyArrayObject *arr;
  
  /* Check that obj is really a 1D array of bytes */
  
  if (!PyArray_Check($input)) {
    PyErr_SetString(PyExc_TypeError,"First argument is not an array");
    return NULL;
  }
  
  /* check type (could also use arr->descr->type_num) */
  
  if (arr->descr->type_num != PyArray_DOUBLE) {
    PyErr_SetString(PyExc_TypeError, "Incorrect array type: we need an array of DOUBLE");
    return NULL;
  }

  arr = PyArray_CONTIGUOUS((PyArrayObject *)$input);

  if (arr->nd != 1) { /* we are really strict ! */
    PyErr_SetString(PyExc_TypeError, "Incorrect number of dims: we want a 1D array");
    return NULL;
  }

  $1 = (double *)arr->data;

  /* Release our local copy of the PyArray */

  Py_DECREF(arr);
}

// The following modified from the SWIG documentation
// Map a Python sequence into any sized C double array

%typemap(python, in) double[ANY](double temp[$1_dim0]) {

  int i;

  // check that we are really getting a sequence (list or tuple)

  if (!PySequence_Check($input)) {
      PyErr_SetString(PyExc_TypeError,"Expecting a sequence");
      return NULL;
  }

  // check that the length is as requested

  if (PyObject_Length($input) != $1_dim0) {
      PyErr_SetString(PyExc_ValueError,"Expecting a sequence with $1_dim0 elements");
      return NULL;
  }

  // convert each element

  for (i =0; i < $1_dim0; i++) {
      PyObject *o = PySequence_GetItem($input,i);

      if(PyFloat_Check(o)) {
         temp[i] = PyFloat_AsDouble(o);
      } else if(PyInt_Check(o)) {
         temp[i] = PyInt_AsLong(o);
      } else {	    
         PyErr_SetString(PyExc_ValueError,"Expecting a sequence of floats");
         return NULL;
      }
  }

  // return pointer the the array

  $1 = &temp[0];
} 

// The following modified from the SWIG documentation
// Map a Python sequence of noise objects into an array of pointers

%typemap(python, in) Noise *[ANY](Noise *temp[$1_dim0]) {
  int i;

  // check that we are really getting a sequence (list or tuple)

  if (!PySequence_Check($input)) {
      PyErr_SetString(PyExc_TypeError,"Expecting a sequence");
      return NULL;
  }

  // check that the length is as requested

  if (PyObject_Length($input) != $1_dim0) {
      PyErr_SetString(PyExc_ValueError,"Expecting a sequence with $1_dim0 elements");
      return NULL;
  }

  // convert each element

  for (i =0; i < $1_dim0; i++) {
      PyObject *o = PySequence_GetItem($input,i);
      
      SWIG_ConvertPtr(o, (void **)&temp[i], $descriptor(Noise *), SWIG_POINTER_EXCEPTION);
  }

  // return pointer the the array

  $1 = &temp[0];
}
