/* $Id$
 * $Date$
 * $Author$
 * $Revision$
 */

%include typemaps.i

// The following from Scott Ransom's numpy.i (ransom@cfa.harvard.edu)
// needs arrayobject.h from the Numeric distribution

%{
#include "numpy/arrayobject.h"

// this was for Numeric
// #define ISCONTIGUOUS(m) ((m)->flags & CONTIGUOUS)

// for numpy
#define ISCONTIGUOUS(m) (PyArray_ISCONTIGUOUS(m))

// if the Numeric array is not contiguous, create a contiguous copy;
// if it is contiguous, use it after increasing its reference count

#define PyArray_CONTIGUOUS(m) (ISCONTIGUOUS(m) ? Py_INCREF(m), m : \
(PyArrayObject *)(PyArray_ContiguousFromObject((PyObject *)(m), (m)->descr->type_num, 0,0)))
%}

%init %{
  import_array();
%}

%typemap(in) double* NUMPY_ARRAY_DOUBLE {
  PyArrayObject *arr;
  
  /* Check that obj is really a 1D array of bytes */
  
  if (!PyArray_Check($input)) {
    PyErr_SetString(PyExc_TypeError,"First argument is not an array");
    return NULL;
  }
  
  /* check type (could also use arr->descr->type_num) */
  
  if (PyArray_ObjectType($input,0) != PyArray_DOUBLE) {
    PyErr_SetString(PyExc_TypeError, \
		    "Incorrect array type: we need an array of DOUBLE");
    return NULL;
  }
  arr = PyArray_CONTIGUOUS((PyArrayObject *)$input);
  if (arr->nd != 1) { /* we are really strict ! */
    PyErr_SetString(PyExc_TypeError, \
		    "Incorrect number of dims: we want a 1D array");
    return NULL;
  }
  $1 = (double *)arr->data;

  // this seems a bit strange; why DECREF it after it's been INCREFed
  // in the PyArray_CONTIGUOUS above? And won't this destroy an
  // eventual contiguous local copy?

  Py_DECREF(arr);  /* Release our local copy of the PyArray */
}

// The following modified from the SWIG documentation
// Map a Python sequence into any sized C double array

%typemap(in) double PYTHON_SEQUENCE_DOUBLE[ANY] (double temp[$1_dim0]) {
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

  // return pointer to the array

  $1 = &temp[0];
} 

// Map a Numeric array into an array of doubles, pass also the number of elements

%typemap(in) (double* numarray, long length) {
	PyArrayObject *arr;
	
	/* Check that obj is really an array */
	
	if (!PyArray_Check($input)) {
		PyErr_SetString(PyExc_TypeError,"First argument is not an array");
		return NULL;
	}
	
	/* check type (could also use arr->descr->type_num) */
	
	if (PyArray_ObjectType($input,0) != PyArray_DOUBLE) {
		PyErr_SetString(PyExc_TypeError, \
			"Incorrect array type: we need an array of DOUBLE");

		return NULL;
	}
	
	arr = PyArray_CONTIGUOUS((PyArrayObject *)$input);

	if (arr->nd == 1) {
		$1 = (double *)arr->data;
		$2 = (long)arr->dimensions[0];
	} else if (arr->nd == 2) {
		/* flatten data for 2-dimensional array */
	
		$1 = (double *)arr->data;
		$2 = (long)arr->dimensions[0] * (long)arr->dimensions[1];
	} else {
		PyErr_SetString(PyExc_TypeError, \
			"Incorrect number of dims: we want a 1D or 2D array");

		return NULL;
	}
	
	// this seems a bit strange; why DECREF it after it's been INCREFed
	// in the PyArray_CONTIGUOUS above? And won't this destroy an
	// eventual contiguous local copy?
	
	Py_DECREF(arr);  /* Release our local copy of the PyArray */
}

// Map a Python sequence into an array of doubles;
// pass also the number of elements

%typemap(in) (double *doublearray, int doublenum) {
	int i;
	
	// check that we are really getting a sequence (list or tuple)
	
	if (!PySequence_Check($input)) {
		PyErr_SetString(PyExc_TypeError,"Expecting a sequence");
		return NULL;
	}
	
	int dim = PySequence_Size($input);
	double *temp = new double[dim];
	
	// convert each element
	
	for (i = 0; i < dim; i++) {
		PyObject *o = PySequence_GetItem($input,i);
		temp[i] = PyFloat_AsDouble(o);
	}
	
	// return pointer to the array
	
	$1 = temp;
	$2 = dim;}

%typemap(freearg) (double *doublearray, int doublenum)  {
   delete [] $1;
}

%typemap(in) (Signal **thesignals, int signals) {
  int i;

  // check that we are really getting a sequence (list or tuple)

  if (!PySequence_Check($input)) {
      PyErr_SetString(PyExc_TypeError,"Expecting a sequence");
      return NULL;
  }

  int dim = PySequence_Size($input);
  Signal **temp = new Signal*[dim];
  
  // convert each element

  for (i=0;i<dim;i++) {
      PyObject *o = PySequence_GetItem($input,i);
      
      SWIG_ConvertPtr(o, (void **)&temp[i], $descriptor(Signal *), SWIG_POINTER_EXCEPTION);
  }

  // return pointer the the array

  $1 = temp;
  $2 = dim;
}

%typemap(freearg) (Signal **thesignals, int signals)  {
   delete [] $1;
}

// The following modified from the SWIG documentation
// Map a Python sequence of noise objects into an array of pointers

%typemap(in) Noise *PYTHON_SEQUENCE_NOISE[ANY] (Noise *temp[$1_dim0]) {
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

      // hack! apparently the ConvertPtr above adds a Python ref. But why?
      Py_DECREF(o);
  }

  // return pointer the the array

  $1 = &temp[0];
}

// convert a list of Wave objects; could probably use a temp variable
// (as above) instead of the freearg typemap

%typemap(in) (Wave **WaveSeq, int WaveNum) {
  int i;

  // check that we are really getting a sequence (list or tuple)

  if (!PySequence_Check($input)) {
      PyErr_SetString(PyExc_TypeError,"Expecting a sequence");
      return NULL;
  }

  int dim = PySequence_Size($input);
  Wave **temp = new Wave*[dim];

  // convert each element

  for (i = 0; i < dim; i++) {
      PyObject *o = PySequence_GetItem($input,i);
      
      SWIG_ConvertPtr(o, (void **)&temp[i], $descriptor(Wave *), SWIG_POINTER_EXCEPTION);
  }

  // return pointer to the array
  
  $1 = temp;
  $2 = dim;
}

%typemap(freearg) (Wave **WaveSeq, int WaveNum)  {
   delete [] $1;
}

// from the SWIG documentation: input a python function

%typemap(in) PyObject* PYTHONFUNC {
  if (!PyCallable_Check($input)) {
      PyErr_SetString(PyExc_TypeError, "Need a callable object!");
      return NULL;
  }
  $1 = $input;
}

// typechecks

%typecheck(SWIG_TYPECHECK_POINTER) double *NUMPY_ARRAY_DOUBLE {
    /* Check that obj is really an array (of something) */
  
    if (!PyArray_Check($input))
	$1 = 0;
    else
	$1 = 1;
}

%typecheck(SWIG_TYPECHECK_POINTER) double *numarray {
    /* Check that obj is really an array (of something) */
  
    if (!PyArray_Check($input))
	$1 = 0;
    else
	$1 = 1;
}


%typecheck(SWIG_TYPECHECK_POINTER) double PYTHON_SEQUENCE_DOUBLE[ANY] {
    if (!PySequence_Check($input)) {
        $1 = 0;
    } else {
	PyObject *o = PySequence_GetItem($input,0);

	if(PyFloat_Check(o) || PyInt_Check(o))
	    $1 = 1;
	else
	    $1 = 0;
    }
}

%typecheck(SWIG_TYPECHECK_POINTER) Noise *PYTHON_SEQUENCE_NOISE[ANY] {
    if (!PySequence_Check($input)) {
	$1 = 0;
    } else {
	PyObject *o = PySequence_GetItem($input,0);
	void *mypointer;

	if(SWIG_ConvertPtr(o, (void **)&mypointer, $descriptor(Noise *), SWIG_POINTER_EXCEPTION) != -1)
	    $1 = 1;
	else {
	    PyErr_Clear();
	    $1 = 0;
	}
    }
}

// return a C++ Vector as a tuple

// This tells SWIG to treat a Vector & argument with name outvector as
// an output value.  We'll append the value to the current result which 
// is guaranteed to be a List object by SWIG.

%typemap(in,numinputs=0) Vector &outvector {
    Vector *a = new Vector();
    $1 = a;
}

%typemap(argout) Vector &outvector {
    PyObject *t;

    t = PyTuple_New(3);

    PyObject *p0, *p1, *p2;

    p0 = PyFloat_FromDouble((*$1)[0]);
    p1 = PyFloat_FromDouble((*$1)[1]);
    p2 = PyFloat_FromDouble((*$1)[2]);

    PyTuple_SetItem(t,0,p0);
    PyTuple_SetItem(t,1,p1);
    PyTuple_SetItem(t,2,p2);

    $result = t;
}

%typemap(freearg) Vector &outvector {
    delete $1;
}

%typemap(out) Vector {
    PyObject *t;

    t = PyTuple_New(3);

    PyObject *p0, *p1, *p2;

    p0 = PyFloat_FromDouble($1[0]);
    p1 = PyFloat_FromDouble($1[1]);
    p2 = PyFloat_FromDouble($1[2]);

    PyTuple_SetItem(t,0,p0);
    PyTuple_SetItem(t,1,p1);
    PyTuple_SetItem(t,2,p2);

    $result = t;
}

// return a C++ Tensor as a tuple of tuples

%typemap(in,numinputs=0) Tensor &outtensor {
    Tensor *a = new Tensor();
    $1 = a;
}

%typemap(argout) Tensor &outtensor {
    PyObject *t;

    t = PyTuple_New(3);

    PyObject *p0, *p1, *p2;

    p0 = PyTuple_New(3);
    p1 = PyTuple_New(3);
    p2 = PyTuple_New(3);

    PyTuple_SetItem(p0,0,PyFloat_FromDouble((*$1)[0][0]));
    PyTuple_SetItem(p0,1,PyFloat_FromDouble((*$1)[0][1]));
    PyTuple_SetItem(p0,2,PyFloat_FromDouble((*$1)[0][2]));

    PyTuple_SetItem(p1,0,PyFloat_FromDouble((*$1)[1][0]));
    PyTuple_SetItem(p1,1,PyFloat_FromDouble((*$1)[1][1]));
    PyTuple_SetItem(p1,2,PyFloat_FromDouble((*$1)[1][2]));

    PyTuple_SetItem(p2,0,PyFloat_FromDouble((*$1)[2][0]));
    PyTuple_SetItem(p2,1,PyFloat_FromDouble((*$1)[2][1]));
    PyTuple_SetItem(p2,2,PyFloat_FromDouble((*$1)[2][2]));

    PyTuple_SetItem(t,0,p0);
    PyTuple_SetItem(t,1,p1);
    PyTuple_SetItem(t,2,p2);

    $result = t;
}

%typemap(out) Tensor {
    PyObject *t;

    t = PyTuple_New(3);

    PyObject *p0, *p1, *p2;

    p0 = PyTuple_New(3);
    p1 = PyTuple_New(3);
    p2 = PyTuple_New(3);

    PyTuple_SetItem(p0,0,PyFloat_FromDouble($1[0][0]));
    PyTuple_SetItem(p0,1,PyFloat_FromDouble($1[0][1]));
    PyTuple_SetItem(p0,2,PyFloat_FromDouble($1[0][2]));

    PyTuple_SetItem(p1,0,PyFloat_FromDouble($1[1][0]));
    PyTuple_SetItem(p1,1,PyFloat_FromDouble($1[1][1]));
    PyTuple_SetItem(p1,2,PyFloat_FromDouble($1[1][2]));

    PyTuple_SetItem(p2,0,PyFloat_FromDouble($1[2][0]));
    PyTuple_SetItem(p2,1,PyFloat_FromDouble($1[2][1]));
    PyTuple_SetItem(p2,2,PyFloat_FromDouble($1[2][2]));

    PyTuple_SetItem(t,0,p0);
    PyTuple_SetItem(t,1,p1);
    PyTuple_SetItem(t,2,p2);

    $result = t;
}

%typemap(freearg) Tensor &outtensor {
    delete $1;
}

%typemap(out) Noise *[ANY] {
    PyObject *t;

    t = PyTuple_New($1_dim0);

    for(int i=0;i<$1_dim0;i++)
        PyTuple_SetItem(t,i,SWIG_NewPointerObj($1[i],$descriptor(Noise *),0));

    $result = t;
}
