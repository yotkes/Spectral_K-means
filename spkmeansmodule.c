#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "spkmeans.h"

static PyObject* fit(PyObject *self, PyObject *args);
static PyObject* kmeans_pp(PyObject *self, PyObject *args);
static PyObject* Convert_Big_Array(double **array, int N, int d);
PyMODINIT_FUNC PyInit_spkmeansmodule(void);

/* This function fits each goal from Python to its specific output (matrix) - received from C */
static PyObject* fit(PyObject *self, PyObject *args){
    /* Variables Decleration */
    int k;
    int N;
    int d;
    int goal;
    double** obs_c;
    double** c_list;
    PyObject *obs_p;
    PyObject *all;
    PyObject *item;
    PyObject *py_list;
    int i,j;

    /* Checks if input from Python is valid */
    if(!PyArg_ParseTuple(args, "O:C function getting arguments", &all)) {
        return NULL;
    }
    if (!PyList_Check(all)){
        return NULL;
    }

    /* Inputs which are passed from Python */
    k = (int)PyLong_AsLong(PyList_GetItem(all, 0));
    N = (int)PyLong_AsLong(PyList_GetItem(all, 1));
    d = (int)PyLong_AsLong(PyList_GetItem(all, 2));
    goal = (int)PyLong_AsLong(PyList_GetItem(all, 3));
    obs_p = PyList_GetItem(all, 4);

    /* Checks if obs_p is a list object */
    if(!PyList_Check(obs_p)){
        return NULL;
    }

    /* Memory allocation - observations matrix */
    obs_c = (double**)malloc(N*sizeof(double*));
    if (obs_c == NULL) 
    {
        printf("%s\n","An Error Has Occured");
        exit(1);
    }

    for(i = 0; i<N; i++)
    {
        obs_c[i] = (double*)malloc(d*sizeof(double));
        if (obs_c[i] == NULL) 
        {
            printf("%s\n","An Error Has Occured");
            exit(1);
        }
    }

    /* Inserting matrix's elements */
    for(i = 0; i<N; i++){
        item = PyList_GetItem(obs_p, i);
        if (!PyList_Check(item)){
            return NULL;
        }
        for(j = 0; j<d; j++){
            obs_c[i][j] = PyFloat_AsDouble(PyList_GetItem(item, j));
        }
    }

    c_list = spk_linker(&k, goal, obs_c, N, d); /* Recieves a matrix defined by goal's value */
    
    /* Free Memory */
    for(i=0;i<N;i++)
    {
        free(obs_c[i]);
    }
    free(obs_c);
    
    if(goal == 0) /* goal == spk */
    {
        py_list = Convert_Big_Array(c_list, N, k); /* Conversion of C matrix to a Python list */
        for(i=0;i<N;i++)
        {
            free(c_list[i]);
        }
    }
    else if(goal == 4) /* goal == jacobi */
    {
        py_list = Convert_Big_Array(c_list, N+1, N); /* Conversion of C matrix to a Python list */
        for(i=0;i<N+1;i++)
        {
            free(c_list[i]);
        }
    }
    else  /* All the other goals */
    {
        py_list = Convert_Big_Array(c_list, N, N); /* Conversion of C matrix to a Python list */
        for(i=0;i<N;i++)
        {
            free(c_list[i]);
        }
    }

    free(c_list);
    return py_list; /* Returns a Python list */
}


/* K-means ++ Algorithm - creating the k-clusters matrix which is sent from C back to Python */
static PyObject* kmeans_pp(PyObject *self, PyObject *args){
    /* Variables Decleration */
    int k;
    int N;
    int* first_c;
    double** obs_c;
    double** c_list;
    PyObject *first_p;
    PyObject *obs_p;
    PyObject *all;
    PyObject *item;
    PyObject *py_list;
    int i;
    int j;
    
    /* Checks if input from Python is valid */
    if(!PyArg_ParseTuple(args, "O:C function getting arguments", &all)) {
        return NULL;
    }
    if (!PyList_Check(all)){
        return NULL;
    }

    /* Inputs which are passed from Python */
    k = (int)PyLong_AsLong(PyList_GetItem(all, 0));
    N = (int)PyLong_AsLong(PyList_GetItem(all, 1));
    first_p = PyList_GetItem(all, 2);
    obs_p = PyList_GetItem(all, 3);
    
    /* Checks if first_p is a list object */
    if (!PyList_Check(first_p)){
        return NULL;
    }

    /* Memory allocation - k initial centroids */
    first_c = calloc(k, sizeof(int));
    if (first_c == NULL) 
    {
        printf("%s\n","An Error Has Occured");
        exit(1);
    }

    /* Inserting initial centroids' indexes */
    for (i = 0; i<k; i++){
        first_c[i] = (int)PyLong_AsLong(PyList_GetItem(first_p, i));
    }

    /* Checks if obs_p is a list object */
    if (!PyList_Check(obs_p)){
        return NULL;
    }

    /* Memory allocation - observations matrix (in this case: T matrix) */
    obs_c = (double**)malloc(N*sizeof(double*));
    if (obs_c == NULL) 
    {
        printf("%s\n","An Error Has Occured");
        exit(1);
    }

    for(i = 0; i<N; i++)
    {
        obs_c[i] = (double*)malloc(k*sizeof(double));
        if (obs_c[i] == NULL) 
        {
            printf("%s\n","An Error Has Occured");
            exit(1);
        }
    }

    /* Inserting matrix's elements */
    for (i = 0; i<N; i++){
        item = PyList_GetItem(obs_p, i);

        /* Checks if item is a list object */
        if (!PyList_Check(item)){
            return NULL;
        }
        for (j = 0; j<k; j++){
            obs_c[i][j] = PyFloat_AsDouble(PyList_GetItem(item, j));
        }
    }

    c_list = k_means(k, N, k, first_c, obs_c); /* Receives k-clusters matrix */

    /* Free Memory */
    free(first_c);
    for(i=0;i<N;i++)
    {
        free(obs_c[i]);
    }
    free(obs_c);

    py_list = Convert_Big_Array(c_list, k, k); /* Conversion of C matrix into a Python list */

    /* Free Memory */
    for(i=0;i<k;i++)
    {
        free(c_list[i]);
    }
    free(c_list);

    return py_list; /* Returns a Python list */
}


/* This function converts a matrix from C to a Python list */
static PyObject *Convert_Big_Array(double **array, int N, int d){
    PyObject *py_list_of_lists, *py_list, *item;
    int i;
    int j;
    py_list_of_lists = PyList_New(N);
    for(i = 0; i < N; i++) {
        py_list = PyList_New(d);
        for(j = 0; j < d; j++) {
            item = PyFloat_FromDouble(array[i][j]);
            PyList_SetItem(py_list, j, item);
        }
        PyList_SetItem(py_list_of_lists, i, py_list);
    }
    return py_list_of_lists;
}

static PyMethodDef capiMethods[] = {
    {"fit",                   /* the Python method name that will be used */
      (PyCFunction) fit, /* the C-function that implements the Python function and returns static PyObject*  */
      METH_VARARGS,           /* flags indicating parametersaccepted for this function */
      PyDoc_STR("main flow of all the spkmeans procedures")}, /*  The docstring for the function */
    {"kmeans_pp",                   /* the Python method name that will be used */
      (PyCFunction) kmeans_pp, /* the C-function that implements the Python function and returns static PyObject*  */
      METH_VARARGS,           /* flags indicating parametersaccepted for this function */
      PyDoc_STR("spkmeans calculator (based on kmeans ++ algorithm)")}, /*  The docstring for the function */  
    {NULL, NULL, 0, NULL}     /* The last entry must be all NULL as shown to act as a
                                 sentinel. Python looks for this entry to know that all
                                 of the functions for the module have been defined. */
};

/* This initiates the module using the above definitions. */
static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "spkmeansmodule", /* name of module */
    NULL, /* module documentation, may be NULL */
    -1,  /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
    capiMethods /* the PyMethodDef array from before containing the methods of the extension */
};

PyMODINIT_FUNC PyInit_spkmeansmodule(void)
{
    PyObject *m;
    m = PyModule_Create(&moduledef);
    if (!m) {
        return NULL;
    }
    return m;
}