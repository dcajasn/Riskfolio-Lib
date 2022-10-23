// arma_numpy.i
// A SWIG interface file to convert between Armadillo vectors/matrices and Numpy arrays.

%header %{
#include <array>
#include <algorithm>
#include <type_traits>
#include <complex>
#include <armadillo>
#include <stdexcept>
#include <numpy/arrayobject.h>
%}

%include "numpy.i"
%include "std_complex.i"

%init %{
    import_array(); // Required: Imports numpy c-api.
%}

// Define some helpers for converting types to/from typecodes or getting info from a type.
// Define also some enable_if helpers for building arma objects or returning their shapes.

%fragment("arma_numpy_utilities","header"){

#include <armadillo>

    // type to numpy typecode
    template<class T> struct to_typecode {};
    template<> struct to_typecode<int>                  { static const int value = NPY_INT; };
    template<> struct to_typecode<unsigned>             { static const int value = NPY_UINT; };
    template<> struct to_typecode<long long>            { static const int value = NPY_INT64; };
    template<> struct to_typecode<unsigned long long>   { static const int value = NPY_UINT64; };
    template<> struct to_typecode<float>                { static const int value = NPY_FLOAT; };
    template<> struct to_typecode<double>               { static const int value = NPY_DOUBLE; };
    template<> struct to_typecode<std::complex<float>>  { static const int value = NPY_COMPLEX64; };
    template<> struct to_typecode<std::complex<double>> { static const int value = NPY_COMPLEX128; };

    // numpy typecode to type
    template<int TYPE> struct from_typecode {};
    template<> struct from_typecode<NPY_INT>        {  using type = int; };
    template<> struct from_typecode<NPY_UINT>       {  using type = unsigned; };
    template<> struct from_typecode<NPY_INT64>      {  using type = long long; };
    template<> struct from_typecode<NPY_UINT64>     {  using type = unsigned long long; };
    template<> struct from_typecode<NPY_FLOAT>      {  using type = float; };
    template<> struct from_typecode<NPY_DOUBLE>     {  using type = double; };
    template<> struct from_typecode<NPY_COMPLEX64>  {  using type = arma::cx_float; };
    template<> struct from_typecode<NPY_COMPLEX128> {  using type = arma::cx_double; };

    // info from arma type
    template<class T>
    struct arma_info {
        // element type
        using element_t = typename T::elem_type;
        // corresponding typecode
        static constexpr int typecode = to_typecode<element_t>::value;
        // how many dims?
        static constexpr int dims = ( std::is_same<T,arma::Col<element_t>>::value || std::is_same<T,arma::Row<element_t>>::value ? 1 :
                                        ( std::is_same<T,arma::Mat<element_t>>::value ? 2 :
                                            ( std::is_same<T,arma::Cube<element_t>>::value ? 3 : 0 )
                                        )
                                    );
        static_assert( dims, "Cannot wrap Armadillo objects with more dims than Cube");
    };

    // arma_numpy exception
    class ArmaNumpyException : public std::runtime_error {
        using std::runtime_error::runtime_error;
    };

    // build arma object given memptr and dims
    template<class T, typename std::enable_if<arma_info<T>::dims==1,bool>::type = true>
    T* arma_from_ptr( typename arma_info<T>::element_t* data, npy_intp* dims, int ndims) {
        switch(ndims){
            case 1: return new T(data,dims[0],false);
            default: throw ArmaNumpyException("Can't cast Numpy array of len(shape) > 1 to an Armadillo vector");
        }
    }

    template<class T, typename std::enable_if<arma_info<T>::dims==2,bool>::type = true>
    T* arma_from_ptr( typename arma_info<T>::element_t* data, npy_intp* dims, int ndims) {
        switch(ndims){
            case 1: return new T(data,dims[0],1,false);
            case 2: return new T(data,dims[0],dims[1],false);
            default: throw ArmaNumpyException("Can't cast Numpy array of len(shape) > 2 to an Armadillo matrix");
        }
    }

    template<class T, typename std::enable_if<arma_info<T>::dims==3,bool>::type = true>
    T* arma_from_ptr( typename arma_info<T>::element_t* data, npy_intp* dims, int ndims) {
        switch(ndims){
            case 1: return new T(data,dims[0],1,1,false);
            case 2: return new T(data,dims[0],dims[1],1,false);
            case 3: return new T(data,dims[0],dims[1],dims[2],false);
            default: throw ArmaNumpyException("Can't cast Numpy array of len(shape) > 3 to an Armadillo cube");
        }
    }

    // get shape of arma object
    template<class T, typename std::enable_if<arma_info<T>::dims==1,bool>::type = true>
    auto arma_shape( const T* t) -> typename std::array<npy_intp,arma_info<T>::dims> {
        return std::array<npy_intp,arma_info<T>::dims>{ (npy_intp)t->n_elem};
    }

    template<class T, typename std::enable_if<arma_info<T>::dims==2,bool>::type = true>
    auto arma_shape( const T* t) -> typename std::array<npy_intp,arma_info<T>::dims> {
        return std::array<npy_intp,arma_info<T>::dims>{ (npy_intp)t->n_rows, (npy_intp)t->n_cols};
    }

    template<class T, typename std::enable_if<arma_info<T>::dims==3,bool>::type = true>
    auto arma_shape( const T* t) -> typename std::array<npy_intp,arma_info<T>::dims> {
        return std::array<npy_intp,arma_info<T>::dims>{ (npy_intp)t->n_rows, (npy_intp)t->n_cols, (npy_intp)t->n_slices};
    }
}
 
// Define Python wrapper objects for Arma objects

%fragment("arma_wrapper","header"){

    // Define basic Python object that holds an Arma object
    template<class T>
    struct PyArmaObject {
        PyObject_HEAD
        T* data;
    };

    // Define dealloc method
    template<class T>
    static void PyArmaObject_dealloc( PyArmaObject<T>* self){
        delete self->data;
        Py_TYPE(self)->tp_free((PyObject*)self);
    }

    // Define type object
    template<class T>
    struct PyArmaTypeObject {
        static PyTypeObject value; 
    };

    template<class T>
    PyTypeObject PyArmaTypeObject<T>::value = {
        PyObject_HEAD_INIT(NULL)
        .tp_name = "PyArmaObject",
        .tp_basicsize = sizeof(PyArmaObject<T>),
        .tp_dealloc = (destructor) PyArmaObject_dealloc<T>,
        .tp_flags = Py_TPFLAGS_DEFAULT,
        .tp_doc = "arma_numpy wrapper for Armadillo objects",
        .tp_alloc = PyType_GenericAlloc,
        .tp_new = PyType_GenericNew,
    };
}

// Define methods for converting between np.array and arma

%fragment("arma_numpy", "header", fragment="NumPy_Fragments", fragment="arma_numpy_utilities", fragment="arma_wrapper"){

    /* Cast check
     * Tests if the input type can be safely cast to the required type.
     * If 'strict', the input type must exactly match the required type.
     */
    template<class T>
    int arma_numpy_castcheck( PyObject* input, bool strict=false){
        constexpr int typecode = arma_info<T>::typecode;

        PyArrayObject* array = (PyArrayObject*)input;
        int input_typecode = array_type(array);
        
        if( strict ){
            return PyArray_EquivTypenums( input_typecode, typecode);
        } else {
            return PyArray_CanCastSafely( input_typecode, typecode);
        }
    }

    /* Dims check
     * Tests if the input array has compatible dimensions with the required type.
     * If 'strict', we require len(array.shape) == required dims.
     * Otherwise, we will allow upcasting. Downcasting (cube->matrix->vector) is not permitted.
     * examples:
     * If expecting a matrix, and receive shape (5), take it as a matrix of size (5,1)
     * If expecting a cube, and receive shape (5), take it as a cube of size (5,1,1)
     * If expecting a cube, and receive shape (2,3), take it as a cube of size (2,3,1)
     */
    template<class T>
    int arma_numpy_dimscheck( PyObject* input, bool strict=false){
        constexpr int dims = arma_info<T>::dims;

        PyArrayObject* array = (PyArrayObject*)input;
        int input_dims = array_numdims(array);

        if( strict ){
            return input_dims == dims;
        } else {
            return input_dims <= dims;
        }
    }

    /* Order check
     * Tests if the input array is Fortran ordered. Additionally checks that the input array:
     * - owns its own data
     * - is aligned in memory
     * - is writeable
     */
    int arma_numpy_ordercheck(PyObject* input){
        PyArrayObject* array = (PyArrayObject*)input;
        return PyArray_CHKFLAGS(array,(NPY_ARRAY_F_CONTIGUOUS|NPY_ARRAY_OWNDATA|NPY_ARRAY_ALIGNED|NPY_ARRAY_WRITEABLE)); 
    }

    /* Typecheck to ensure the input array has acceptable type and shape.
     * If a numpy array is not contiguous or is not column-ordered, a new array will be constructed for the conversion.
     * Numpy arrays may be cast to a different type, though only if it can be done so safely.
     */
    template<class T>
    int arma_numpy_typecheck(PyObject* input, bool strict=false){
        if( strict ){
            return (is_array(input) && arma_numpy_ordercheck(input) && arma_numpy_castcheck<T>(input,true) && arma_numpy_dimscheck<T>(input,true));
        } else {
            return (is_array(input) && arma_numpy_castcheck<T>(input) && arma_numpy_dimscheck<T>(input));
        }
    }
    
    /* Converts numpy array to Arma container.
     * Calls 'new' internally. Any Arma objects created this way must be deleted. 
     */
    template<class T>
    T* numpy_to_arma(PyObject* input){
        // Determine internal type of Vec and corresponding numpy typecode
        using element_t = typename arma_info<T>::element_t;
        static constexpr int typecode = arma_info<T>::typecode;
        // Convert generic PyObject to Numpy array.
        // If we need to convert it to a Fortran-contiguous copy or change the type, do so.
        // Note that this will not happen if an object passes strict typechecking.
        PyArrayObject* array = NULL;
        int is_new_obj;
        array = obj_to_array_fortran_allow_conversion( input, typecode, &is_new_obj);
        // Get dimensionality of numpy array and pointer to its raw data
        npy_intp* dims = array_dimensions(array);
        int ndims = array_numdims(array);
        element_t* data = reinterpret_cast<element_t*>(array_data(array));
        // Build arma object directly from Numpy
        return arma_from_ptr<T>(data,dims,ndims);
    }

    /* Converts arma container to numpy array
     * If 'copy', new memory is allocated and elements are copied over. This is necessary when returning by value.
     * Otherwise, a new numpy array is created that uses the arma array as its base.
     */
    template<class T>
    PyObject* arma_to_numpy( T* t, bool copy=true){
        // Get element type and corresponding type code
        using element_t = typename arma_info<T>::element_t;
        static constexpr int typecode = arma_info<T>::typecode;
        static constexpr int dims = arma_info<T>::dims;
        PyObject* array;
        // Get shape of arma container
        std::array<npy_intp,dims> shape = arma_shape(t);
        // Set base object to Arma array, or copy elements over.
        if(copy){
            // Create new empty array with same shape as Arma object
            array = PyArray_EMPTY( dims, shape.data(), typecode, /*'fortran' ordering*/ true);
            // Copy elements over, return
            std::copy( t->begin(), t->end(), reinterpret_cast<element_t*>(array_data(array)));
        } else {
            // Create PyArmaObject wrapper
            PyArmaObject<T>* wrap = PyObject_New( PyArmaObject<T>, &PyArmaTypeObject<T>::value);
            wrap->data = t;
            // Create new PyArray using Arma data
            array = PyArray_New(&PyArray_Type,dims,shape.data(),typecode,NULL,wrap->data->memptr(),0,NPY_ARRAY_F_CONTIGUOUS,NULL);
            // Set this wrapped object as the base of the numpy array
            PyArray_SetBaseObject( (PyArrayObject*) array, (PyObject*) wrap);
        }
        return array;
    }
}

// Finalise arma wrapper types

%init %{

#define FINALISE_ARMA_WRAPPER(_T)\
    PyType_Ready(&PyArmaTypeObject<_T<int>>::value);\
    PyType_Ready(&PyArmaTypeObject<_T<unsigned>>::value);\
    PyType_Ready(&PyArmaTypeObject<_T<long long>>::value);\
    PyType_Ready(&PyArmaTypeObject<_T<unsigned long long>>::value);\
    PyType_Ready(&PyArmaTypeObject<_T<float>>::value);\
    PyType_Ready(&PyArmaTypeObject<_T<double>>::value);\
    PyType_Ready(&PyArmaTypeObject<_T<std::complex<float>>>::value);\
    PyType_Ready(&PyArmaTypeObject<_T<std::complex<double>>>::value);

    FINALISE_ARMA_WRAPPER(arma::Col)
    FINALISE_ARMA_WRAPPER(arma::Row)
    FINALISE_ARMA_WRAPPER(arma::Mat)
    FINALISE_ARMA_WRAPPER(arma::Cube)
%}

// Create typemaps for armadillo typedefs
// (assumes ARMA_64BIT_WORD)
%apply long long { arma::sword};
%apply long long unsigned { arma::uword};
%apply std::complex<float> { arma::cx_float};
%apply std::complex<double> { arma::cx_double};

// Create macro for arma container typemaps
%define %gen_typemaps(T,prec)

    // Typecheck
    %typemap(typecheck, precedence=prec, fragment="arma_numpy") T, const T, T*, const T*, T&, const T& {
        $1 = arma_numpy_typecheck<T>($input);
    }

    // In by value
    // Always copies and quietly converts if it needs to.
    %typemap(in,  fragment="arma_numpy") T (T* temp), const T (T* temp) {
        temp = NULL;
        try{
            temp = numpy_to_arma<T>($input);
            $1 = *temp;
        } catch (const ArmaNumpyException& e){
            PyErr_SetString( PyExc_TypeError, e.what());
            return NULL;
        }
    }

    // Clean-up by value
    %typemap(freearg) T, const T {
        if(temp$argnum != NULL) delete temp$argnum;
    }

    // In by reference
    // Only copies if the type required by arma doesn't match the type provided by numpy, or if the input isn't Fortran ordered.
    // See argout for handling in case copying was performed.
    %typemap(in,  fragment="arma_numpy") T&, T*, const T&, const T* {
        try{
            $1 = numpy_to_arma<T>($input);
        } catch (const ArmaNumpyException& e){
            PyErr_SetString( PyExc_TypeError, e.what());
            return NULL;
        }
    }

    // Argout by reference:
    // The input may be converted before passing to arma.
    // To maintain the illusion of pass-by-reference, update input numpy object to reflect any changes.
    %typemap(argout, fragment="arma_numpy") T&, T* {
        // Check if the creation of an Arma object resulted in a conversion/copy.
        // - If no copy was made, and we may freely delete the Arma object.
        //   This shouldn't free the array data, as the Arma object was created by passing a pointer.
        // - If a copy was made when converting to arma, $input must be updated to reflect any changes made.
        bool strict_typecheck = arma_numpy_typecheck<T>($input,true);
        bool arma_reallocated = (array_data($input) == (char*) $1->memptr());
        if(!strict_typecheck || arma_reallocated ){
            // convert arma back to numpy
            PyObject* np = arma_to_numpy<T>($1,false);
            // move back into original array
            int err = PyArray_MoveInto((PyArrayObject*)$input,(PyArrayObject*)np);
            if(err==-1){
                PyErr_Format( PyExc_TypeError, "ArmaNumpyError: Conversion error when passing by reference in function %s", "$symname");
                return NULL;
            }
        } else {
            delete $1;
        }
    }

    // Argout by const reference:
    // Something passed by const ref/ptr can't be modified by the function, so there's no need to convert back again.
    %typemap(argout, fragment="arma_numpy") const T&, const T* {
        delete $1;
    }

    // Out by value
    %typemap(out, optimal="1", fragment="arma_numpy") T {
        T result = $1;
        $result = arma_to_numpy<T>(&result);
    }

    // Out by reference
    %typemap(out, fragment="arma_numpy") T*, T& {
        // Returns Numpy object with an Arma object as its base.
        // The dealloc method of the Arma wrapper calls 'delete', so it can only take Arma objects created by 'new'.
        // Therefore must create a new Arma object using the same memptr (this shouldn't involve any copying, besides nrows, ncols etc)
        static constexpr int dims = arma_info<T>::dims;
        std::array<npy_intp,dims> shape = arma_shape<T>($1);
        T* new_arma = arma_from_ptr<T>( (typename arma_info<T>::element_t*)($1)->memptr(), shape.data(), dims);
        $result = arma_to_numpy<T>(new_arma,false);
    }

%enddef

// Some preprocessor magic...
#define GET_NTH_ARG( _1, _2, _3, _4, _5, N, ...) N
#define APPLY_SYMBOLS_1(_pre,_post,_a) (_pre _a _post)
#define APPLY_SYMBOLS_2(_pre,_post,_a,...) (_pre _a _post) , APPLY_SYMBOLS_1(_pre,_post,__VA_ARGS__)
#define APPLY_SYMBOLS_3(_pre,_post,_a,...) (_pre _a _post) , APPLY_SYMBOLS_2(_pre,_post,__VA_ARGS__)
#define APPLY_SYMBOLS_4(_pre,_post,_a,...) (_pre _a _post) , APPLY_SYMBOLS_3(_pre,_post,__VA_ARGS__)
#define APPLY_SYMBOLS_5(_pre,_post,_a,...) (_pre _a _post) , APPLY_SYMBOLS_4(_pre,_post,__VA_ARGS__)
#define APPLY_SYMBOLS(_pre,_post,...) \
    GET_NTH_ARG(__VA_ARGS__,APPLY_SYMBOLS_5,APPLY_SYMBOLS_4,APPLY_SYMBOLS_3,APPLY_SYMBOLS_2,APPLY_SYMBOLS_1)(_pre,_post,##__VA_ARGS__)

// Create macro to apply typemaps to all variations on the type T (e.g. vec, colvec*, const dvec&, dcolvec&, const Col<double>)
%define %apply_typemaps(T,...)
    %apply T  { APPLY_SYMBOLS(,,__VA_ARGS__) };
    %apply T& { APPLY_SYMBOLS(,&,__VA_ARGS__) };
    %apply T* { APPLY_SYMBOLS(,*,__VA_ARGS__) };
    %apply const T  { APPLY_SYMBOLS(const,,__VA_ARGS__) };
    %apply const T& { APPLY_SYMBOLS(const,&,__VA_ARGS__) };
    %apply const T* { APPLY_SYMBOLS(const,*,__VA_ARGS__) };
%enddef

// Define some new typechecking precedence values
// Hopefully the swig devs won't invalidate this anytime soon!

%define TYPECHECK_UINT32_1ARRAY         1044 %enddef // SWIG_TYPECHECK_INT32_ARRAY-1
%define TYPECHECK_INT32_1ARRAY          1045 %enddef // SWIG_TYPECHECK_INT32_ARRAY
%define TYPECHECK_UINT32_2ARRAY         1046 %enddef // SWIG_TYPECHECK_INT32_ARRAY+1
%define TYPECHECK_INT32_2ARRAY          1047 %enddef // SWIG_TYPECHECK_INT32_ARRAY+2
%define TYPECHECK_UINT32_3ARRAY         1048 %enddef // SWIG_TYPECHECK_INT32_ARRAY+3
%define TYPECHECK_INT32_3ARRAY          1049 %enddef // SWIG_TYPECHECK_INT32_ARRAY+4

%define TYPECHECK_UINT64_1ARRAY         1054 %enddef // SWIG_TYPECHECK_INT64_ARRAY-1
%define TYPECHECK_INT64_1ARRAY          1055 %enddef // SWIG_TYPECHECK_INT64_ARRAY
%define TYPECHECK_UINT64_2ARRAY         1056 %enddef // SWIG_TYPECHECK_INT64_ARRAY+1
%define TYPECHECK_INT64_2ARRAY          1057 %enddef // SWIG_TYPECHECK_INT64_ARRAY+2
%define TYPECHECK_UINT64_3ARRAY         1058 %enddef // SWIG_TYPECHECK_INT64_ARRAY+3
%define TYPECHECK_INT64_3ARRAY          1059 %enddef // SWIG_TYPECHECK_INT64_ARRAY+4

%define TYPECHECK_FLOAT_1ARRAY          1080 %enddef // SWIG_TYPECHECK_FLOAT_ARRAY
%define TYPECHECK_FLOAT_2ARRAY          1081 %enddef // SWIG_TYPECHECK_FLOAT_ARRAY+1
%define TYPECHECK_FLOAT_3ARRAY          1082 %enddef // SWIG_TYPECHECK_FLOAT_ARRAY+2
%define TYPECHECK_DOUBLE_1ARRAY         1090 %enddef // SWIG_TYPECHECK_DOUBLE_ARRAY
%define TYPECHECK_DOUBLE_2ARRAY         1091 %enddef // SWIG_TYPECHECK_DOUBLE_ARRAY+1
%define TYPECHECK_DOUBLE_3ARRAY         1092 %enddef // SWIG_TYPECHECK_DOUBLE_ARRAY+2

%define TYPECHECK_COMPLEX_FLOAT_1ARRAY  1093 %enddef // SWIG_TYPECHECK_DOUBLE_ARRAY+3
%define TYPECHECK_COMPLEX_FLOAT_2ARRAY  1094 %enddef // SWIG_TYPECHECK_DOUBLE_ARRAY+4
%define TYPECHECK_COMPLEX_FLOAT_3ARRAY  1095 %enddef // SWIG_TYPECHECK_DOUBLE_ARRAY+5
%define TYPECHECK_COMPLEX_DOUBLE_1ARRAY 1096 %enddef // SWIG_TYPECHECK_DOUBLE_ARRAY+6
%define TYPECHECK_COMPLEX_DOUBLE_2ARRAY 1097 %enddef // SWIG_TYPECHECK_DOUBLE_ARRAY+7
%define TYPECHECK_COMPLEX_DOUBLE_3ARRAY 1098 %enddef // SWIG_TYPECHECK_DOUBLE_ARRAY+8

// Generate typemaps for arma::Col

%gen_typemaps(arma::Col<int>,TYPECHECK_INT32_1ARRAY);
%gen_typemaps(arma::Col<unsigned>,TYPECHECK_UINT32_1ARRAY);
%gen_typemaps(arma::ivec,TYPECHECK_INT64_1ARRAY);
%gen_typemaps(arma::uvec,TYPECHECK_UINT64_1ARRAY);
%gen_typemaps(arma::fvec,TYPECHECK_FLOAT_1ARRAY);
%gen_typemaps(arma::dvec,TYPECHECK_DOUBLE_1ARRAY);
%gen_typemaps(arma::cx_fvec,TYPECHECK_COMPLEX_FLOAT_1ARRAY);
%gen_typemaps(arma::cx_dvec,TYPECHECK_COMPLEX_DOUBLE_1ARRAY);

%apply_typemaps( arma::ivec, arma::icolvec, arma::Col<long long>, arma::Col<arma::sword>);
%apply_typemaps( arma::uvec, arma::ucolvec, arma::Col<unsigned long long>, arma::Col<arma::uword>);
%apply_typemaps( arma::fvec, arma::fcolvec, arma::Col<float>);
%apply_typemaps( arma::dvec, arma::dcolvec, arma::colvec, arma::vec, arma::Col<double>);
%apply_typemaps( arma::cx_fvec, arma::cx_fcolvec, arma::Col<std::complex<float>>, arma::Col<arma::cx_float>);
%apply_typemaps( arma::cx_dvec, arma::cx_dcolvec, arma::cx_vec, arma::cx_colvec, arma::Col<std::complex<double>>, arma::Col<arma::cx_double>);

// Generate typemaps for arma::Row

%gen_typemaps(arma::Row<int>,TYPECHECK_INT32_1ARRAY);
%gen_typemaps(arma::Row<unsigned>,TYPECHECK_UINT32_1ARRAY);
%gen_typemaps(arma::irowvec,TYPECHECK_INT64_1ARRAY);
%gen_typemaps(arma::urowvec,TYPECHECK_UINT64_1ARRAY);
%gen_typemaps(arma::frowvec,TYPECHECK_FLOAT_1ARRAY);
%gen_typemaps(arma::drowvec,TYPECHECK_DOUBLE_1ARRAY);
%gen_typemaps(arma::cx_frowvec,TYPECHECK_COMPLEX_FLOAT_1ARRAY);
%gen_typemaps(arma::cx_drowvec,TYPECHECK_COMPLEX_DOUBLE_1ARRAY);

%apply_typemaps( arma::irowvec, arma::Row<long long>, arma::Row<arma::sword>);
%apply_typemaps( arma::urowvec, arma::Row<unsigned long long>, arma::Row<arma::uword>);
%apply_typemaps( arma::frowvec, arma::Row<float>);
%apply_typemaps( arma::drowvec, arma::rowvec, arma::Row<double>);
%apply_typemaps( arma::cx_frowvec, arma::Row<std::complex<float>>, arma::Row<arma::cx_float>);
%apply_typemaps( arma::cx_drowvec, arma::cx_rowvec, arma::Row<std::complex<double>>, arma::Row<arma::cx_double>);

// Generate typemaps for arma::Mat

%gen_typemaps(arma::Mat<int>,TYPECHECK_INT32_2ARRAY);
%gen_typemaps(arma::Mat<unsigned>,TYPECHECK_UINT32_2ARRAY);
%gen_typemaps(arma::imat,TYPECHECK_INT64_2ARRAY);
%gen_typemaps(arma::umat,TYPECHECK_UINT64_2ARRAY);
%gen_typemaps(arma::fmat,TYPECHECK_FLOAT_2ARRAY);
%gen_typemaps(arma::dmat,TYPECHECK_DOUBLE_2ARRAY);
%gen_typemaps(arma::cx_fmat,TYPECHECK_COMPLEX_FLOAT_2ARRAY);
%gen_typemaps(arma::cx_dmat,TYPECHECK_COMPLEX_DOUBLE_2ARRAY);

%apply_typemaps( arma::imat, arma::Mat<long long>, arma::Mat<arma::sword>);
%apply_typemaps( arma::umat, arma::Mat<unsigned long long>, arma::Mat<arma::uword>);
%apply_typemaps( arma::fmat, arma::Mat<float>);
%apply_typemaps( arma::dmat, arma::mat, arma::Mat<double>);
%apply_typemaps( arma::cx_fmat, arma::Mat<std::complex<float>>, arma::Mat<arma::cx_float>);
%apply_typemaps( arma::cx_dmat, arma::cx_mat, arma::Mat<std::complex<double>>, arma::Mat<arma::cx_double>);

// Generate typemaps for arma::Cube

%gen_typemaps(arma::Cube<int>,TYPECHECK_INT32_3ARRAY);
%gen_typemaps(arma::Cube<unsigned>,TYPECHECK_UINT32_3ARRAY);
%gen_typemaps(arma::icube,TYPECHECK_INT64_3ARRAY);
%gen_typemaps(arma::ucube,TYPECHECK_UINT64_3ARRAY);
%gen_typemaps(arma::fcube,TYPECHECK_FLOAT_3ARRAY);
%gen_typemaps(arma::dcube,TYPECHECK_DOUBLE_3ARRAY);
%gen_typemaps(arma::cx_fcube,TYPECHECK_COMPLEX_FLOAT_3ARRAY);
%gen_typemaps(arma::cx_dcube,TYPECHECK_COMPLEX_DOUBLE_3ARRAY);

%apply_typemaps( arma::icube, arma::Cube<long long>, arma::Cube<arma::sword>);
%apply_typemaps( arma::ucube, arma::Cube<unsigned long long>, arma::Cube<arma::uword>);
%apply_typemaps( arma::fcube, arma::Cube<float>);
%apply_typemaps( arma::dcube, arma::cube, arma::Cube<double>);
%apply_typemaps( arma::cx_fcube, arma::Cube<std::complex<float>>, arma::Cube<arma::cx_float>);
%apply_typemaps( arma::cx_dcube, arma::cx_cube, arma::Cube<std::complex<double>>, arma::Cube<arma::cx_double>);