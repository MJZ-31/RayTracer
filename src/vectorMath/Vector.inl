#include <cmath>

///* VECTOR *///

// CONSTRUCTORS AND DESTRUCTOR //

template <typename T>
Vector<T>::Vector(unsigned int dimension) {

    m_dimension = dimension;
    if (m_dimension == 0)
        arr_data = NULL;
    else
        arr_data = new T[m_dimension];
    for (unsigned int i = 0; i < m_dimension; ++i)
        arr_data[i] = 0;    // Default to a zero vector.
    return;
}

template <typename T>
Vector<T>::Vector(unsigned int dimension, const T* values) {

    m_dimension = dimension;
    if (m_dimension == 0)
        arr_data = NULL;
    else
        arr_data = new T[m_dimension];
    for (unsigned int i = 0; i < m_dimension; ++i)
        arr_data[i] = values[i];
    return;
}

template <typename T>
Vector<T>::Vector(const Vector<T>& vector) {

    m_dimension = vector.getDimension();
    if (m_dimension == 0)
        arr_data = NULL;
    else
        arr_data = new T[m_dimension];
    for (unsigned int i = 0; i < m_dimension; ++i)
        arr_data[i] = vector[i];
    return;
}

template <typename T>
Vector<T>::~Vector() {

    if (arr_data != NULL)
        delete [] arr_data;
    return;
}

// ACCESSORS //

template <typename T>
unsigned int Vector<T>::getDimension() const {

    return m_dimension;
}

template <typename T>
double Vector<T>::getLength() const {

    /*
      v
    |[x]|
    |[y]| = sqrt(x^2 + y^2 + z^2)
    |[z]|
    */

    double output = 0;
    for (unsigned int i = 0; i < m_dimension; ++i)
        output += pow(arr_data[i], 2);
    return sqrt(output);
}

template <typename T>
Vector<T> Vector<T>::getNormalizedVector() const {

    /* 
    n = v / |v|

    Scaling a Vector {v} by the inverse of its length produces a normalized Vector {n} with length 1.
    */

    return scale(1.0 / getLength());
}

template <typename T>
Vector<T> Vector<T>::scale(double scalar) const {

    /*
     u           v
    [x]       [x * s]
    [y] * s = [y * s]
    [z]       [z * s]
    */

    Vector output(m_dimension);
    for (unsigned int i = 0; i < m_dimension; ++i)
        output[i] = arr_data[i] * scalar;
    return output;
}

template <typename T>
double Vector<T>::dotProduct(const Vector<T>& vector) const {

    /*
     u      v
    [x1]   [x2]
    [y1] * [y2] = x1 * x2 + y1 * y2 + z1 * z2
    [z1]   [z2]
    */

    if (m_dimension != vector.getDimension())
        throw dimension_mismatch("Cannot calculate a dot product between a " + std::to_string(m_dimension) + "-dimensional Vector and a " + std::to_string(vector.getDimension()) + "-dimensional Vector. The Vectors must be the same dimension.");

    double output = 0;
    for (unsigned int i = 0; i < m_dimension; ++i)
        output += arr_data[i] * vector[i];
    return output;
}

template <typename T>
Vector<T> Vector<T>::crossProduct(const Vector<T>& vector) const {

    /*
     u      v              c
    [x1]   [x2]   [y1 * z2 - z1 * y2]
    [y1] x [y2] = [z1 * x2 - x1 * z2]
    [z1]   [z2]   [x1 * y2 - y1 * x2]
    */

    // A cross product can only be calculated with two 3-dimensional Vectors.
    if (m_dimension != 3 || vector.getDimension() != 3)
        throw dimension_mismatch("Cannot calculate a cross product between a " + std::to_string(m_dimension) + "-dimensional Vector and a " + std::to_string(vector.getDimension()) + "-dimensional Vector. Both Vectors must be 3-dimensional.");

    Vector output(m_dimension);
    for (unsigned int i = 0; i < m_dimension; ++i)
        output[i] = arr_data[(i + 1) % 3] * vector[(i + 2) % 3] - arr_data[(i + 2) % 3] * vector[(i + 1) % 3];
    return output;
}

template <typename T>
Vector<T> Vector<T>::add(const Vector<T>& vector) const {

    /*
     u      v         w
    [x1]   [x2]   [x1 + x2]
    [y1] + [y2] = [y1 + y2]
    [z1]   [z2]   [z1 + z2]

    The addition is simply performed on the values of the Vectors {u, v} on each independent axis.
    */

    if (m_dimension != vector.getDimension())
        throw dimension_mismatch("Cannot add a " + std::to_string(m_dimension) + "-dimensional Vector and a " + std::to_string(vector.getDimension()) + "-dimensional Vector together. The Vectors must be the same dimension.");

    Vector output(m_dimension);
    for (unsigned int i = 0; i < m_dimension; ++i)
        output[i] = arr_data[i] + vector[i];
    return output;
}

template <typename T>
Vector<T> Vector<T>::subtract(const Vector<T>& vector) const {

    /*
     u      v         w
    [x1]   [x2]   [x1 - x2]
    [y1] - [y2] = [y1 - y2]
    [z1]   [z2]   [z1 - z2]

    The subtraction is simply performed on the values of the Vectors {u, v} on each independent axis.
    */

    if (m_dimension != vector.getDimension())
        throw dimension_mismatch("Cannot subtract a " + std::to_string(vector.getDimension()) + "-dimensional Vector from a " + std::to_string(m_dimension) + "-dimensional Vector. The Vectors must be the same dimension.");

    Vector output(m_dimension);
    for (unsigned int i = 0; i < m_dimension; ++i)
        output[i] = arr_data[i] - vector[i];
    return output;
}

template <typename T>
Vector<T> Vector<T>::applyMatrix(const Matrix<T>& matrix) const {

    /*
     u        A       v
    [x1]   [a b c]   [x2]   [a * x2 + b * y2 + c * z2]
    [x1] = [d e f] * [y2] = [d * x2 + e * y2 + f * z2]
    [x1]   [g h i]   [z2]   [g * x2 + h * y2 + i * z2]

    Each value in the resulting Vector {u} is a dot product between the given Vector {v} and a row in the tranformation
    Matrix {A} at the same index of said value.
    */

    if (m_dimension != matrix.getWidth())
        throw dimension_mismatch("Cannot apply a " + std::to_string(matrix.getWidth()) + " value wide Matrix to a " + std::to_string(m_dimension) + "-dimensional Vector. The Vector dimension and Matrix width must be the same.");

    Vector output(matrix.getHeight());
    for (unsigned int i = 0; i < m_dimension; ++i)
        output[i] = dotProduct(matrix.getRow(i));
    return output;
}