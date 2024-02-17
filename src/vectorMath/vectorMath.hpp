#include <array>
#include <cmath>
#include <exception>
#include <stdexcept>
#include <string>

class dimension_mismatch : public std::logic_error {

    public:

        dimension_mismatch(std::string msg) : std::logic_error(msg) {}
};

template <typename T>
class Vector;

template <typename T>
class Matrix;

template <typename T>
class Vector {

    private:

        unsigned int m_dimension;
        T* arr_data;

    public:

        // CONSTRUCTORS AND DESTRUCTOR //

        // Constructs a Vector with the given dimension and all values set to 0.
        Vector(unsigned int dimension);

        // Constructs a Vector with the given dimension and values.
        // PRE: values must be an array defined as T[dimension] - May result in memory errors if false.
        Vector(unsigned int dimension, const T* values);

        // Constructs a copy of the given Vector with a deep copy of arr_data.
        Vector(const Vector<T>& vector);

        // Detroys the calling Vector by deallocating arr_data.
        ~Vector();

        // ACCESSORS //

        // Returns the m_dimension of the calling Vector.
        unsigned int getDimension() const;

        // Calculates and returns the length of the calling Vector.
        double getLength() const;

        // Calculates and returns the normalized Vector of the calling Vector.
        Vector<T> getNormalizedVector() const;

        // Calculates and returns a Vector which results from scaling the calling Vector by the given scalar.
        Vector<T> scale(double scalar) const;

        // Calculates a vector dot product between the calling Vector and the argument Vector and returns the resulting double value.
        // PRE: this->m_dimension == vector.m_dimension - Throws dimension_mismatch if false.
        double dotProduct(const Vector<T>& vector) const;

        // Calculates a vector cross product between the calling Vector and the argument Vector and returns the resulting Vector.
        // PRE: this->m_dimension == 3 && vector.m_dimension == 3 - Throws dimension_mismatch if false.
        Vector<T> crossProduct(const Vector<T>& vector) const;

        // Adds together the calling Vector and the argument Vector and returns the resulting Vector.
        // PRE: this->m_dimension == vector.m_dimension - Throws dimension_mismatch if false.
        Vector<T> add(const Vector<T>& vector) const;

        // Subtracts the argument Vector from the calling Vector and returns the resulting Vector.
        // PRE: this->m_dimension == vector.m_dimension - Throws dimension_mismatch if false.
        Vector<T> subtract(const Vector<T>& vector) const;

        // Applies the argument Matrix to the calling Vector and returns the resulting Vector.
        // PRE: this->m_dimension = matrix.m_width - Throws dimension_mismatch if false.
        Vector<T> applyMatrix(const Matrix<T>& matrix) const;
        
        // MEMBER OPERATOR OVERLOADS //

        // Assignment operator. Reconstructs the calling Vector to have the same m_dimension as the argument Vector and a deep copy of the m_data in the argument Vector.
        Vector<T>& operator=(const Vector<T>& vector) {

            // Protect against self assignment.
            if (this == &vector)
                return *this;

            if (arr_data != NULL)
                delete [] arr_data;

            m_dimension = vector.getDimension();
            if (m_dimension == 0)
                arr_data = NULL;
            else
                arr_data = new T[m_dimension];
            for (unsigned int i = 0; i < m_dimension; ++i)
                arr_data[i] = vector[i];

            return *this;
        }

        // Array subscript operator. Returns a reference to a value in arr_data.
        // PRE: index < this->dimension - Throws std::out_of_range if false.
        T& operator[](unsigned int index) const {
            
            // Protect against out of range indecies.
            if (index >= m_dimension)
                throw std::out_of_range("An index of " + std::to_string(index) + " is out of range of a " + std::to_string(m_dimension) + "-dimensional Vector.");
            return arr_data[index];
        }
};

template <typename T>
class Matrix {

    private:

        unsigned int m_width;
        unsigned int m_height;
        Vector<T>* arr_data;       // Each Vector is a column in the Matrix.

    public:

        // CONSTRUCTORS AND DESTRUCTORS //

        // Constructs a Matrix with the given width and height and all values in arr_data set to 0.
        Matrix(unsigned int width, unsigned int height);

        // Constructs a Matrix with the given width, height, and values.
        // PRE: values must be an 2-dimensional array defined as T[width][height] - May result in memory errors if false.
        Matrix(unsigned int width, unsigned int height, const T** values);

        // Constructs a copy of the given Vector with a deep copy of arr_data.
        Matrix(const Matrix<T>& matrix);

        // Detroys the calling Matrix by deallocating arr_data.
        ~Matrix();

        // ACCESSORS //

        // Returns the m_width of the calling Matrix.
        unsigned int getWidth() const;

        // Returns the m_height of the calling Matrix.
        unsigned int getHeight() const;

        // Calculates the transpose of the calling Matrix and returns the resulting Matrix.
        Matrix<T> getTransposedMatrix() const;

        Matrix<T> getReducedMatrix() const;

        Matrix<T> getInverseMatrix() const;

        double getDeterminant() const;

        Vector<T> getRow(unsigned int index) const;

        Vector<T> getColumn(unsigned int index) const;

        Matrix<T> scale(unsigned int scalar) const;

        Vector<T> applyToVector(const Vector<T>& vector) const;

        Matrix<T> compose(const Matrix<T>& matrix) const;

        Matrix<T> add(const Matrix<T>& matrix) const;

        Matrix<T> subtract(const Matrix<T>& matrix) const;

        // MEMBER OPERATOR OVERLOADS //

        Matrix<T>& operator=(const Matrix<T>& matrix) {

            return *this;
        }

        Vector<T>& operator[](unsigned int index) {

            // Protect against out of range indecies.
            if (index >= m_width)
                throw std::out_of_range("An index of " + std::to_string(index) + " is out of range of a " + std::to_string(m_width) + " value wide Matrix.");
            return arr_data[index];
        }
};
#include "Vector.inl"
#include "Matrix.inl"
