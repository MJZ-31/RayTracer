#include <array>
#include <concepts>
#include <initializer_list>

template<typename T>
concept Numeric = std::signed_integral<T> || std::floating_point<T>;

template <Numeric T, std::size_t D>
class Vector : public std::array<T, D> {

    public:

        Vector() {

            for (unsigned int i = 0; i < D; ++i)
                (*this)[i] = static_cast<T>(0);
        }
        Vector(const std::initializer_list<T>& list) {

            for (unsigned int i = 0; i < D; ++i)
                (*this)[i] = *(data(list) + i);
        }
};

template <Numeric T, std::size_t W, std::size_t H>
class Matrix : public std::array<Vector<T, H>, W> {

    public:

        Matrix() {

            for (unsigned int i = 0; i < W; ++i) {

                for (unsigned int u = 0; u < H; ++u)
                    (*this)[i][u] = static_cast<T>(0);
            }
        }
        Matrix(const std::initializer_list<std::initializer_list<T>>& list) {

            for (unsigned int i = 0; i < W; ++i) {

                for (unsigned int u = 0; u < H; ++u)
                    (*this)[i][u] = *(data(*(data(list) + i)) + u);
            }
        }
};

// Calculates the length of the given Vector.
template <Numeric T, std::size_t D>
double length(const Vector<T, D>& vector);

// Adds together Vector a and Vector b.
template <Numeric T, std::size_t D>
Vector<T, D> add(const Vector<T, D>& a, const Vector<T, D>& b);

template <Numeric T, std::size_t D>
double operator+(const Vector<T, D>& a, const Vector<T, D>& b) { return add(a, b); }

// Subtracts Vector b from Vector a.
template <Numeric T, std::size_t D>
Vector<T, D> subtract(const Vector<T, D>& a, const Vector<T, D>& b);

template <Numeric T, std::size_t D>
Vector<T, D> operator-(const Vector<T, D>& a, const Vector<T, D>& b) { return subtract(a, b); }

// Scales the length of the given Vector by the given scalar.
template <Numeric T, std::size_t D>
Vector<T, D> scale(const Vector<T, D>& vector, double scalar);

template <Numeric T, std::size_t D>
double operator*(const Vector<T, D>& vector, double scalar) { return scale(vector, scalar); }
template <Numeric T, std::size_t D>
double operator*(double scalar, const Vector<T, D>& vector) { return scale(vector, scalar); }

// Calculates the dot product between Vector a and Vector b.
template <Numeric T, std::size_t D>
T dotProduct(const Vector<T, D>& a, const Vector<T, D>& b);

template <Numeric T, std::size_t D>
T operator*(const Vector<T, D>& a, const Vector<T, D>& b) { return dotProduct(a, b); };

// Calculates the cross product between Vector a and Vector b.
template <Numeric T>
Vector<T, 3> crossProduct(const Vector<T, 3>& a, const Vector<T, 3>& b);

template <Numeric T>
Vector<T, 3> operator%(const Vector<T, 3>& a, const Vector<T, 3>& b) { return crossProduct(a, b); };

// Calculates the transformed Vector by applying the given Matrix to the given Vector.
template <Numeric T, std::size_t D_input, std::size_t D_output>
Vector<T, D_output> transform(const Vector<T, D_input>& vector, const Matrix<T, D_input, D_output>& matrix);

template <Numeric T, std::size_t D_input, std::size_t D_output>
Vector<T, D_output> operator*(const Vector<T, D_input>& vector, const Matrix<T, D_input, D_output>& matrix) { return transform(vector, matrix); };
template <Numeric T, std::size_t D_input, std::size_t D_output>
Vector<T, D_output> operator*(const Matrix<T, D_input, D_output>& matrix, const Vector<T, D_input>& vector) { return transform(vector, matrix); };

// Calculates the determinant of the given Matrix.
template <Numeric T, std::size_t W, std::size_t H>
double determinent(const Matrix<T, W, H>& matrix);

// Retreives a row from the given Matrix.
template <Numeric T, std::size_t W, std::size_t H>
Vector<T, W> row(const Matrix<T, W, H>& matrix, unsigned int index);

// Retreives a column from the given Matrix.
template <Numeric T, std::size_t W, std::size_t H>
Vector<T, H> column(const Matrix<T, W, H>& matrix, unsigned int index);

// Calculates the transpose Matrix of the given Matrix.
template <Numeric T, std::size_t W, std::size_t H>
Matrix<T, H, W> transpose(const Matrix<T, W, H>& matrix);

// Calculates the row reduced Matrix of the given Matrix.
template <Numeric T, std::size_t W, std::size_t H>
Matrix<T, W, H> reduce(const Matrix<T, W, H>& matrix);

// Calculates the inverse Matrix of the given Matrix.
template <Numeric T, std::size_t W, std::size_t H>
Matrix<T, W, H> inverse(const Matrix<T, W, H>& matrix);

// Adds together Matrix a and Matrix b.
template <Numeric T, std::size_t W, std::size_t H>
Matrix<T, W, H> add(const Matrix<T, W, H>& a, const Matrix<T, W, H>& b);

template <Numeric T, std::size_t W, std::size_t H>
Matrix<T, W, H> operator+(const Matrix<T, W, H>& a, const Matrix<T, W, H>& b) { return add(a, b); }

// Subtracts Matrix b from Matrix a.
template <Numeric T, std::size_t W, std::size_t H>
Matrix<T, W, H> subtract(const Matrix<T, W, H>& a, const Matrix<T, W, H>& b);

template <Numeric T, std::size_t W, std::size_t H>
Matrix<T, W, H> operator-(const Matrix<T, W, H>& a, const Matrix<T, W, H>& b) { return subtract(a, b); }

// Composes Matrix a and Matrix b together.
template <Numeric T, std::size_t W, std::size_t H, std::size_t D>
Matrix<T, W, H> compose(const Matrix<T, D, H>& a, const Matrix<T, W, D>& b);

template <Numeric T, std::size_t W, std::size_t H, std::size_t D>
Matrix<T, W, H> operator*(const Matrix<T, D, H>& a, const Matrix<T, W, D>& b) { return compose(a, b); }

#include "vectorMath.inl"