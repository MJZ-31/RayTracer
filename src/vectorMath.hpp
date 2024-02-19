#include <array>
#include <concepts>
#include <initializer_list>

template <typename T>
concept RealNumberTypes = std::integral<T> || std::floating_point<T>;

template <RealNumberTypes T, std::size_t D>
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

template <RealNumberTypes T, std::size_t W, std::size_t H>
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
template <RealNumberTypes T, std::size_t D>
double length(const Vector<T, D>& vector);

// Adds together Vector a and Vector b.
template <RealNumberTypes T1, RealNumberTypes T2, std::size_t D>
Vector<std::common_type_t<T1, T2>, D> add(const Vector<T1, D>& a, const Vector<T2, D>& b);

template <RealNumberTypes T1, RealNumberTypes T2, std::size_t D>
Vector<std::common_type_t<T1, T2>, D> operator+(const Vector<T1, D>& a, const Vector<T2, D>& b) { return add(a, b); };

// Subtracts Vector b from Vector a.
template <RealNumberTypes T1, RealNumberTypes T2, std::size_t D>
Vector<std::common_type_t<T1, T2>, D> subtract(const Vector<T1, D>& a, const Vector<T2, D>& b);

template <RealNumberTypes T1, RealNumberTypes T2, std::size_t D>
Vector<std::common_type_t<T1, T2>, D> operator-(const Vector<T1, D>& a, const Vector<T2, D>& b) { return subtract(a, b); }

// Scales the length of the given Vector by the given scalar.
template <RealNumberTypes T1, RealNumberTypes T2, std::size_t D>
Vector<std::common_type_t<T1, T2>, D> scale(const Vector<T1, D>& vector, T2 scalar);

template <RealNumberTypes T1, RealNumberTypes T2, std::size_t D>
Vector<std::common_type_t<T1, T2>, D> operator*(const Vector<T1, D>& vector, T2 scalar) { return scale(vector, scalar); }
template <RealNumberTypes T1, RealNumberTypes T2, std::size_t D>
Vector<std::common_type_t<T1, T2>, D> operator*(T1 scalar, const Vector<T2, D>& vector) { return scale(vector, scalar); }

// Calculates the dot product between Vector a and Vector b.
template <RealNumberTypes T1, RealNumberTypes T2, std::size_t D>
std::common_type_t<T1, T2> dotProduct(const Vector<T1, D>& a, const Vector<T2, D>& b);

template <RealNumberTypes T1, RealNumberTypes T2, std::size_t D>
std::common_type_t<T1, T2> operator*(const Vector<T1, D>& a, const Vector<T2, D>& b) { return dotProduct(a, b); };

// Calculates the cross product between Vector a and Vector b.
template <RealNumberTypes T1, RealNumberTypes T2>
Vector<std::common_type_t<T1, T2>, 3> crossProduct(const Vector<T1, 3>& a, const Vector<T2, 3>& b);

template <RealNumberTypes T1, RealNumberTypes T2>
Vector<std::common_type_t<T1, T2>, 3> operator%(const Vector<T1, 3>& a, const Vector<T2, 3>& b) { return crossProduct(a, b); };

// Calculates the transformed Vector by applying the given Matrix to the given Vector.
template <RealNumberTypes T1, RealNumberTypes T2, std::size_t D1, std::size_t D2>
Vector<std::common_type_t<T1, T2>, D2> transform(const Vector<T1, D1>& vector, const Matrix<T2, D1, D2>& matrix);

template <RealNumberTypes T1, RealNumberTypes T2, std::size_t D1, std::size_t D2>
Vector<std::common_type_t<T1, T2>, D2> operator*(const Vector<T1, D1>& vector, const Matrix<T2, D1, D2>& matrix) { return transform(vector, matrix); };
template <RealNumberTypes T1, RealNumberTypes T2, std::size_t D1, std::size_t D2>
Vector<std::common_type_t<T1, T2>, D2> operator*(const Matrix<T1, D1, D2>& matrix, const Vector<T2, D1>& vector) { return transform(vector, matrix); };

// Calculates the determinant of the given Matrix.
template <RealNumberTypes T, std::size_t W, std::size_t H>
T determinant(const Matrix<T, W, H>& matrix);

// Retreives a row from the given Matrix.
template <RealNumberTypes T, std::size_t W, std::size_t H>
Vector<T, W> row(const Matrix<T, W, H>& matrix, unsigned int index);

// Retreives a column from the given Matrix.
template <RealNumberTypes T, std::size_t W, std::size_t H>
Vector<T, H> column(const Matrix<T, W, H>& matrix, unsigned int index);

// Calculates the transpose Matrix of the given Matrix.
template <RealNumberTypes T, std::size_t W, std::size_t H>
Matrix<T, H, W> transpose(const Matrix<T, W, H>& matrix);

// Calculates the row reduced Matrix of the given Matrix.
template <RealNumberTypes T, std::size_t W, std::size_t H>
Matrix<T, W, H> reduce(const Matrix<T, W, H>& matrix);

// Calculates the inverse Matrix of the given Matrix.
template <RealNumberTypes T, std::size_t W, std::size_t H>
Matrix<T, W, H> inverse(const Matrix<T, W, H>& matrix);

// Adds together Matrix a and Matrix b.
template <RealNumberTypes T1, RealNumberTypes T2, std::size_t W, std::size_t H>
Matrix<std::common_type_t<T1, T2>, W, H> add(const Matrix<T1, W, H>& a, const Matrix<T2, W, H>& b);

template <RealNumberTypes T1, RealNumberTypes T2, std::size_t W, std::size_t H>
Matrix<std::common_type_t<T1, T2>, W, H> operator+(const Matrix<T1, W, H>& a, const Matrix<T2, W, H>& b) { return add(a, b); }

// Subtracts Matrix b from Matrix a.
template <RealNumberTypes T1, RealNumberTypes T2, std::size_t W, std::size_t H>
Matrix<std::common_type_t<T1, T2>, W, H> subtract(const Matrix<T1, W, H>& a, const Matrix<T2, W, H>& b);

template <RealNumberTypes T1, RealNumberTypes T2, std::size_t W, std::size_t H>
Matrix<std::common_type_t<T1, T2>, W, H> operator-(const Matrix<T1, W, H>& a, const Matrix<T2, W, H>& b) { return subtract(a, b); }

// Composes Matrix a and Matrix b together.
template <RealNumberTypes T1, RealNumberTypes T2, std::size_t W, std::size_t H, std::size_t D>
Matrix<std::common_type_t<T1, T2>, W, H> compose(const Matrix<T1, D, H>& a, const Matrix<T2, W, D>& b);

template <RealNumberTypes T1, RealNumberTypes T2, std::size_t W, std::size_t H, std::size_t D>
Matrix<std::common_type_t<T1, T2>, W, H> operator*(const Matrix<T1, D, H>& a, const Matrix<T2, W, D>& b) { return compose(a, b); }

#include "vectorMath.inl"