#include <array>

template <std::size_t D>
using Vector = std::array<double, D>;

template <std::size_t W, std::size_t H>
using Matrix = std::array<std::array<double, H>, W>;

// Calculates the length of the given Vector.
template <std::size_t D>
double length(const Vector<D>& vector);

// Adds together Vector a and Vector b.
template <std::size_t D>
Vector<D> add(const Vector<D>& a, const Vector<D>& b);

template <std::size_t D>
double operator+(const Vector<D>& a, const Vector<D>& b) { return add(a, b); }

// Subtracts Vector b from Vector a.
template <std::size_t D>
Vector<D> subtract(const Vector<D>& a, const Vector<D>& b);

template <std::size_t D>
Vector<D> operator-(const Vector<D>& a, const Vector<D>& b) { return subtract(a, b); }

// Scales the length of the given Vector by the given scalar.
template <std::size_t D>
Vector<D> scale(const Vector<D>& vector, double scalar);

template <std::size_t D>
double operator*(const Vector<D>& vector, double scalar) { return scale(vector, scalar); }
template <std::size_t D>
double operator*(double scalar, const Vector<D>& vector) { return scale(vector, scalar); }

// Calculates the dot product between Vector a and Vector b.
template <std::size_t D>
double dotProduct(const Vector<D>& a, const Vector<D>& b);

template <std::size_t D>
double operator*(const Vector<D>& a, const Vector<D>& b) { return dotProduct(a, b); };

// Calculates the cross product between Vector a and Vector b.
Vector<3> crossProduct(const Vector<3>& a, const Vector<3>& b);

Vector<3> operator%(const Vector<3>& a, const Vector<3>& b) { return crossProduct(a, b); };

// Calculates the transformed Vector by applying the given Matrix to the given Vector.
template <std::size_t D_input, std::size_t D_output>
Vector<D_output> transform(const Vector<D_input>& vector, const Matrix<D_input, D_output>& matrix);

template <std::size_t D_input, std::size_t D_output>
Vector<D_output> operator*(const Vector<D_input>& vector, const Matrix<D_input, D_output>& matrix) { return transform(vector, matrix); };
template <std::size_t D_input, std::size_t D_output>
Vector<D_output> operator*(const Matrix<D_input, D_output>& matrix, const Vector<D_input>& vector) { return transform(vector, matrix); };

// Calculates the determinant of the given Matrix.
template <std::size_t W, std::size_t H>
double determinent(const Matrix<W, H>& matrix);


template <std::size_t W, std::size_t H>
Vector<W> row(const Matrix<W, H>& matrix, unsigned int index);

template <std::size_t W, std::size_t H>
Vector<H> column(const Matrix<W, H>& matrix, unsigned int index);

// Calculates the transpose Matrix of the given Matrix.
template <std::size_t W, std::size_t H>
Matrix<H, W> transpose(const Matrix<W, H>& matrix);

// Calculates the row reduced Matrix of the given Matrix.
template <std::size_t W, std::size_t H>
Matrix<W, H> reduce(const Matrix<W, H>& matrix);

// Calculates the inverse Matrix of the given Matrix.
template <std::size_t W, std::size_t H>
Matrix<W, H> inverse(const Matrix<W, H>& matrix);

// Adds together Matrix a and Matrix b.
template <std::size_t W, std::size_t H>
Matrix<W, H> add(const Matrix<W, H>& a, const Matrix<W, H>& b);

template <std::size_t W, std::size_t H>
Matrix<W, H> operator+(const Matrix<W, H>& a, const Matrix<W, H>& b) { return add(a, b); }

// Subtracts Matrix b from Matrix a.
template <std::size_t W, std::size_t H>
Matrix<W, H> subtract(const Matrix<W, H>& a, const Matrix<W, H>& b);

template <std::size_t W, std::size_t H>
Matrix<W, H> operator-(const Matrix<W, H>& a, const Matrix<W, H>& b) { return subtract(a, b); }

// Composes Matrix a and Matrix b together.
template <std::size_t W, std::size_t H, std::size_t D>
Matrix<W, H> compose(const Matrix<D, H>& a, const Matrix<W, D>& b);

template <std::size_t W, std::size_t H, std::size_t D>
Matrix<W, H> operator*(const Matrix<D, H>& a, const Matrix<W, D>& b) { return compose(a, b); }

#include "vectorMath.inl"