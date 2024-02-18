#include <cmath>

template <std::size_t D>
double length(const Vector<D>& vector) {

    double output = 0;
    for (unsigned int i = 0; i < D; ++i)
        output += pow(vector[i], 2);
    return sqrt(output);
}

// Adds together Vector a and Vector b.
template <std::size_t D>
Vector<D> add(const Vector<D>& a, const Vector<D>& b) {

    Vector<D> output;
    for (unsigned int i = 0; i < D; ++i)
        output[i] = a[i] + b[i];
    return output;
}

// Subtracts Vector b from Vector a.
template <std::size_t D>
Vector<D> subtract(const Vector<D>& a, const Vector<D>& b) {

    Vector<D> output;
    for (unsigned int i = 0; i < D; ++i)
        output[i] = a[i] - b[i];
    return output;
}

// Scales the length of the given Vector by the given scalar.
template <std::size_t D>
Vector<D> scale(const Vector<D>& vector, double scalar) {

    Vector<D> output;
    for (unsigned int i = 0; i < D; ++i)
        output[i] = vector[i] * scalar;
    return output;
}

template <std::size_t D>
double dotProduct(const Vector<D>& a, const Vector<D>& b) {

    double output = static_cast<double>(0);
    for (unsigned int i = 0; i < D; ++i)
        output += a[i] * b[i];
    return output;
}

Vector<3> crossProduct(const Vector<3>& a, const Vector<3>& b) {

    Vector<3> output;
    for (unsigned int i = 0; i < 3; ++i)
        output[i] = a[(i + 1) % 3] * b[(i + 2) % 3] - a[(i + 2) % 3] * b[(i + 1) % 3];
    return output;
}

template <std::size_t D_input, std::size_t D_output>
Vector<D_output> transform(const Vector<D_input>& vector, const Matrix<D_input, D_output>& matrix) {

    Vector<D_output> output;
    for (unsigned int i = 0; i < D_output; ++i)
        output[i] = dotProduct(vector, row(matrix, i));
    return output;
}

template <std::size_t W, std::size_t H>
double determinent(const Matrix<W, H>& matrix) {

    return 0;
}

template <std::size_t W, std::size_t H>
Vector<W> row(const Matrix<W, H>& matrix, unsigned int index) {

    Vector<W> output;
    for (unsigned int i = 0; i < W; ++i)
        output[i] = matrix[i][index];
    return output;
}

template <std::size_t W, std::size_t H>
Vector<H> column(const Matrix<W, H>& matrix, unsigned int index) {

    Vector<H> output;
    for (unsigned int i = 0; i < H; ++i)
        output[i] = matrix[index][i];
    return output;
}

template <std::size_t W, std::size_t H>
Matrix<H, W> transpose(const Matrix<W, H>& matrix) {

    Matrix<H, W> output;
    for (unsigned int i = 0; i < W; ++i) {

        for (unsigned int u = 0; u < H; ++u)
            output[u][i] = matrix[i][u];
    }
    return output;
}

template <std::size_t W, std::size_t H>
Matrix<W, H> reduce(const Matrix<W, H>& matrix) {

    Matrix<W, H> output;
    
    return output;
}

template <std::size_t W, std::size_t H>
Matrix<W, H> inverse(const Matrix<W, H>& matrix) {

    Matrix<W, H> output;
    Matrix<2 * W, H> extended;

    return output;
}

template <std::size_t W, std::size_t H>
Matrix<W, H> add(const Matrix<W, H>& a, const Matrix<W, H>& b) {

    Matrix<W, H> output;
    for (unsigned int i = 0; i < W; ++i) {

        for (unsigned int u = 0; u < H; ++u)
            output[i][u] = a[i][u] + b[i][u];
    }
    return output;
}

template <std::size_t W, std::size_t H>
Matrix<W, H> subtract(const Matrix<W, H>& a, const Matrix<W, H>& b) {

    Matrix<W, H> output;
    for (unsigned int i = 0; i < W; ++i) {

        for (unsigned int u = 0; u < H; ++u)
            output[i][u] = a[i][u] - b[i][u];
    }
    return output;
}

template <std::size_t W_output, std::size_t H_output, std::size_t D>
Matrix<W_output, H_output> compose(const Matrix<D, H_output>& a, const Matrix<W_output, D>& b) {

    Matrix<W_output, H_output> output;
    for (unsigned int i = 0; i < W_output; ++i) {

        for (unsigned int u = 0; u < H_output; ++u)
            output[i][u] = dotProduct(column(a, i), row(b, u));
    }
    return output;
}