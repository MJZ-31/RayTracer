#include <cmath>

template <Numeric T, std::size_t D>
double length(const Vector<T, D>& vector) {

    double output = 0;
    for (unsigned int i = 0; i < D; ++i)
        output += pow(vector[i], 2);
    return sqrt(output);
}

template <Numeric T, std::size_t D>
Vector<T, D> add(const Vector<T, D>& a, const Vector<T, D>& b) {

    Vector<T, D> output;
    for (unsigned int i = 0; i < D; ++i)
        output[i] = a[i] + b[i];
    return output;
}

template <Numeric T, std::size_t D>
Vector<T, D> subtract(const Vector<T, D>& a, const Vector<T, D>& b) {

    Vector<T, D> output;
    for (unsigned int i = 0; i < D; ++i)
        output[i] = a[i] - b[i];
    return output;
}

template <Numeric T, std::size_t D>
Vector<T, D> scale(const Vector<T, D>& vector, double scalar) {

    Vector<T, D> output;
    for (unsigned int i = 0; i < D; ++i)
        output[i] = vector[i] * scalar;
    return output;
}

template <Numeric T, std::size_t D>
T dotProduct(const Vector<T, D>& a, const Vector<T, D>& b) {

    double output = static_cast<double>(0);
    for (unsigned int i = 0; i < D; ++i)
        output += a[i] * b[i];
    return output;
}

template <Numeric T>
Vector<T, 3> crossProduct(const Vector<T, 3>& a, const Vector<T, 3>& b) {

    Vector<T, 3> output;
    for (unsigned int i = 0; i < 3; ++i)
        output[i] = a[(i + 1) % 3] * b[(i + 2) % 3] - a[(i + 2) % 3] * b[(i + 1) % 3];
    return output;
}

template <Numeric T, std::size_t D_input, std::size_t D_output>
Vector<T, D_output> transform(const Vector<T, D_input>& vector, const Matrix<T, D_input, D_output>& matrix) {

    Vector<T, D_output> output;
    for (unsigned int i = 0; i < D_output; ++i)
        output[i] = dotProduct(vector, row(matrix, i));
    return output;
}

template <Numeric T, std::size_t W, std::size_t H>
double determinent(const Matrix<T, W, H>& matrix) {

    return 0;
}

template <Numeric T, std::size_t W, std::size_t H>
Vector<T, W> row(const Matrix<T, W, H>& matrix, unsigned int index) {

    Vector<T, W> output;
    for (unsigned int i = 0; i < W; ++i)
        output[i] = matrix[i][index];
    return output;
}

template <Numeric T, std::size_t W, std::size_t H>
Vector<T, H> column(const Matrix<T, W, H>& matrix, unsigned int index) {

    Vector<T, H> output;
    for (unsigned int i = 0; i < H; ++i)
        output[i] = matrix[index][i];
    return output;
}

template <Numeric T, std::size_t W, std::size_t H>
Matrix<T, H, W> transpose(const Matrix<T, W, H>& matrix) {

    Matrix<T, H, W> output;
    for (unsigned int i = 0; i < W; ++i) {

        for (unsigned int u = 0; u < H; ++u)
            output[u][i] = matrix[i][u];
    }
    return output;
}

template <Numeric T, std::size_t W, std::size_t H>
Matrix<T, W, H> reduce(const Matrix<T, W, H>& matrix) {

    Matrix<T, W, H> output;
    
    return output;
}

template <Numeric T, std::size_t W, std::size_t H>
Matrix<T, W, H> inverse(const Matrix<T, W, H>& matrix) {

    Matrix<T, W, H> output;
    Matrix<T, 2 * W, H> extended;

    return output;
}

template <Numeric T, std::size_t W, std::size_t H>
Matrix<T, W, H> add(const Matrix<T, W, H>& a, const Matrix<T, W, H>& b) {

    Matrix<T, W, H> output;
    for (unsigned int i = 0; i < W; ++i) {

        for (unsigned int u = 0; u < H; ++u)
            output[i][u] = a[i][u] + b[i][u];
    }
    return output;
}

template <Numeric T, std::size_t W, std::size_t H>
Matrix<T, W, H> subtract(const Matrix<T, W, H>& a, const Matrix<T, W, H>& b) {

    Matrix<T, W, H> output;
    for (unsigned int i = 0; i < W; ++i) {

        for (unsigned int u = 0; u < H; ++u)
            output[i][u] = a[i][u] - b[i][u];
    }
    return output;
}

template <Numeric T, std::size_t W_output, std::size_t H_output, std::size_t D>
Matrix<T, W_output, H_output> compose(const Matrix<T, D, H_output>& a, const Matrix<T, W_output, D>& b) {

    Matrix<T, W_output, H_output> output;
    for (unsigned int i = 0; i < W_output; ++i) {

        for (unsigned int u = 0; u < H_output; ++u)
            output[i][u] = dotProduct(column(a, i), row(b, u));
    }
    return output;
}