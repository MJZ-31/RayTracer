#include <cmath>

template <RealNumberTypes T, std::size_t D>
double length(const Vector<T, D>& vector) {

    double output = 0;
    for (unsigned int i = 0; i < D; ++i)
        output += pow(vector[i], 2);
    return sqrt(output);
}

template <RealNumberTypes T1, RealNumberTypes T2, std::size_t D>
Vector<std::common_type_t<T1, T2>, D> add(const Vector<T1, D>& a, const Vector<T2, D>& b) {

    Vector<std::common_type_t<T1, T2>, D> output;
    for (unsigned int i = 0; i < D; ++i)
        output[i] = a[i] + b[i];
    return output;
}

template <RealNumberTypes T1, RealNumberTypes T2, std::size_t D>
Vector<std::common_type_t<T1, T2>, D> subtract(const Vector<T1, D>& a, const Vector<T2, D>& b) {

    Vector<std::common_type_t<T1, T2>, D> output;
    for (unsigned int i = 0; i < D; ++i)
        output[i] = a[i] - b[i];
    return output;
}

template <RealNumberTypes T1, RealNumberTypes T2, std::size_t D>
Vector<std::common_type_t<T1, T2>, D> scale(const Vector<T1, D>& vector, T2 scalar) {

    Vector<std::common_type_t<T1, T2>, D> output;
    for (unsigned int i = 0; i < D; ++i)
        output[i] = vector[i] * scalar;
    return output;
}

template <RealNumberTypes T1, RealNumberTypes T2, std::size_t D>
std::common_type_t<T1, T2> dotProduct(const Vector<T1, D>& a, const Vector<T2, D>& b) {

    std::common_type_t<T1, T2> output = static_cast<std::common_type_t<T1, T2>>(0);
    for (unsigned int i = 0; i < D; ++i)
        output += a[i] * b[i];
    return output;
}

template <RealNumberTypes T1, RealNumberTypes T2>
Vector<std::common_type_t<T1, T2>, 3> crossProduct(const Vector<T1, 3>& a, const Vector<T2, 3>& b) {

    Vector<std::common_type_t<T1, T2>, 3> output;
    for (unsigned int i = 0; i < 3; ++i)
        output[i] = a[(i + 1) % 3] * b[(i + 2) % 3] - a[(i + 2) % 3] * b[(i + 1) % 3];
    return output;
}

template <RealNumberTypes T, std::size_t D_input, std::size_t D_output>
Vector<T, D_output> transform(const Vector<T, D_input>& vector, const Matrix<T, D_input, D_output>& matrix) {

    Vector<T, D_output> output;
    for (unsigned int i = 0; i < D_output; ++i)
        output[i] = dotProduct(vector, row(matrix, i));
    return output;
}

template <RealNumberTypes T, std::size_t W, std::size_t H>
T determinant(const Matrix<T, W, H>& matrix) {

    return 0;
}

template <RealNumberTypes T, std::size_t W, std::size_t H>
Vector<T, W> row(const Matrix<T, W, H>& matrix, unsigned int index) {

    Vector<T, W> output;
    for (unsigned int i = 0; i < W; ++i)
        output[i] = matrix[i][index];
    return output;
}

template <RealNumberTypes T, std::size_t W, std::size_t H>
Vector<T, H> column(const Matrix<T, W, H>& matrix, unsigned int index) {

    Vector<T, H> output;
    for (unsigned int i = 0; i < H; ++i)
        output[i] = matrix[index][i];
    return output;
}

template <RealNumberTypes T, std::size_t W, std::size_t H>
Matrix<T, H, W> transpose(const Matrix<T, W, H>& matrix) {

    Matrix<T, H, W> output;
    for (unsigned int i = 0; i < W; ++i) {

        for (unsigned int u = 0; u < H; ++u)
            output[u][i] = matrix[i][u];
    }
    return output;
}

template <RealNumberTypes T, std::size_t W, std::size_t H>
Matrix<T, W, H> reduce(const Matrix<T, W, H>& matrix) {

    Matrix<T, W, H> output;
    
    return output;
}

template <RealNumberTypes T, std::size_t W, std::size_t H>
Matrix<T, W, H> inverse(const Matrix<T, W, H>& matrix) {

    Matrix<T, W, H> output;
    Matrix<T, 2 * W, H> extended;

    return output;
}

template <RealNumberTypes T1, RealNumberTypes T2, std::size_t W, std::size_t H>
Matrix<std::common_type_t<T1, T2>, W, H> add(const Matrix<T1, W, H>& a, const Matrix<T2, W, H>& b) {

    Matrix<std::common_type_t<T1, T2>, W, H> output;
    for (unsigned int i = 0; i < W; ++i) {

        for (unsigned int u = 0; u < H; ++u)
            output[i][u] = a[i][u] + b[i][u];
    }
    return output;
}

template <RealNumberTypes T1, RealNumberTypes T2, std::size_t W, std::size_t H>
Matrix<std::common_type_t<T1, T2>, W, H> subtract(const Matrix<T1, W, H>& a, const Matrix<T2, W, H>& b) {

    Matrix<std::common_type_t<T1, T2>, W, H> output;
    for (unsigned int i = 0; i < W; ++i) {

        for (unsigned int u = 0; u < H; ++u)
            output[i][u] = a[i][u] - b[i][u];
    }
    return output;
}

template <RealNumberTypes T1, RealNumberTypes T2, std::size_t W, std::size_t H, std::size_t D>
Matrix<std::common_type_t<T1, T2>, W, H> compose(const Matrix<T1, D, H>& a, const Matrix<T2, W, D>& b) {

    Matrix<std::common_type_t<T1, T2>, W, H> output;
    for (unsigned int i = 0; i < W; ++i) {

        for (unsigned int u = 0; u < H; ++u)
            output[i][u] = dotProduct(row(a, i), column(b, u));
    }
    return output;
}