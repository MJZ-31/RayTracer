// Vector Member Function Definitions //

template <ScalarTypes T, std::size_t D>
Vector<T, D>::Vector() {

    for (std::size_t i_dimension = 0; i_dimension < D; ++i_dimension)
        (*this)[i_dimension] = static_cast<T>(0);
}

template <ScalarTypes T, std::size_t D> template <std::convertible_to<T> U>
Vector<T, D>::Vector(const std::initializer_list<U>& list) : Vector() {

    if (list.size() == 0)
        return;
    else if (list.size() == 1) {

        for (std::size_t i_dimension = 0; i_dimension < D; ++i_dimension)
            (*this)[i_dimension] = static_cast<T>(*list.begin());
    }
    else {

        for (std::size_t i_dimension = 0; i_dimension < D && i_dimension < list.size(); ++i_dimension)
            (*this)[i_dimension] = static_cast<T>(*(list.begin() + i_dimension));
    }
}

template <ScalarTypes T, std::size_t D> template <std::convertible_to<T> U>
Vector<T, D>::Vector(const Vector<U, D>& vector) {

    for (std::size_t i_dimension = 0; i_dimension < D; ++i_dimension)
        (*this)[i_dimension] = static_cast<T>(vector[i_dimension]);
}

// Matrix Member Function Definitions

template <ScalarTypes T, std::size_t W, std::size_t H>
Matrix<T, W, H>::Matrix() {

    for (std::size_t i_columns = 0; i_columns < W; ++i_columns)
        (*this)[i_columns] = Vector<T, H>();
}

template <ScalarTypes T, std::size_t W, std::size_t H> template <std::convertible_to<T> U>
Matrix<T, W, H>::Matrix(const std::initializer_list<std::initializer_list<U>>& list) : Matrix() {

    if (list.size() == 0)
        return;
    else if (list.size() == 1) {

        for (std::size_t i_columns = 0; i_columns < W; ++i_columns)
            (*this)[i_columns] = Vector<T, H>(*list.begin());
    }
    else {

        for (std::size_t i_columns = 0; i_columns < W && i_columns < list.size(); ++i_columns)
            (*this)[i_columns] = Vector<T, H>(*(list.begin() + i_columns));
    }
}

template <ScalarTypes T, std::size_t W, std::size_t H> template <std::convertible_to<T> U>
Matrix<T, W, H>::Matrix(const Matrix<U, W, H>& matrix) {

    for (std::size_t i_columns = 0; i_columns < W; ++i_columns)
        (*this)[i_columns] = matrix[i_columns];
}

// Calculation Function Definitions //

template <ScalarTypes T, std::size_t D>
double length(const Vector<T, D>& vector) {

    double output = 0;
    for (std::size_t i_dimension = 0; i_dimension < D; ++i_dimension)
        output += pow(vector[i_dimension], 2);
    return sqrt(output);
}

template <std::convertible_to<double> T, std::size_t D>
Vector<double, D> normalize(const Vector<T, D>& vector) {

    return scale(vector, 1 / length(vector));
}

template <ScalarTypes T1, ScalarTypes T2, std::size_t D>
Vector<std::common_type_t<T1, T2>, D> add(const Vector<T1, D>& a, const Vector<T2, D>& b) {

    Vector<std::common_type_t<T1, T2>, D> output;
    for (std::size_t i_dimension = 0; i_dimension < D; ++i_dimension)
        output[i_dimension] = a[i_dimension] + b[i_dimension];
    return output;
}

template <ScalarTypes T1, ScalarTypes T2, std::size_t D>
Vector<std::common_type_t<T1, T2>, D> subtract(const Vector<T1, D>& a, const Vector<T2, D>& b) {

    Vector<std::common_type_t<T1, T2>, D> output;
    for (std::size_t i_dimension = 0; i_dimension < D; ++i_dimension)
        output[i_dimension] = a[i_dimension] - b[i_dimension];
    return output;
}

template <ScalarTypes T1, ScalarTypes T2, std::size_t D>
Vector<std::common_type_t<T1, T2>, D> scale(const Vector<T1, D>& vector, T2 factor) {

    Vector<std::common_type_t<T1, T2>, D> output;
    for (std::size_t i_dimension = 0; i_dimension < D; ++i_dimension)
        output[i_dimension] = vector[i_dimension] * factor;
    return output;
}

template <ScalarTypes T1, ScalarTypes T2, std::size_t D>
std::common_type_t<T1, T2> dotProduct(const Vector<T1, D>& a, const Vector<T2, D>& b) {

    std::common_type_t<T1, T2> output = static_cast<std::common_type_t<T1, T2>>(0);
    for (std::size_t i_dimension = 0; i_dimension < D; ++i_dimension)
        output += a[i_dimension] * b[i_dimension];
    return output;
}

template <ScalarTypes T1, ScalarTypes T2>
Vector<std::common_type_t<T1, T2>, 3> crossProduct(const Vector<T1, 3>& a, const Vector<T2, 3>& b) {

    Vector<std::common_type_t<T1, T2>, 3> output;
    for (std::size_t i_dimension = 0; i_dimension < 3; ++i_dimension)
        output[i_dimension] = a[(i_dimension + 1) % 3] * b[(i_dimension + 2) % 3] - a[(i_dimension + 2) % 3] * b[(i_dimension + 1) % 3];
    return output;
}

template <ScalarTypes T1, ScalarTypes T2, std::size_t D1, std::size_t D2>
Vector<std::common_type<T1, T2>, D2> transform(const Vector<T1, D1>& vector, const Matrix<T2, D1, D2>& matrix) {

    Vector<std::common_type<T1, T2>, D2> output;
    for (std::size_t i_dimension = 0; i_dimension < D2; ++i_dimension)
        output[i_dimension] = dotProduct(vector, row(matrix, i_dimension));
    return output;
}

template <ScalarTypes T, std::size_t W, std::size_t H>
Vector<T, W> row(const Matrix<T, W, H>& matrix, std::size_t index) {

    Vector<T, W> output;
    for (std::size_t i_columns = 0; i_columns < W; ++i_columns)
        output[i_columns] = matrix[i_columns][index];
    return output;
}

template <ScalarTypes T, std::size_t W, std::size_t H>
Vector<T, H> column(const Matrix<T, W, H>& matrix, std::size_t index) {

    Vector<T, H> output;
    for (std::size_t i_rows = 0; i_rows < H; ++i_rows)
        output[i_rows] = matrix[index][i_rows];
    return output;
}

template <ScalarTypes T, std::size_t W, std::size_t H>
Matrix<T, H, W> transpose(const Matrix<T, W, H>& matrix) {

    Matrix<T, H, W> output;
    for (std::size_t i_columns = 0; i_columns < W; ++i_columns) {

        for (std::size_t i_rows = 0; i_rows < H; ++i_rows)
            output[i_rows][i_columns] = matrix[i_columns][i_rows];
    }
    return output;
}

// template <ScalarTypes T, std::size_t W, std::size_t H>
// Matrix<double, W, H> reduce(const Matrix<T, W, H>& matrix) {

//     Matrix<T, W, H> output = matrix;

//     for (std::size_t i_pivot = 0; i_pivot < W && i_pivot < H; ++i_pivot) {

//         for (std::size_t i_rows = i_pivot; i_rows < H; ++i_rows) {

//             float pivot = output[i_pivot][i_rows];

//             if (pivot != 0) {

//                 for (std::size_t i_columns = 0; i_columns < W; ++i_columns) {

//                     output[i_columns][i_rows] /= pivot;
//                 }
//             }
//         }

//         for (std::size_t i_rows = 0; i_rows < H; ++i_rows) {

//             float factor = output[i_pivot][i_rows];

//             if (i_rows != i_pivot) {

//                 for (std::size_t i_columns = 0; i_columns < W; ++i_columns) {

//                     output[i_columns][i_rows] -= output[i_columns][i_pivot] * factor;
//                 }
//             }
//         }

//         for (std::size_t i_columns = 0; i_columns < W; ++i_columns) {

//             for (std::size_t i_rows = 0; i_rows < H; ++i_rows) {

//                 if (output[i_columns][i_rows] < 0.0001 && output[i_columns][i_rows] > -0.0001)
//                     output[i_columns][i_rows] = 0;
//             }
//         }
//     }

//     return output;
// }

template <ScalarTypes T1, ScalarTypes T2, std::size_t W, std::size_t H>
Matrix<std::common_type_t<T1, T2>, W, H> add(const Matrix<T1, W, H>& a, const Matrix<T2, W, H>& b) {

    Matrix<std::common_type_t<T1, T2>, W, H> output;
    for (std::size_t i_columns = 0; i_columns < W; ++i_columns) {

        for (std::size_t i_rows = 0; i_rows < H; ++i_rows)
            output[i_columns][i_rows] = a[i_columns][i_rows] + b[i_columns][i_rows];
    }
    return output;
}

template <ScalarTypes T1, ScalarTypes T2, std::size_t W, std::size_t H>
Matrix<std::common_type_t<T1, T2>, W, H> subtract(const Matrix<T1, W, H>& a, const Matrix<T2, W, H>& b) {

    Matrix<std::common_type_t<T1, T2>, W, H> output;
    for (std::size_t i_columns = 0; i_columns < W; ++i_columns) {

        for (std::size_t i_rows = 0; i_rows < H; ++i_rows)
            output[i_columns][i_rows] = a[i_columns][i_rows] - b[i_columns][i_rows];
    }
    return output;
}

template <ScalarTypes T1, ScalarTypes T2, std::size_t W, std::size_t H, std::size_t D>
Matrix<std::common_type_t<T1, T2>, W, H> compose(const Matrix<T1, D, H>& a, const Matrix<T2, W, D>& b) {

    Matrix<std::common_type_t<T1, T2>, W, H> output;
    for (std::size_t i_columns = 0; i_columns < W; ++i_columns) {

        for (std::size_t i_rows = 0; i_rows < H; ++i_rows)
            output[i_columns][i_rows] = dotProduct(row(a, i_rows), column(b, i_columns));
    }
    return output;
}