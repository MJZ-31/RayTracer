#include <array>
#include <concepts>
#include <initializer_list>
#include <cmath>

// Define the concept of real numbers, which includes any integral or floating point numerical type.
template <typename T>
concept ScalarTypes = std::integral<T> || std::floating_point<T>;

/*
    Vector: A set of scalars which define a position in a Euclidian Space. Represented as a 1-dimensional arrays of real numbers.

    T: The type of the scalars. Must lie within the constraints of ScalarTypes.
    D: The number of scalars, or the number of dimensions the Vector has.
*/
template <ScalarTypes T, std::size_t D>
class Vector : public std::array<T, D> {

    public:

        /*
            Constructs a zero Vector (a Vector with all scalars set to 0).

                                       [0]
            Vector<int, 3> vector; --> [0]
                                       [0]
        */
        Vector();

        /*
            Constructs a Vector with its scalars set to the values in the given initializer_list.

            list: The values which the scalars are set to.

                - The size of the list is expected to be equal to D, except in either cases of
                    - An empty list, which will result in a zero Vector [B.].
                    - A list with only one value, which  will result in all scalars being set to that value [C.].
                
                - If they are not equal, missing values will be initialized to 0 [D.] and excess values will be ignored [E.].

                                                      [1]
            A. Vector<int, 3> vector { 1, 2, 3 }; --> [2]
                                                      [3]

                                             [0]
            B. Vector<int, 3> vector {}; --> [0]
                                             [0]

                                                [1]
            C. Vector<int, 3> vector { 1 }; --> [1]
                                                [1]

                                                   [2]
            D. Vector<int, 3> vector { 1, 2 }; --> [2]
                                                   [0]

                                                        [1]
            E. Vector<int, 3> vector { 1, 2, 3, 4 } --> [2]
                                                        [3]
        */
        template <std::convertible_to<T> U>
        Vector(const std::initializer_list<U>& list);

        /*
            Constructs a copy of the given Vector.

            vector: Vector to copy.
        */
        template <std::convertible_to<T> U>
        Vector(const Vector<U, D>& vector);
};

/*
    Matrix: A rectangular array of scalars, used to perform tranformations on Vectors.

    T: The type of the scalars. Must lie within the constraints of ScalarTypes.
    W: The width of the Matrix.
    H: The height of the Matrix.
*/
template <ScalarTypes T, std::size_t W, std::size_t H>
class Matrix : public std::array<Vector<T, H>, W> {

    public:

        /*
            Constructs a zero Matrix (a Matrix with all scalars set to 0).

                                          [0 0 0]
            Matrix<int, 3, 3> matrix; --> [0 0 0]
                                          [0 0 0]
        */
        Matrix();

        /*
            Constructs a Matrix with its scalars set to the values in the given initializer_list.

            list: The values which the scalars are set to. Each Vector within the list is a column in the Matrix.

                - The size of the list is expected to be equal to W [A.], except in either case of
                    - An empty list, which will result in a zero Matrix [B.].
                    - A list with only one column, which  will result in all columns being set to those values [C. and E.].

                - The size of every column is expected to be equal to H [A.], except in either case of
                    - An empty column, which will result in a column of all zeroes [D.].
                    - A column with only one value, which will result in all values in that column being set to that value [E. and F.].

                - Otherwise, missing values will be initialized to 0 [G.] and excess values will be ignored [H.].

                                                                                       [1 4 7]
            A. Matrix<int, 3, 3> matrix { { 1, 2, 3 }, { 4, 5, 6 }, { 7, 8, 9 } }; --> [2 5 8]
                                                                                       [3 6 9]

                                                [0 0 0]
            B. Matrix<int, 3, 3> matrix {}; --> [0 0 0]
                                                [0 0 0]

                                                             [1 1 1]
            C. Matrix<int, 3, 3> matrix { { 1, 2, 3 } }; --> [2 2 2]
                                                             [3 3 3]

                                                                              [0 1 4]
            D. Matrix<int, 3, 3> matrix { {}, { 1, 2, 3 }, { 4, 5, 6 } }; --> [0 2 5]
                                                                              [0 3 6]
             
                                                       [2 2 2]
            E. Matrix<int, 3, 3> matrix { { 2 } }; --> [2 2 2]
                                                       [2 2 2]
        
                                                                     [1 2 3]
            F. Matrix<int, 3, 3> matrix { { 1 }, { 2 }, { 3 } }; --> [1 2 3]
                                                                     [1 2 3]

                                                                    [1 3 0]
            G. Matrix<int, 3, 3> matrix { { 1, 2 }, { 3, 4 } }; --> [2 4 0]
                                                                    [0 0 0]

                                                                                [1 5]
            H. Matrix<int, 2, 3> matrix { { 1, 2, 3, 4 }, { 5, 6, 7, 8 } }; --> [2 6]
                                                                                [3 7]
        */
        template <std::convertible_to<T> U>
        Matrix(const std::initializer_list<std::initializer_list<U>>& list);

        /*
            Constructs a copy of the given Matrix.

            matrix: Matrix to copy.
        */
        template <std::convertible_to<T> U>
        Matrix(const Matrix<U, W, H>& matrix);
};

// Calculation Function Prototypes //

/*
    Calculates the length of the given Vector.

    T: The scalar type of the given Vector.
    D: The dimension of the given Vector.

    vector: The Vector to calculate the length of.
    return: The length of the Vector.

    No equivalent operator overloads.
*/
template <ScalarTypes T, std::size_t D>
double length(const Vector<T, D>& vector);

/*
    Normalizes a given Vector.

    T: The scalar type of the given Vector.
    D: The dimension of the given Vector and the returned Vector.

    vector: The Vector to normalize.
    return: The normalized Vector.

    No equivalent operator overloads.
*/
template <std::convertible_to<double> T, std::size_t D>
Vector<double, D> normalize(const Vector<T, D>& vector);

/*
    Adds together two given Vectors.

    T1: The scalar type of given Vector a.
    T2: The scalar type of given Vector b.
    D: The dimension of a, b, and the returned Vector.

    a, b: The Vectors to add.
    return: The resulting Vector of adding a and b.

    Equivalent operator overload: operator+(a, b)
*/
template <ScalarTypes T1, ScalarTypes T2, std::size_t D>
Vector<std::common_type_t<T1, T2>, D> add(const Vector<T1, D>& a, const Vector<T2, D>& b);

/*
    Subtracts two given Vectors.

    T1: The scalar type of given Vector a.
    T2: The scalar type of given Vector b.
    D: The dimension of a, b, and the returned Vector.

    a, b: The Vectors to subtract.
    return: The resulting Vector of subtracting b from a.

    Equivalent operator overload: operator-(a, b)
*/
template <ScalarTypes T1, ScalarTypes T2, std::size_t D>
Vector<std::common_type_t<T1, T2>, D> subtract(const Vector<T1, D>& a, const Vector<T2, D>& b);

/*
    Scales the length of a given Vector by a given factor.

    T1: The scalar type of the given Vector.
    T2: The type of the given factor.
    D: The dimension of the given Vector and the returned Vector.

    vector: The Vector to scale.
    factor: The factor to scale the Vector's length by.
    return: The scaled Vector.

    Equivalent operator overloads: operator*(vector, factor) or operator*(factor, vector)
*/
template <ScalarTypes T1, ScalarTypes T2, std::size_t D>
Vector<std::common_type_t<T1, T2>, D> scale(const Vector<T1, D>& vector, T2 factor);

/*
    Calculates the dot product between two given Vectors.

    T1: The scalar type of given Vector a.
    T2: The scalar type of given Vector b.
    D: The dimension of a and b.

    a, b: The Vectors to calculate the dot product of.
    return: The result of the dot product.

    Equivalent operator overload: operator*(a, b)
*/
template <ScalarTypes T1, ScalarTypes T2, std::size_t D>
std::common_type_t<T1, T2> dotProduct(const Vector<T1, D>& a, const Vector<T2, D>& b);

/*
    Calculates the cross product between two given Vectors.

    T1: The scalar type of given Vector a.
    T2: The scalar type of given Vector b.

    a, b: The Vectors to calculate the dot product of.
    return: The result of the cross product, a Vector which is perpendicular to both a and b.

        - All three must be 3-dimensional Vectors.

    Equivalent operator overload: operator%(a, b)
*/
template <ScalarTypes T1, ScalarTypes T2>
Vector<std::common_type_t<T1, T2>, 3> crossProduct(const Vector<T1, 3>& a, const Vector<T2, 3>& b);

/*
    Transforms a given Vector using a given Matrix.

    T1: The scalar type of the given Vector.
    T2: The scalar type of the given Matrix.
    D1: The dimension of the given Vector and the width of the given Matrix.
    D2: The dimension of the returned Vector and the height of the given Matrix.

    vector: The Vector to be transformed.
    matrix: The tranformation Matrix to be used.
    return: The transformed Vector.

    Equivalent operator overloads: operator*(vector, matrix) or operator*(matrix, vector)
*/
template <ScalarTypes T1, ScalarTypes T2, std::size_t D1, std::size_t D2>
Vector<std::common_type<T1, T2>, D2> transform(const Vector<T1, D1>& vector, const Matrix<T2, D1, D2>& matrix);

/*
    Retrieves a row from a given Matrix at a given row index.

    T: The scalar type of the given Matrix and the returned Vector.
    W: The width of the given Matrix and the dimension of the returned Vector.
    H: The height of the given Matrix.

    matrix: The Matrix to retrieve the row from.
    index: The index of the row, counting from the top to the bottom of the given Matrix.
    return: A Vector with the row of scalars.

    No equivalent operator overloads.
*/
template <ScalarTypes T, std::size_t W, std::size_t H>
Vector<T, W> row(const Matrix<T, W, H>& matrix, std::size_t index);

/*
    Retrieves a column from a given Matrix at a given column index.

    T: The scalar type of the given Matrix and the returned Vector.
    W: The width of the given Matrix.
    H: The height of the given Matrix and the dimension of the returned Vector.

    matrix: The Matrix to retrieve the column from.
    index: The index of the column, counting from the left to the right of the given Matrix.
    return: A Vector with the column of scalars.

    No equivalent operator overloads.
*/
template <ScalarTypes T, std::size_t W, std::size_t H>
Vector<T, H> column(const Matrix<T, W, H>& matrix, std::size_t index);

/*
    Calculates the transpose of a given Matrix.

    T: The scalar type of the given Matrix and the returned Matrix.
    W: The width of the given Matrix and the height of the returned Matrix.
    H: The height of the given Matrix and the width of the returned Matrix.

    matrix: The Matrix to calculate the transpose of.
    return: The transposed Matrix.

    No equivalent operator overloads.
*/
template <ScalarTypes T, std::size_t W, std::size_t H>
Matrix<T, H, W> transpose(const Matrix<T, W, H>& matrix);

// template <ScalarTypes T, std::size_t W, std::size_t H>
// Matrix<double, W, H> reduce(const Matrix<T, W, H>& matrix);

/*
    Adds together two given Matrices.

    T1: The scalar type of given Matrix a.
    T2: The scalar type of given Matrix b.
    W: The width of a, b, and the returned Matrix.
    H: The height of a, b, and the returned Matrix.

    a, b: The Matrices to add.
    return: The resulting Matrix of adding a and b.

    Equivalent operator overload: operator+(a, b)
*/
template <ScalarTypes T1, ScalarTypes T2, std::size_t W, std::size_t H>
Matrix<std::common_type_t<T1, T2>, W, H> add(const Matrix<T1, W, H>& a, const Matrix<T2, W, H>& b);

/*
    Subtracts two given Matrices.

    T1: The scalar type of given Matrix a.
    T2: The scalar type of given Matrix b.
    W: The width of a, b, and the returned Matrix.
    H: The height of a, b, and the returned Matrix.

    a, b: The Matrices to subtract.
    return: The resulting Matrix of subtracting b from a.

    Equivalent operator overload: operator-(a, b)
*/
template <ScalarTypes T1, ScalarTypes T2, std::size_t W, std::size_t H>
Matrix<std::common_type_t<T1, T2>, W, H> subtract(const Matrix<T1, W, H>& a, const Matrix<T2, W, H>& b);

/*
    Composes two given Matrices together.

    T1: The scalar type of given Matrix a.
    T2: The scalar type of given Matrix b.
    W: The width of b and the returned Matrix.
    H: The height of a and the returned Matrix.
    D: The width of a and the height of b.

    a, b: The Matrices to compose.
    return: The resulting composed Matrix.

    Equivalent operator overload: operator*(a, b)
*/
template <ScalarTypes T1, ScalarTypes T2, std::size_t W, std::size_t H, std::size_t D>
Matrix<std::common_type_t<T1, T2>, W, H> compose(const Matrix<T1, D, H>& a, const Matrix<T2, W, D>& b);

// Operators Overloads //

template <ScalarTypes T1, ScalarTypes T2, std::size_t D>
Vector<std::common_type_t<T1, T2>, D> operator+(const Vector<T1, D>& a, const Vector<T2, D>& b) { return add(a, b); };

template <ScalarTypes T1, ScalarTypes T2, std::size_t W, std::size_t H>
Matrix<std::common_type_t<T1, T2>, W, H> operator+(const Matrix<T1, W, H>& a, const Matrix<T2, W, H>& b) { return add(a, b); }

template <ScalarTypes T1, ScalarTypes T2, std::size_t D>
Vector<std::common_type_t<T1, T2>, D> operator-(const Vector<T1, D>& a, const Vector<T2, D>& b) { return subtract(a, b); }

template <ScalarTypes T1, ScalarTypes T2, std::size_t W, std::size_t H>
Matrix<std::common_type_t<T1, T2>, W, H> operator-(const Matrix<T1, W, H>& a, const Matrix<T2, W, H>& b) { return subtract(a, b); }

template <ScalarTypes T1, ScalarTypes T2, std::size_t D>
Vector<std::common_type_t<T1, T2>, D> operator*(const Vector<T1, D>& vector, T2 scalar) { return scale(vector, scalar); }

template <ScalarTypes T1, ScalarTypes T2, std::size_t D>
Vector<std::common_type_t<T1, T2>, D> operator*(T1 scalar, const Vector<T2, D>& vector) { return scale(vector, scalar); }

template <ScalarTypes T1, ScalarTypes T2, std::size_t D>
std::common_type_t<T1, T2> operator*(const Vector<T1, D>& a, const Vector<T2, D>& b) { return dotProduct(a, b); };

template <ScalarTypes T1, ScalarTypes T2, std::size_t D1, std::size_t D2>
Vector<std::common_type_t<T1, T2>, D2> operator*(const Vector<T1, D1>& vector, const Matrix<T2, D1, D2>& matrix) { return transform(vector, matrix); };

template <ScalarTypes T1, ScalarTypes T2, std::size_t D1, std::size_t D2>
Vector<std::common_type_t<T1, T2>, D2> operator*(const Matrix<T1, D1, D2>& matrix, const Vector<T2, D1>& vector) { return transform(vector, matrix); };

template <ScalarTypes T1, ScalarTypes T2, std::size_t W, std::size_t H, std::size_t D>
Matrix<std::common_type_t<T1, T2>, W, H> operator*(const Matrix<T1, D, H>& a, const Matrix<T2, W, D>& b) { return compose(a, b); }

template <ScalarTypes T1, ScalarTypes T2>
Vector<std::common_type_t<T1, T2>, 3> operator%(const Vector<T1, 3>& a, const Vector<T2, 3>& b) { return crossProduct(a, b); };

#include "vectorMath.inl"