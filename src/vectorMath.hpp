#include <array>
#include <concepts>
#include <initializer_list>

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
        Vector() {

            for (unsigned int i = 0; i < D; ++i)
                (*this)[i] = static_cast<T>(0);
        }

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
        template<std::convertible_to<T> U>
        Vector(const std::initializer_list<U>& list) : Vector() {

            if (list.size() == 0)
                return;
            else if (list.size() == 1) {

                for (unsigned int i = 0; i < D; ++i)
                    (*this)[i] = static_cast<T>(*list.begin());
            }
            else {

                for (unsigned int i = 0; i < D && i < list.size(); ++i)
                    (*this)[i] = static_cast<T>(*(list.begin() + i));
            }
        }

        /*
            Constructs a copy of the given Vector.

            vector: Vector to copy.
        */
        template<std::convertible_to<T> U>
        Vector(const Vector<U, D>& vector) {

            for (unsigned int i = 0; i < D; ++i)
                (*this)[i] = static_cast<T>(vector[i]);
        }
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
        Matrix() {

            for (unsigned int i = 0; i < W; ++i)
                (*this)[i] = Vector<T, H>();
        }

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
        template<std::convertible_to<T> U>
        Matrix(const std::initializer_list<std::initializer_list<U>>& list) : Matrix() {

            if (list.size() == 0)
                return;
            else if (list.size() == 1) {

                for (unsigned int i = 0; i < W; ++i)
                    (*this)[i] = Vector<T, H>(*list.begin());
            }
            else {

                for (unsigned int i = 0; i < W && i < list.size(); ++i)
                    (*this)[i] = Vector<T, H>(*(list.begin() + i));
            }
        }

        /*
            Constructs a copy of the given Matrix.

            matrix: Matrix to copy.
        */
        template<std::convertible_to<T> U>
        Matrix(const Matrix<U, W, H>& matrix) {

            for (unsigned int i = 0; i < W; ++i)
                (*this)[i] = matrix[i];
        }
};

// Calculates the length of the given Vector.
template <ScalarTypes T, std::size_t D>
double length(const Vector<T, D>& vector);

// Adds together Vector a and Vector b.
template <ScalarTypes T1, ScalarTypes T2, std::size_t D>
Vector<std::common_type_t<T1, T2>, D> add(const Vector<T1, D>& a, const Vector<T2, D>& b);

// Subtracts Vector b from Vector a.
template <ScalarTypes T1, ScalarTypes T2, std::size_t D>
Vector<std::common_type_t<T1, T2>, D> subtract(const Vector<T1, D>& a, const Vector<T2, D>& b);

// Scales the length of the given Vector by the given scalar.
template <ScalarTypes T1, ScalarTypes T2, std::size_t D>
Vector<std::common_type_t<T1, T2>, D> scale(const Vector<T1, D>& vector, T2 scalar);

// Calculates the dot product between Vector a and Vector b.
template <ScalarTypes T1, ScalarTypes T2, std::size_t D>
std::common_type_t<T1, T2> dotProduct(const Vector<T1, D>& a, const Vector<T2, D>& b);

// Calculates the cross product between Vector a and Vector b.
template <ScalarTypes T1, ScalarTypes T2>
Vector<std::common_type_t<T1, T2>, 3> crossProduct(const Vector<T1, 3>& a, const Vector<T2, 3>& b);

// Calculates the transformed Vector by applying the given Matrix to the given Vector.
template <ScalarTypes T1, ScalarTypes T2, std::size_t D1, std::size_t D2>
Vector<std::common_type_t<T1, T2>, D2> transform(const Vector<T1, D1>& vector, const Matrix<T2, D1, D2>& matrix);

// Calculates the determinant of the given Matrix.
template <ScalarTypes T, std::size_t W, std::size_t H>
T determinant(const Matrix<T, W, H>& matrix);

// Retreives a row from the given Matrix.
template <ScalarTypes T, std::size_t W, std::size_t H>
Vector<T, W> row(const Matrix<T, W, H>& matrix, unsigned int index);

// Retreives a column from the given Matrix.
template <ScalarTypes T, std::size_t W, std::size_t H>
Vector<T, H> column(const Matrix<T, W, H>& matrix, unsigned int index);

// Calculates the transpose Matrix of the given Matrix.
template <ScalarTypes T, std::size_t W, std::size_t H>
Matrix<T, H, W> transpose(const Matrix<T, W, H>& matrix);

// Calculates the row reduced Matrix of the given Matrix.
template <ScalarTypes T, std::size_t W, std::size_t H>
Matrix<T, W, H> reduce(const Matrix<T, W, H>& matrix);

// Calculates the inverse Matrix of the given Matrix.
template <ScalarTypes T, std::size_t W, std::size_t H>
Matrix<T, W, H> inverse(const Matrix<T, W, H>& matrix);

// Adds together Matrix a and Matrix b.
template <ScalarTypes T1, ScalarTypes T2, std::size_t W, std::size_t H>
Matrix<std::common_type_t<T1, T2>, W, H> add(const Matrix<T1, W, H>& a, const Matrix<T2, W, H>& b);

// Subtracts Matrix b from Matrix a.
template <ScalarTypes T1, ScalarTypes T2, std::size_t W, std::size_t H>
Matrix<std::common_type_t<T1, T2>, W, H> subtract(const Matrix<T1, W, H>& a, const Matrix<T2, W, H>& b);

// Composes Matrix a and Matrix b together.
template <ScalarTypes T1, ScalarTypes T2, std::size_t W, std::size_t H, std::size_t D>
Matrix<std::common_type_t<T1, T2>, W, H> compose(const Matrix<T1, D, H>& a, const Matrix<T2, W, D>& b);

// Operators

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