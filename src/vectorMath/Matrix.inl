///* MATRIX *///

// CONSTRUCTORS AND DESTRUCTOR //

template <typename T>
Matrix<T>::Matrix(unsigned int width, unsigned int height) {

    m_width = width;
    m_height = height;
    if (m_width == 0)
        arr_data = NULL;
    else
        arr_data = new Vector<T>[m_width];
    for (unsigned int i_col = 0; i_col < m_width; ++i_col)
        arr_data[i_col] = Vector<T>(m_height);
    return;
}

template <typename T>
Matrix<T>::Matrix(unsigned int width, unsigned int height, const T** values) {

    m_width = width;
    m_height = height;
    if (m_width == 0)
        arr_data = NULL;
    else
        arr_data = new Vector<T>[m_width];
    for (unsigned int i_col = 0; i_col < m_width; ++i_col)
        arr_data[i_col] = Vector<T>(m_height);
    
    for (unsigned int i_col = 0; i_col < m_width; ++i_col) {

        for (unsigned int i_row = 0; i_row < m_height; ++i_row)
            arr_data[i_col][i_row] = values[i_col][i_row];
    }
    return;
}

template <typename T>
Matrix<T>::Matrix(const Matrix& matrix) {

    m_width = matrix.getWidth();
    m_height = matrix.getHeight();
    if (m_width == 0)
        arr_data = NULL;
    else
        arr_data = new Vector<T>[m_width];
    for (unsigned int i_col = 0; i_col < m_width; ++i_col)
        arr_data[i_col] = Vector<T>(m_height);

    for (unsigned int i_col = 0; i_col < m_width; ++i_col) {

        for (unsigned int i_row = 0; i_row < m_height; ++i_row)
            arr_data[i_col][i_row] = matrix[i_col][i_row];
    }
    return;
}

template <typename T>
Matrix<T>::~Matrix() {

    if (arr_data != NULL)
        delete [] arr_data;
    return;
}

// ACCESSORS //

template <typename T>
unsigned int Matrix<T>::getWidth() const {

    return m_width;
}

template <typename T>
unsigned int Matrix<T>::getHeight() const {

    return m_height;
}

template <typename T>
Matrix<T> Matrix<T>::getTransposedMatrix() const {

    Matrix<T> output(m_height, m_width);
    for (unsigned int i_col = 0; i_col < m_width; ++i_col) {

        for (unsigned int i_row = 0; i_row < m_height; ++i_row)
            output[i_row][i_col] = arr_data[i_col][i_row];
    }
    return output;
}

template <typename T>
Matrix<T> Matrix<T>::getReducedMatrix() const {

    Matrix<T> output(*this);

    unsigned pivotRow = 0;
    unsigned pivotColumn = 0;

    for (unsigned int i_row = 0; i_row < m_height; ++i_row) {

        for (unsigned int i_col = 0; i_col < m_width; ++i_col) {

            if (i_row == pivotRow) {
            
            
            }
            else {
            
            
            }
        }
    }

    return output;
}

template <typename T>
Matrix<T> Matrix<T>::getInverseMatrix() const {


}

template <typename T>
double Matrix<T>::getDeterminant() const {


}

template <typename T>
Vector<T> Matrix<T>::getRow(unsigned int index) const {

    Vector<T> output(m_width);
    for (unsigned int i_col = 0; i_col < m_width; ++i_col)
        output[i_col] = arr_data[i_col][index];
    return output;
}

template <typename T>
Vector<T> Matrix<T>::getColumn(unsigned int index) const {

    return arr_data[index];
}

template <typename T>
Matrix<T> Matrix<T>::scale(unsigned int scalar) const {

    Matrix<T> output(m_width, m_height);
    for (unsigned int i_col = 0; i_col < m_width; ++i_col) {

        for (unsigned int i_row = 0; i_row < m_height; ++i_row)
            output[i_col][i_row] = arr_data[i_col][i_row] * scalar;
    }
    return output;
}

template <typename T>
Vector<T> Matrix<T>::applyToVector(const Vector<T>& vector) const {

    return vector.applyMatrix(*this);
}


template <typename T>
Matrix<T> Matrix<T>::compose(const Matrix<T>& matrix) const {


}

template <typename T>
Matrix<T> Matrix<T>::add(const Matrix<T>& matrix) const {

    Matrix<T> output(m_width, m_height);
    for (unsigned int i_col = 0; i_col < m_width; ++i_col) {

        for (unsigned int i_row = 0; i_row < m_height; ++i_row)
            output[i_col][i_row] = arr_data[i_col][i_row] + matrix[i_col][i_row];
    }
    return output;
}

template <typename T>
Matrix<T> Matrix<T>::subtract(const Matrix<T>& matrix) const {

    Matrix<T> output(m_width, m_height);
    for (unsigned int i_col = 0; i_col < m_width; ++i_col) {

        for (unsigned int i_row = 0; i_row < m_height; ++i_row)
            output[i_col][i_row] = arr_data[i_col][i_row] - matrix[i_col][i_row];
    }
    return output;
}