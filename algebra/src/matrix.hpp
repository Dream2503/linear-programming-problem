#pragma once
#include "matrix.hpp"

template <typename T>
class algebra::Matrix {
public:
    int row, column;
    std::vector<std::vector<T>> matrix;

    Matrix() = default;

    Matrix(const int row, const int column, const T& value = T()) : row(row), column(column), matrix(std::vector(row, std::vector(column, value))) {}

    Matrix(const std::vector<std::vector<T>>& matrix) : row(matrix.size()), column(matrix[0].size()), matrix(matrix) {}

    template <typename U>
    Matrix(const Matrix<U>& other) : Matrix(other.row, other.column) {
        for (int i = 0; i < row; i++) {
            for (int j = 0; j < column; j++) {
                matrix[i][j] = static_cast<T>(other[i, j]);
            }
        }
    }

    T determinant() const {
        assert(row == column);
        Matrix mat = echelon_form();
        T res = T(1);

        for (int i = 0; i < row; i++) {
            res *= mat[i, i];
        }
        return res;
    }

    Matrix echelon_form() const {
        Matrix res = *this;
        int pivot = 0;

        for (int cnt = 0; cnt < column && pivot < row; cnt++) {
            int pivot_row = pivot;

            while (pivot_row < row && res[pivot_row, cnt] == 0) {
                pivot_row++;
            }
            if (pivot_row == row) {
                continue;
            }
            if (pivot_row != pivot) {
                for (int j = 0; j < column; j++) {
                    std::swap(res[pivot_row, j], res[pivot, j]);
                }
            }
            for (int i = pivot + 1; i < row; i++) {
                T factor = res[i, cnt] / res[pivot, cnt];

                for (int j = cnt; j < column; j++) {
                    res[i, j] -= factor * res[pivot, j];
                }
            }
            pivot++;
        }
        return res;
    }

    Matrix inverse() const {
        assert(determinant() != 0);
        Matrix res(row, column), aug_matrix = make_identity(row);

        for (int i = 0; i < row; i++) {
            aug_matrix[i].insert(aug_matrix[i].begin(), matrix[i].begin(), matrix[i].end());
        }
        aug_matrix.column *= 2;
        aug_matrix = aug_matrix.echelon_form();

        for (int i = row - 1; i >= 0; i--) {
            T pivot = aug_matrix[i, i];
            assert(pivot != 0);

            for (int j = 0; j < row * 2; j++) {
                aug_matrix[i, j] /= pivot;
            }
            for (int k = 0; k < i; k++) {
                T factor = aug_matrix[k, i];

                for (int j = 0; j < row * 2; j++) {
                    aug_matrix[k, j] -= factor * aug_matrix[i, j];
                }
            }
        }
        for (int i = 0; i < row; i++) {
            for (int j = 0; j < row; j++) {
                res[i][j] = aug_matrix[i, j + row];
            }
        }
        return res;
    }

    static Matrix make_identity(const int row) {
        Matrix res(row, row);

        for (int i = 0; i < row; i++) {
            res[i, i] = T(1);
        }
        return res;
    }

    Matrix transpose() const {
        Matrix res(row, column);

        for (int i = 0; i < row; i++) {
            for (int j = 0; j < column; j++) {
                res[i, j] = matrix[j][i];
            }
        }
        return res;
    }

    template <typename U>
    Matrix& operator+=(const Matrix<U>& other) {
        assert(row == other.row && column == other.column);

        for (int i = 0; i < row; i++) {
            for (int j = 0; j < column; j++) {
                this->matrix[i][j] += other[i, j];
            }
        }
        return *this;
    }

    template <typename U, typename R = decltype(T() + U())>
    Matrix<R> operator+(const Matrix<U>& other) const {
        Matrix<R> res = *this;
        res += other;
        return res;
    }

    template <typename U>
        requires(!std::is_same_v<U, Matrix>)
    Matrix& operator+=(const U& other) {
        for (int i = 0; i < row; i++) {
            for (int j = 0; j < column; j++) {
                this->matrix[i][j] += other;
            }
        }
        return *this;
    }

    template <typename U, typename R = decltype(T() + U())>
    Matrix<R> operator+(const U& other) const {
        Matrix<R> res = *this;
        res += other;
        return res;
    }

    template <typename U>
    Matrix& operator-=(const Matrix<U>& other) {
        assert(row == other.row && column == other.column);

        for (int i = 0; i < row; i++) {
            for (int j = 0; j < column; j++) {
                this->matrix[i][j] -= other[i, j];
            }
        }
        return *this;
    }

    template <typename U, typename R = decltype(T() - U())>
    Matrix<R> operator-(const Matrix<U>& other) const {
        Matrix<R> res = *this;
        res -= other;
        return res;
    }

    template <typename U>
        requires(!std::is_same_v<U, Matrix>)
    Matrix& operator-=(const U& other) {
        for (int i = 0; i < row; i++) {
            for (int j = 0; j < column; j++) {
                this->matrix[i][j] -= other;
            }
        }
        return *this;
    }

    template <typename U, typename R = decltype(T() - U())>
    Matrix<R> operator-(const U& other) const {
        Matrix<R> res = *this;
        res -= other;
        return res;
    }

    template <typename U>
    Matrix& operator*=(const Matrix<U>& other) {
        *this = *this * other;
        return *this;
    }

    template <typename U, typename R = decltype(T() * U() + T() * U())>
    Matrix<R> operator*(const Matrix<U>& other) const {
        assert(column == other.row);
        Matrix<R> res(row, other.column);

        for (int i = 0; i < row; i++) {
            for (int j = 0; j < other.column; j++) {
                for (int k = 0; k < column; k++) {
                    res[i, j] += matrix[i][k] * other[k, j];
                }
            }
        }
        return res;
    }

    template <typename U>
        requires(!std::is_same_v<U, Matrix>)
    Matrix& operator*=(const U& other) {
        for (int i = 0; i < row; i++) {
            for (int j = 0; j < column; j++) {
                this->matrix[i][j] *= other;
            }
        }
        return *this;
    }

    template <typename U, typename R = decltype(T() * U())>
    Matrix<R> operator*(const U& other) const {
        Matrix<R> res = *this;
        res *= other;
        return res;
    }

    template <typename U>
        requires(!std::is_same_v<U, Matrix>)
    Matrix& operator/=(const U& other) {
        for (int i = 0; i < row; i++) {
            for (int j = 0; j < column; j++) {
                this->matrix[i][j] /= other;
            }
        }
        return *this;
    }

    template <typename U, typename R = decltype(T() / U(1))>
    Matrix<R> operator/(const U& other) const {
        Matrix<R> res = *this;
        res /= other;
        return res;
    }

    bool operator==(const Matrix& other) const { return matrix == other.matrix; }

    friend std::ostream& operator<<(std::ostream& out, const Matrix& other) {
        for (int i = 0; i < other.row; i++) {
            for (int j = 0; j < other.column; j++) {
                out << other[i, j] << ' ';
            }
            out << std::endl;
        }
        return out;
    }

    std::vector<T>& operator[](const int i) { return matrix[i]; }

    const std::vector<T>& operator[](const int i) const { return matrix[i]; }

    T& operator[](const int i, const int j) { return matrix[i][j]; }

    const T& operator[](const int i, const int j) const { return matrix[i][j]; }
};
