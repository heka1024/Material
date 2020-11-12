#include <cmath>
#include <iostream>
#include <vector>
#include "Vector.h"
#include "Matrix.h"

Matrix::Matrix(int rows, int cols) {
    for (int i = 0; i < rows; i++) {
        this->elem.push_back(Vector(cols));
    }
}

Matrix::Matrix(std::initializer_list<Vector> l) {
    for(auto x : l) {
        this->elem.push_back(x);
    }
}

Vector& Matrix::operator[](int n) {
    if (n >= this->rows() || n < 0) {
        std::cout << "ERROR: Out of Bound\n";
        exit(1);
    }

    return this->elem[n];
}


Matrix::Matrix( Matrix& origin) {
    for (int i = 0; i < origin.rows(); i++) {
        this->elem.push_back(origin[i]);
    }
}

Matrix& Matrix::operator+(Matrix& other) {
    if((this->rows() != other.rows()) || (this->cols() != other.cols())) {
        std::cout << "ERROR: Matrix size mismatch" << '\n';
        exit(1);
    }

    Matrix *pnew = new Matrix;
    for (int i = 0; i < this->rows(); i++) {
        Vector add = (*this)[i] + other[i];
        pnew->elem.push_back(add);
    }

    return *pnew;
}

Matrix& Matrix::operator-(Matrix& other) {
    if((this->rows() != other.rows()) || (this->cols() != other.cols())) {
        std::cout << "ERROR: Matrix size mismatch" << '\n';
        exit(1);
    }

    Matrix *pnew = new Matrix;
    for (int i = 0; i < this->rows(); i++) {
        Vector sub = (*this)[i] - other[i];
        pnew->elem.push_back(sub);
    }

    return *pnew;
}

Vector& Matrix::operator*( Vector& rhs) {
    Vector *ans = new Vector;

    for(auto row : this->elem) {
        double mult = row * rhs;
        ans->elem.push_back(mult);
    }

    return *ans;
}

Matrix& Matrix::transpose() {
    Matrix *trans = new Matrix(this->cols(), this->rows());

    for(int i = 0; i < this->rows(); i++) {
        for(int j = 0; j < this->cols(); j++) {
            (*trans)[j][i] = (*this)[i][j];
        }
    }
    return *trans;
}

Matrix& Matrix::operator*(Matrix& other) {
    if(this->cols() != other.rows()) {
        std::cout << "ERROR: Matrix mult error. Please check size of matrix.\n";
        exit(1);
    }
    Matrix *mult = new Matrix;
    Matrix rhs = other.transpose();
    for(auto a : this->elem) {
        Vector row;
        for(auto b : rhs.elem) {
            row.elem.push_back(a * b);
        }
        mult->elem.push_back(row);
    }

    return *mult;
}

void swap(double& a, double& b) {
    double tmp = a;
    a = b;
    b = tmp;
}

void swap(Vector& a, Vector& b) {
    a.elem.swap(b.elem);
}

Vector& Matrix::solve(Vector& colVec) {
    if (this->cols() != this->rows()) {
        std::cout << "ERROR: Matrix is not sqaure";
        exit(1);
    } else if (this->cols() != colVec.size()) {
        std::cout << "ERROR: Matrix and Vector size does not match.";
        exit(1);
    }

    int n = this->cols();
    Matrix A = (*this);
    Vector b = colVec;
    for (int k = 0; k < n; k++) {
        int m = k;
        for (int j = k + 1; j < n; j++) {
            if(std::abs(A[m][k]) < std::abs(A[j][k])) {
                m = j;
            }
        }
        
        if(A[m][k] == 0) {
            std::cout << "ERROR: No unique solutions" << '\n';
        } else {
            swap(A[m], A[k]);
            swap(b[m], b[k]);
        }
        if(A[n-1][n-1] == 0) {
            std::cout << "ERROR: No unique solutions" << '\n';
        } else {
            for (int j = k + 1; j < n; j++) {
                double factor = A[j][k] / A[k][k];
                A[j] = A[j] - A[k] * factor;
                b[j] = b[j] - factor * b[k];
            }
        }
    }

    Vector *ans = new Vector(n);
    for (int i = n - 1; i >= 0; i--) {
        Vector coe = A[i]; coe[i] = 0;
        (*ans)[i] = (b[i] - coe * (*ans)) / A[i][i];
    }

    return *ans;
}

Matrix& Matrix::inverse() {
    Matrix inv;
    int n = this->cols();
    for (int i = 0; i < n; i++) {
        Vector e(n); e[i] = 1;
        Vector tmp = this->solve(e);
        inv.elem.push_back(tmp);
    }
    Matrix *pnew = new Matrix(inv.transpose());
    return *pnew;
}