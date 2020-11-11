#ifndef __MATRIX__H
#define __MATRIX__H
#include <iostream>
#include <vector>
#include "Vector.h"

class Matrix {
public:
    std::vector<Vector> elem;

    // Constructors
    Matrix() : elem(0) {}
    Matrix(int rows, int cols);
    Matrix(std::initializer_list<Vector> l);
    Matrix(const Matrix& origin);
    // Desctructor
    ~Matrix() {};

    // public methods
    double rows() { return elem.size(); }
    const double rows() const { return elem.size(); };
    double cols() { return elem[0].size(); }
    const double cols() const { return elem[0].size(); };
    Matrix& pivoting(Matrix& P);
    double determinant();
    bool isSingular();
    Matrix& transpose();
    Matrix& inverse();
    Vector& solve(Vector& colVec);

    Matrix& operator+(Matrix& other);
    Matrix& operator-(Matrix& other);
    Matrix& operator*(Matrix& other);
    Vector& operator*(const Vector& rhs);
    Vector& operator[](int n);
    const Vector& operator[](int n) const;
};
Matrix& operator*(const int& factor, const Matrix& m);
Matrix& operator*(const double& factor, const Matrix& m);
std::ostream& operator<<(std::ostream& os, const Matrix& M);
std::ostream& operator<<(std::ostream& os, Matrix& M);
#endif
