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
    Matrix(Matrix& origin);
    // Desctructor
    ~Matrix() {};

    // public methods
    int rows() { return elem.size(); }
    int cols() { return elem[0].size(); }
    Matrix& transpose();
    Matrix& inverse();
    Vector& solve(Vector& colVec);

    Matrix& operator+(Matrix& other);
    Matrix& operator-(Matrix& other);
    Matrix& operator*(Matrix& other);
    Matrix& operator*(double factor) {
        Matrix *pnew = new Matrix();
        for(auto x : this->elem) {
            pnew->elem.push_back(x * factor);
        }
        return *pnew;
    }
    Vector& operator*(Vector& rhs);
    Vector& operator[](int n);
};
#endif
