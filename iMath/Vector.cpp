#include "Vector.h"
#include <cmath>
#include <stdlib.h>
#include <iostream>
#include <vector>
using namespace std;

Vector::Vector(std::initializer_list<double> l) {
    for(auto x : l) {
        this->elem.push_back(x);
    }
}

Vector::Vector(int len) {
    this->elem.assign(len, 0);
}

Vector::Vector(double *values, int len) {
    for(int i = 0; len > i; i++) {
        this->elem.push_back(values[i]);
        this->elem[i] = values[i];
    }
}

Vector::~Vector() {
}

size_t Vector::size() {
    return this->elem.size();
}

const size_t Vector::size() const{
    return this->elem.size();
}

double& Vector::operator[](int n) {
    if (n >= this->size() || n < 0) {
        std::cout << "ERROR: Out of Bound\n";
        exit(1);
    }

    return this->elem[n];
}

const double& Vector::operator[](int n) const {
    if (n >= this->size() || n < 0) {
        std::cout << "ERROR: Out of Bound\n";
        exit(1);
    }

    return this->elem[n];
}

std::ostream& operator<<(std::ostream& os, Vector& v) {
    os << '[';
    for (int i = 0; i < v.size() - 1; i++) {
        os << v[i] << ", ";
    }

    os << v[v.size() - 1] << ']';
    return os;
}

std::ostream& operator<<(std::ostream& os, const Vector& v) {
    os << '[';
    for (int i = 0; i < v.size() - 1; i++) {
        os << v[i] << ", ";
    }

    os << v[v.size() - 1] << ']';
    return os;
}


double Vector::operator*(const Vector &other) {
    if(this->size() != other.size()) {
        std::cout << "ERROR: Vector size mismatch" << '\n';
        exit(1);
    }
    double ans = 0;
    for (size_t i = 0; i < this->size(); i++) {
        ans += (*this)[i] * other[i];
    }
    return ans;
}

Vector& Vector::operator+(const Vector &other) {
    Vector *ans = new Vector;
    if(this->size() != other.size()) {
        std::cout << "ERROR: Vector size mismatch" << '\n';
        exit(1);
    }
    double tmp = 0;
    for (size_t i = 0; i < this->size(); i++) {
        tmp = (*this)[i] + other[i];
        ans->elem.push_back(tmp);
    }
    return *ans;
}

Vector& Vector::operator-(const Vector &other) {
    Vector *ans = new Vector;
    if(this->size() != other.size()) {
        std::cout << "ERROR: Vector size mismatch" << '\n';
        exit(1);
    }
    double tmp = 0;
    for (size_t i = 0; i < this->size(); i++) {
        tmp = (*this)[i] - other[i];
        ans->elem.push_back(tmp);
    }
    return *ans;
}

Vector& Vector::operator^(const Vector &other) {
    if(this->size() != 3 || other.size() != 3) {
        std::cout << "ERROR: Vector size mismatch" << '\n';
        exit(1);
    }
    Vector *ans = new Vector;
    ans->elem.push_back((*this)[1] * other[2] - (*this)[2] * other[1]);
    ans->elem.push_back((*this)[2] * other[0] - (*this)[0] * other[2]);
    ans->elem.push_back((*this)[0] * other[1] - (*this)[1] * other[0]);

    return *ans;
}

Vector& operator*(const int& mult, Vector& v) {
    Vector *ans = new Vector;

    for(auto i : v.elem) {
        ans->elem.push_back(mult * i);
    }

    return *ans;
}

Vector& operator*(const double& mult, Vector& v) {
    Vector *ans = new Vector;

    for(auto i : v.elem) {
        ans->elem.push_back(mult * i);
    }

    return *ans;
}
