#ifndef __iVECTOR__H
#define __iVECTOR__H
#include <iostream>
#include <vector>


class Vector {
public:
    std::vector<double> elem;

    // Constructors
    Vector() : elem(0) {};
    Vector(std::initializer_list<double> l) {
        for (auto x : l) {
            this->elem.push_back(x);
        }
    }
    Vector(int len) {
        this->elem.assign(len, 0);
    }

    // Desctructor
    ~Vector() {}

    int size() { 
        return this->elem.size();
    }


    double& operator[](int n) {
        return this->elem[n];
    }

    Vector& operator+(Vector &other) {
        Vector *ans = new Vector;
        double tmp = 0;
        for (int i = 0; i < this->size(); i++) {
            tmp = (*this)[i] + other[i];
            ans->elem.push_back(tmp);
        }
        return *ans;
    }
    Vector& operator-(Vector &other) {
        Vector *ans = new Vector;
        double tmp = 0;
        for (int i = 0; i < this->size(); i++) {
            tmp = (*this)[i] - other[i];
            ans->elem.push_back(tmp);
        }
        return *ans;
    }
    Vector& operator^(Vector &other) {
        Vector *ans = new Vector;
        ans->elem.push_back((*this)[1] * other[2] - (*this)[2] * other[1]);
        ans->elem.push_back((*this)[2] * other[0] - (*this)[0] * other[2]);
        ans->elem.push_back((*this)[0] * other[1] - (*this)[1] * other[0]);

        return *ans;
    }
    
    double operator*(Vector &other) {
        double ans = 0;
        for (int i = 0; i < this->size(); i++) {
            ans += (*this)[i] * other[i];
        }
        return ans;
    }

    Vector& operator*(const double mult) {
        Vector *ans = new Vector;

        for(auto i : this->elem) {
            ans->elem.push_back(mult * i);
        }

        return *ans;
    }
};



#endif
