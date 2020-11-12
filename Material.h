#ifndef __Material__H
#define __Material__H
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <utility>
#include "Matrix.h"
#include "Vector.h"
using namespace std;

struct Node {
    double x = 0;
    double moment = 0, react = 0;
};

class Moment {
public:
    double u(double x);
    Moment() : m(0), x(0) {}
    Moment(double m, double x) : m(m), x(x) {}
    double magnitude();
    double bodyMoment(double x);
    double m, x;
};

class cLoad {
public:
    double u(double x);
    cLoad() : q(0), x(0) {}
    cLoad(double q, double x) : q(q), x(x) {}
    double q, x;
    double magnitude();
    double shearForce(double x);
};

class dLoad {
public:
    double q1, q2, x1, x2;
    double val(double x);
    double magnitude() ;
    double shearForce(double x);
};

struct Beam {
public:
    Node *from, *to;
    double E, I;
    double len;
    pair<double, double> angle{0.0, 0.0};
};

struct Mesh {
    double x = 0;
    double sfd = 0, bmd = 0;
};

class Structure {
public:
    ~Structure() {}
    void moment3();
    int nodeNum, beamNum, cfNum, dfNum, mNum, meshNum;
    double totalLen;
    vector<Node> node;
    vector<Beam> beam;
    vector<cLoad> cf;
    vector<dLoad> df;
    vector<Moment> moment;
    vector<cLoad> reaction;

};
#endif
