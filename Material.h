#ifndef __Material__H
#define __Material__H
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <utility>
#include "Matrix.h"
#include "Vector.h"

struct Node {
    double x = 0;
    double moment = 0, react = 0;
};
std::istream& operator>>(std::istream& is, Node& n);

class Force {
    virtual double magnitude() {};
    virtual double shearForce(double x) {};
};

class Moment {
private:
    double u(double x);
public:
    Moment() : m(0), x(0) {}
    Moment(double m, double x) : m(m), x(x) {}
    double magnitude();
    double bodyMoment(double x);
    double m, x;
};
std::istream& operator>>(std::istream& is, Moment& m);

class cForce : Force {
public:
    cForce() : q(0), x(0) {}
    cForce(double q, double x) : q(q), x(x) {}
    double q, x;
    double magnitude() override;
    double shearForce(double x) override;
private:
    double u(double x);
};
std::istream& operator>>(std::istream& is, cForce& cf);
std::ostream& operator<<(std::ostream& os, cForce& cf);

class dForce : Force {
public:
    double q1, q2, x1, x2;
    double val(double x);
    double magnitude() override;
    double shearForce(double x) override;

};
std::istream& operator>>(std::istream& is, dForce& df);
std::ostream& operator<<(std::ostream& os, dForce& df);

struct Beam {
    Node *from, *to;
    double E, I;
    double len;
    std::pair<double, double> angle{0.0, 0.0};
};
std::istream& operator>>(std::istream& is, Beam& b);

struct Mesh {
    double x = 0;
    double sfd = 0, bmd = 0;
};
std::ostream& operator<<(std::ostream& os, Mesh& m);

class Structure {
public:
    ~Structure() {}
    void ls();
    void moment3();
    void calcGraph(std::vector<Mesh>& mesh);
    int nodeNum, beamNum, cfNum, dfNum, mNum, meshNum;
    double totalLen;
    std::vector<Node> node;
    std::vector<Beam> beam;
    std::vector<cForce> cf;
    std::vector<dForce> df;
    std::vector<Moment> moment;
    std::vector<cForce> reaction;

};
std::istream& operator>>(std::istream& is, Structure& st);
#endif
