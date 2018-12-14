#ifndef __Material__H
#define __Material__H
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "iMath.h"

class Force {
    virtual double magnitude() {};
    virtual double shearForce(double x) {};
};

struct Node {
    double x = 0;
    double moment = 0, angle = 0;
};
std::istream& operator>>(std::istream& is, Node& n);

class cForce : Force {
public:
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
    double magnitude() override;
    double shearForce(double x) override;

};
std::istream& operator>>(std::istream& is, dForce& df);
std::ostream& operator<<(std::ostream& os, dForce& df);

struct Beam {
    double E, I;
    double len;

};
std::istream& operator>>(std::istream& is, Beam& b);

struct Mesh {
    double x = 0;
    double sfd = 0, bmd = 0;

};
std::ostream& operator<<(std::ostream& os, dForce& df);

class Structure {
public:
    void ls();
    int nodeNum, beamNum, cfNum, dfNum, meshNum;
    double totalLen;
    std::vector<Node> node;
    std::vector<Beam> beam;
    std::vector<cForce> cf;
    std::vector<dForce> df;
};
std::istream& operator>>(std::istream& is, Structure& st);
#endif
