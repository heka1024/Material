#ifndef __Material__H
#define __Material__H
#include <iostream>
#include <fstream>
#include <string>
#include "iMath.h"
using namespace std;

class Force {
public:
    double valFrom, valTo;
    double from, to;
    bool isDistributed;

    double magnitude();
    double u(double x);
};

struct Beam {
    double E, I;
};

#endif
