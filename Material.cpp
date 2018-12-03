#include <iostream>
#include <fstream>
#include <string>
#include "Material.h"
#include "iMath.h"

double Force::magnitude() {
    if(this->isDistributed) {
        return (this->to - this->from) * (this->valFrom + this->valTo) / 2;
    } else {    // concentrated
        return this->valFrom;
    }
}

double Force::u(double x) {
    if(isDistributed) {
        double a = this->from, b = this->to;
        double fa = this->valFrom, fb = this->valTo;
        if(x < a) {
            return 0;
        } else if(x < b) {
            double ans = fa * (x - a);
            ans += ((fb - fa) / (b - a)) * (x - a) * (x - a) / 2;
            return ans;
        } else {
            return (b - a) * (fa + fb) / 2;
        }
    } else {    // concentrated
        double a = this->from;
        if(x >= a) {
            return this->valFrom;
        } else {
            return 0;
        }
    }
}
