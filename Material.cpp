#include <cmath>
#include <vector>
#include <utility>
#include <algorithm>
#include "Material.h"
#include "Matrix.hpp"
#include "Vector.h"

using namespace std;

double cLoad::magnitude() {
    return this->q;
}

double cLoad::u(double x) {
    return (x >= this->x);
}

double cLoad::shearForce(double x) {
    return this->u(x) * this->magnitude();
}

double Moment::magnitude() {
    return this->m;
}

double Moment::bodyMoment(double x) {
    return this->u(x) * this->magnitude();
}

double Moment::u(double x) {
    return (x >= this->x);
}

double dLoad::val(double x) {
    return this->q1 + (this->q2 - this->q1) * ((x - this->x1) / (this->x2 - this->x1));
}

double dLoad::magnitude() {
    double avgVal = (this->q1 + this->q2) / 2;
    return avgVal * (this->x2 - this->x1);
}

double dLoad::shearForce(double x) {
    if(x < this->x1) {
        return 0;
    } else if (x < this->x2) {
        return (x - this->x1) * (this->val(x) + q1) / 2;
    } else {
        return this->magnitude();
    }
}

template <typename T,typename U>
pair<T,U> operator+(const pair<T,U> & l,const pair<T,U> & r) {
    return {l.first + r.first, l.second + r.second};
}

template <typename T,typename U>
pair<T,U> operator-(const pair<T,U> & l,const pair<T,U> & r) {
    return {l.first - r.first, l.second - r.second};
}

pair<double, double> calcSquareAngle(Beam &beam, double q, double k1, double k2) {
    double c = beam.to->x - k2, a = k1 - beam.from->x, b = k2 - k1;
    double l = beam.to->x - beam.from->x, d = a + b;
    pair<double, double> ans{0, 0};

    double tmp = (double) (q * b) / (48 * beam.E * beam.I * l);
    ans.first -= tmp * (b + 2 * c) * (
        4 * pow(l, 2) - pow(b, 2) - pow(b + 2 * c, 2)
    );
    ans.second += tmp * (2 * d -b) * (
        2 * b * (a + d) + 4 * c * (l + d)
    );
    return ans;
}

pair<double, double> calcTriAngle(Beam &beam, double q, double k1, double k2, bool leftHigher) {
    double c = beam.to->x - k2, a = k1 - beam.from->x, b = k2 - k1;
    double l = beam.to->x - beam.from->x, d = a + b;
    double tmp = (q * b) / (12 * beam.E * beam.I * l);
    pair<double, double> ans(0, 0);
    if(leftHigher) {
        ans = calcSquareAngle(beam, q, k1, k2);
        ans = ans - calcTriAngle(beam, q, k1, k2, false);
    } else { // right is higher
        double e = l - a - ((double) 2 / 3) * b;
        ans.first = tmp * (-1 * e * pow(l, 2) + pow(e, 3)
            + (pow(b, 2) / 6) * (c + ((double) 17 / 45) * b)
        );
        ans.second = l * tmp * ( 2 * e * l + pow(e, 2) * ((e / l) - 3)
            + (pow(b, 2) / (6 * l)) * (((double) 17 / 45) * b - d)
        );
    }
    return ans;
}

void calcAngle(Beam& beam, dLoad& df, double k1, double k2) {
    double c = beam.to->x - k2, a = k1 - beam.from->x, b = k2 - k1;
    double q1 = df.val(k1), q2 = df.val(k2);
    pair<double, double> tmpAng{0, 0};
    if(q1 == q2) {
        tmpAng = calcSquareAngle(beam, q1, k1, k2);
    } else {
        if(q1 > q2) {
            tmpAng = calcSquareAngle(beam, q2, k1, k2);
            tmpAng = tmpAng + calcTriAngle(beam, q1 - q2, k1, k2, true);

        } else {
            tmpAng = calcSquareAngle(beam, q1, k1, k2);
            tmpAng = tmpAng + calcTriAngle(beam, q2 - q1, k1, k2, false);
        }
    }
    beam.angle = beam.angle + tmpAng;

    double l = a + b + c;
    double f = b * (q1 + q2) / 2;
    double e = (b *(q1 + 2 * q2)) / (3 * (q1 + q2)) + a;
    beam.from->react += (1 - e / l) * f;
    beam.to->react += (e / l) * f;
}

void calcAngle(Beam& beam, Moment& m) {
    double a = m.x - beam.from->x, b = beam.to->x - m.x, l = a + b;
    beam.from->react -= (m.m / beam.len);
    beam.to->react += (m.m / beam.len);

    double tmp = (m.m * a * b) / (6 * beam.E * beam.I * l);
    pair<double, double> tmpAng(-1 * tmp * (a + 2 * b), tmp * (2 * a + b));
    tmpAng.first = -1 * tmp * (
        6 * a * l - 3 * pow(a, 2) - 2 * pow(l, 2)
    );
    tmpAng.second = tmp * (
        3 * pow(a, 2) - pow(l, 2)
    );
    beam.angle = beam.angle + tmpAng;
}

void calcAngle(Beam& beam, cLoad& f) {
    double a = f.x - beam.from->x, b = beam.to->x - f.x, l = a + b;
    beam.from->react += (b / l) * f.q;
    beam.to->react += (a / l) * f.q;

    double tmp = (f.q * a * b) / (6 * beam.E * beam.I * l);
    pair<double, double> tmpAng(-1 * tmp * (a + 2 * b), tmp * (2 * a + b));
    beam.angle = beam.angle + tmpAng;
}

void Structure::moment3() {
    for(auto f : this->cf) {
        for(auto& beam : this->beam) {
            if(f.x >= beam.from->x && f.x < beam.to->x) {
                calcAngle(beam, f);
            }
        }
    }
    for(auto f : this->df) {
        for(auto& beam : this->beam) {
            double k1 = max(beam.from->x, f.x1), k2 = min(beam.to->x, f.x2);
            double b = k2 - k1;
            if(b > 0) {
                calcAngle(beam, f, k1, k2);
            }
        }
    }

    for(auto m : this->moment) {
        for(auto& beam : this->beam) {
            if(m.x >= beam.from->x && m.x < beam.to->x) {
                calcAngle(beam, m);
            }
        }
    }

    int n = nodeNum;

    // Build M
    Matrix M;
    Vector e1(n); e1[0] = 1;
    M.elem.push_back(e1);
    for (int i = 0; i < beamNum - 1; i++) {
        Vector tmp(n);
        tmp[i] = beam[i].len / (beam[i].I * beam[i].E);
        tmp[i+1] = 2 * (beam[i].len / (beam[i].I * beam[i].E) + beam[i+1].len / (beam[i+1].I * beam[i+1].E));
        tmp[i+2] = beam[i+1].len / (beam[i+1].I * beam[i+1].E);
        M.elem.push_back(tmp);
    }
    Vector en(n); en[n-1] = 1;
    M.elem.push_back(en);

    // Build d
    Vector d(n);
    for (size_t i = 1; i < n - 1; i++) {
        d[i] = 6 * (beam[i].angle.first - beam[i-1].angle.second);
    }

    Vector moment = M.solve(d);

    for (size_t i = 0; i < this->nodeNum; i++) {
        node[i].moment += moment[i];
    }

    for(auto x : this->beam) {
        x.from->react += (x.to->moment - x.from->moment) / x.len;
        x.to->react -= (x.to->moment - x.from->moment) / x.len;
    }

    for(auto x : this->node) {
        cLoad tmp(-1 * x.react, x.x);
        this->reaction.push_back(tmp);
    }
}
