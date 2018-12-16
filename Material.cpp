#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <string>
#include <utility>
#include <algorithm>
#include "Material.h"
#include "Matrix.h"
#include "Vector.h"

std::istream& operator>>(std::istream& is, Node& n) {
    is >> n.x;
    return is;
}

std::istream& operator>>(std::istream& is, Beam& b) {
    is >> b.E >> b.I;
    return is;
}

double cForce::magnitude() {
    return this->q;
}

double cForce::u(double x) {
    return (x >= this->x);
}

double cForce::shearForce(double x) {
    return this->u(x) * this->magnitude();
}

std::istream& operator>>(std::istream& is, cForce& cf) {
    is >> cf.q;
    is >> cf.x;
    return is;
}

std::ostream& operator<<(std::ostream& os, cForce& cf) {
    os << cf.q << ", " << cf.x;
}

double dForce::val(double x) {
    return this->q1 + (this->q2 - this->q1) * ((x - this->x1) / (this->x2 - this->x1));
}

double dForce::magnitude() {
    double avgVal = (this->q1 + this->q2) / 2;
    return avgVal * (this->x2 - this->x1);
}

double dForce::shearForce(double x) {
    if(x < this->x1) {
        return 0;
    } else if (x < this->x2) {
        return (x - this->x1) * (this->val(x) + q1) / 2;
    } else {
        return this->magnitude();
    }
}

std::istream& operator>>(std::istream& is, dForce& df) {
    is >> df.q1 >> df.q2;
    is >> df.x1 >> df.x2;
    return is;
}

std::ostream& operator<<(std::ostream& os, dForce& df) {
    os << df.q1 << ", " << df.q2 << ", ";
    os << df.x1 << ", " << df.x2;
}

std::ostream& operator<<(std::ostream& os, Mesh& m) {
    os << m.x << ", " << m.sfd << ", " << m.bmd << '\n';
}

std::istream& operator>>(std::istream& is, Structure& st) {
    is >> st.nodeNum >> st.cfNum >> st.dfNum >> st.totalLen >> st.meshNum;
    st.beamNum = st.nodeNum - 1;

    int idx = 0;
    for(int i = 0; i < st.nodeNum; i++) {
        Node nodeForInput;
        is >> idx >> nodeForInput;
        if(idx != (i + 1)) {
            std::cout << "error!" << '\n';
            exit(1);
        } else {
            st.node.push_back(nodeForInput);
        }
    }
    for(int i = 0; i < st.beamNum; i++) {
        Beam beamForInput;
        is >> idx >> beamForInput;
        if(idx != (i + 1)) {
            std::cout << "error!" << '\n';
            exit(1);
        } else {
            st.beam.push_back(beamForInput);
        }
    }
    for(int i = 0; i < st.cfNum; i++) {
        cForce cfForInput;
        is >> idx >> cfForInput;
        if(idx != (i + 1)) {
            std::cout << "error!" << '\n';
            exit(1);
        } else {
            st.cf.push_back(cfForInput);
        }
    }
    for(int i = 0; i < st.dfNum; i++) {
        dForce dfForInput;
        is >> idx >> dfForInput;if(idx != (i + 1)) {
            std::cout << "error!" << '\n';
            exit(1);
        } else {
            st.df.push_back(dfForInput);
        }
    }
    for (size_t i = 0; i < st.beamNum; i++) {
        st.beam[i].from = &st.node[i];
        st.beam[i].to = &st.node[i+1];
        st.beam[i].len = st.node[i+1].x - st.node[i].x;
    }

    return is;
}

void Structure::calcGraph(std::vector<Mesh>& mesh) {
    std::cout << "calculating...!" << '\n';
    double ITER = ((double) totalLen / meshNum);
    double pos = 0;
    for (size_t i = 0; i < meshNum; i++) {
        Mesh tmpMesh;
        tmpMesh.x = pos;
        for(auto x : this->cf) {
            tmpMesh.sfd += x.shearForce(pos);
        }
        for(auto x : this->df) {
            tmpMesh.sfd += x.shearForce(pos);
        }
        for(auto x : this->reaction) {
            tmpMesh.sfd += x.shearForce(pos);
        }
        tmpMesh.sfd *= -1;
        for (auto x : mesh) {
            tmpMesh.bmd += x.sfd * ITER;
        }
        mesh.push_back(tmpMesh);
        pos += ITER;
    }
    std::cout << "calc end" << '\n';
}

template <typename T,typename U>
std::pair<T,U> operator+(const std::pair<T,U> & l,const std::pair<T,U> & r) {
    return {l.first + r.first, l.second + r.second};
}

template <typename T,typename U>
std::pair<T,U> operator-(const std::pair<T,U> & l,const std::pair<T,U> & r) {
    return {l.first - r.first, l.second - r.second};
}

std::pair<double, double> calcSquareAngle(Beam &beam, double q, double k1, double k2) {
    double c = beam.to->x - k2, a = k1 - beam.from->x, b = k2 - k1;
    double l = beam.to->x - beam.from->x, d = a + b;
    std::pair<double, double> ans{0, 0};

    double tmp = (double) (q * b) / (48 * beam.E * beam.I * l);
    ans.first -= tmp * (b + 2 * c) * (
        4 * std::pow(l, 2) - std::pow(b, 2) - std::pow(b + 2 * c, 2)
    );
    ans.second += tmp * (2 * d -b) * (
        2 * b * (a + d) + 4 * c * (l + d)
    );
    return ans;
}

std::pair<double, double> calcTriAngle(Beam &beam, double q, double k1, double k2, bool leftHigher) {
    double c = beam.to->x - k2, a = k1 - beam.from->x, b = k2 - k1;
    double l = beam.to->x - beam.from->x, d = a + b;
    double tmp = (q * b) / (12 * beam.E * beam.I * l);
    std::pair<double, double> ans(0, 0);
    if(leftHigher) {
        ans = calcSquareAngle(beam, q, k1, k2);
        ans = ans - calcTriAngle(beam, q, k1, k2, false);
    } else { // right is higher
        double e = l - a - ((double) 2 / 3) * b;
        ans.first = tmp * (-1 * e * std::pow(l, 2) + std::pow(e, 3)
            + (std::pow(b, 2) / 6) * (c + ((double) 17 / 45) * b)
        );
        ans.second = l * tmp * ( 2 * e * l + std::pow(e, 2) * ((e / l) - 3)
            + (std::pow(b, 2) / (6 * l)) * (((double) 17 / 45) * b - d)
        );
    }
    return ans;
}

void calcAngle(Beam& beam, dForce& df, double k1, double k2) {
    double c = beam.to->x - k2, a = k1 - beam.from->x, b = k2 - k1;
    double q1 = df.val(k1), q2 = df.val(k2);
    std::pair<double, double> tmpAng{0, 0};
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

void calcAngle(Beam& beam, cForce& f) {
    double a = f.x - beam.from->x, b = beam.to->x - f.x, l = a + b;
    beam.from->react += (b / l) * f.q;
    beam.to->react += (a / l) * f.q;

    double tmp = (f.q * a * b) / (6 * beam.E * beam.I * l);
    std::pair<double, double> tmpAng(-1 * tmp * (a + 2 * b), tmp * (2 * a + b));
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
            double k1 = std::max(beam.from->x, f.x1), k2 = std::min(beam.to->x, f.x2);
            double b = k2 - k1;
            if(b > 0) {
                calcAngle(beam, f, k1, k2);
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
    M.rows = n; M.cols = n;

    // Build d
    Vector d(n);
    for (size_t i = 1; i < n - 1; i++) {
        d[i] = 6 * (beam[i].angle.first - beam[i-1].angle.second);
    }
    std::cout << "----------- Matrix ---------" << '\n';
    std::cout << M << '\n';
    std::cout << "-------- angle vector -------" << '\n';
    std::cout << d << '\n';
    std::cout << "----------- moment ---------" << '\n';
    Vector moment = M.solve(d);
    std::cout << moment << '\n';

    for (size_t i = 0; i < this->nodeNum; i++) {
        node[i].moment += moment[i];
    }

    for(auto x : this->beam) {
        x.from->react += (x.to->moment - x.from->moment) / x.len;
        x.to->react -= (x.to->moment - x.from->moment) / x.len;
    }

    for(auto x : this->node) {
        cForce tmp(-1 * x.react, x.x);
        this->reaction.push_back(tmp);
    }
}

void Structure::ls() {
    std::cout << '\n';
    std::cout << "Num of Nodes :: " << this->nodeNum << '\n';
    std::cout << "Num of Beams :: " << this->beamNum << '\n';
    std::cout << "concentrated Force :: " << this->cfNum << '\n';
    std::cout << "distributed Force :: " << this->dfNum << '\n';
    std::cout << "totalLen :: " << this->totalLen << '\n';
    std::cout << "how many meshes :: " << this->meshNum << "\n\n";

    std::cout << "----- nodes -----" << '\n';
    for(int i = 0; i < this->nodeNum; i++) {
        printf("%2d %9.3lf\n", i, this->node[i].x);
    }
    std::cout << "----- beams -----" << '\n';
    for(int i = 0; i < this->beamNum; i++) {
        printf("%2d %9.3lf %9.3lf %9.3lf| ", i, this->beam[i].E, this->beam[i].I, this->beam[i].len);
        printf("from :: %2d %9.3lf | ", i, this->beam[i].from->x);
        printf("to   :: %2d %9.3lf\n", i + 1, this->beam[i].to->x);
    }

    std::cout << "----- cForce -----" << '\n';
    for(int i = 0; i < this->cfNum; i++) {
        printf("%2d %9.3lf %9.3lf\n", i, this->cf[i].q, this->cf[i].x);
    }

    std::cout << "----- dForce -----" << '\n';
    for(int i = 0; i < this->dfNum; i++) {
        printf("%2d %9.3lf %9.3lf %9.3lf %9.3lf\n", i, this->df[i].q1, this->df[i].q2, this->df[i].x1, this->df[i].x2);
    }

    std::cout << "----- reaction -----" << '\n';
    for (int i = 0; i < this->reaction.size(); i++) {
        printf("%2d %9.3lf %9.3lf\n", i, this->reaction[i].q, this->reaction[i].x);
    }
}
