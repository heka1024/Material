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
using namespace std;

istream& operator>>(istream& is, Node& n) {
    is >> n.x;
    return is;
}

istream& operator>>(istream& is, Beam& b) {
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

istream& operator>>(istream& is, cForce& cf) {
    is >> cf.q;
    is >> cf.x;
    return is;
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

istream& operator>>(istream& is, Moment& m) {
    is >> m.m;
    is >> m.x;
    return is;
}


ostream& operator<<(ostream& os, cForce& cf) {
    os << cf.q << ", " << cf.x;
    return os;
}

double dForce::val(double x) {
    return this->q1 + (this->q2 - this->q1) * ((x - this->x1) / (this->x2 - this->x1));
}

double dForce::magnitude() {
    const double avgVal = (this->q1 + this->q2) / 2;
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

istream& operator>>(istream& is, dForce& df) {
    is >> df.q1 >> df.q2;
    is >> df.x1 >> df.x2;
    return is;
}

ostream& operator<<(ostream& os, dForce& df) {
    os << df.q1 << ", " << df.q2 << ", ";
    os << df.x1 << ", " << df.x2;
    return os;
}

ostream& operator<<(ostream& os, Mesh& m) {
    os << m.x << ", " << m.sfd << ", " << m.bmd << '\n';
    return os;
}

istream& operator>>(istream& is, Structure& st) {
    cout << "input!\n";
    is >> st.nodeNum >> st.cfNum >> st.dfNum >> st.mNum >> st.totalLen >> st.meshNum ;
    st.beamNum = st.nodeNum - 1;
    cout << "noo\n";


    int idx = 0;
    for(int i = 0; i < st.nodeNum; i++) {
        Node nodeForInput;
        is >> idx >> nodeForInput;
        if(idx != (i + 1)) {
            cout << "error! 1" << '\n';
            exit(1);
        } else {
            st.node.push_back(nodeForInput);
        }
    }
    for(int i = 0; i < st.beamNum; i++) {
        Beam beamForInput;
        is >> idx >> beamForInput;
        if(idx != (i + 1)) {
            cout << "error! 2" << '\n';
            exit(1);
        } else {
            st.beam.push_back(beamForInput);
        }
    }
    for(int i = 0; i < st.cfNum; i++) {
        cForce cfForInput;
        is >> idx >> cfForInput;
        if(idx != (i + 1)) {
            cout << "error! 3" << '\n';
            exit(1);
        } else {
            st.cf.push_back(cfForInput);
        }
    }
    for(int i = 0; i < st.dfNum; i++) {
        dForce dfForInput;
        is >> idx >> dfForInput;
        if(idx != (i + 1)) {
            cout << "error! 4" << '\n';
            exit(1);
        } else {
            st.df.push_back(dfForInput);
        }
    }

    
    for(size_t i = 0; i < st.mNum; i++)
    {
        Moment momentForInput;
        is >> idx >> momentForInput;
        if(idx != (i + 1)) {
            cout << "error!" << '\n';
            exit(1);
        } else {
            st.moment.push_back(momentForInput);
        }

    }

    for (size_t i = 0; i < st.beamNum; i++) {
        st.beam[i].from = &st.node[i];
        st.beam[i].to = &st.node[i+1];
        st.beam[i].len = st.node[i+1].x - st.node[i].x;
    }

    return is;
}

void Structure::calcGraph(vector<Mesh>& mesh) {
    cout << "calculating...!" << '\n';
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
        for (auto x : this->moment) {
            tmpMesh.bmd += x.bodyMoment(pos);
        }
        mesh.push_back(tmpMesh);
        pos += ITER;
    }
    cout << "calc end" << '\n';
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

void calcAngle(Beam& beam, dForce& df, double k1, double k2) {
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

void calcAngle(Beam& beam, cForce& f) {
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
    cout << "----------- Matrix ---------" << '\n';
    cout << M << '\n';
    cout << "-------- angle vector -------" << '\n';
    cout << d << '\n';
    cout << "----------- moment ---------" << '\n';
    Vector moment = M.solve(d);
    cout << moment << '\n';

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
    cout << '\n';
    cout << "Num of Nodes :: " << this->nodeNum << '\n';
    cout << "Num of Beams :: " << this->beamNum << '\n';
    cout << "concentrated Force :: " << this->cfNum << '\n';
    cout << "distributed Force :: " << this->dfNum << '\n';
    cout << "totalLen :: " << this->totalLen << '\n';
    cout << "how many meshes :: " << this->meshNum << "\n\n";

    cout << "----- nodes -----" << '\n';
    for(int i = 0; i < this->nodeNum; i++) {
        printf("%2d %9.3lf\n", i, this->node[i].x);
    }
    cout << "----- beams -----" << '\n';
    for(int i = 0; i < this->beamNum; i++) {
        printf("%2d %9.3lf %9.3lf %9.3lf| ", i, this->beam[i].E, this->beam[i].I, this->beam[i].len);
        printf("from :: %2d %9.3lf | ", i, this->beam[i].from->x);
        printf("to   :: %2d %9.3lf\n", i + 1, this->beam[i].to->x);
    }

    cout << "----- cForce -----" << '\n';
    for(int i = 0; i < this->cfNum; i++) {
        printf("%2d %9.3lf %9.3lf\n", i, this->cf[i].q, this->cf[i].x);
    }

    cout << "----- dForce -----" << '\n';
    for(int i = 0; i < this->dfNum; i++) {
        printf("%2d %9.3lf %9.3lf %9.3lf %9.3lf\n", i, this->df[i].q1, this->df[i].q2, this->df[i].x1, this->df[i].x2);
    }

    cout << "----- moment -----" << '\n';
    for (int i = 0; i < this->mNum; i++) {
        printf("%2d %9.3lf %9.3lf\n", i, this->moment[i].m, this->moment[i].x);
    }

    cout << "----- reaction -----" << '\n';
    for (int i = 0; i < this->reaction.size(); i++) {
        printf("%2d %9.3lf %9.3lf\n", i, this->reaction[i].q, this->reaction[i].x);
    }
}
