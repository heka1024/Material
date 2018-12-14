#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "Material.h"
#include "iMath.h"

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

double dForce::magnitude() {
    double avgVal = (this->q1 + this->q2) / 2;
    return avgVal * (this->x2 - this->x1);
}

double dForce::shearForce(double x) {
    if(x < this->x1) {
        return 0;
    } else if (x < this->x2) {
        double qx = (this->q1 + ((this->q2 - this->q1) * ((x - this->x1) / (this->x2 - this->x1))));
        return (x - this->x1) * (this->q1 + qx) / 2;
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

    return is;
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
        printf("%2d %9.3lf %9.3lf\n", i, this->beam[i].E, this->beam[i].I);
    }

    std::cout << "----- cForce -----" << '\n';
    for(int i = 0; i < this->cfNum; i++) {
        printf("%2d %9.3lf %9.3lf\n", i, this->cf[i].q, this->cf[i].x);
    }

    std::cout << "----- dForce -----" << '\n';
    for(int i = 0; i < this->dfNum; i++) {
        printf("%2d %9.3lf %9.3lf %9.3lf %9.3lf\n", i, this->df[i].q1, this->df[i].q2, this->df[i].x1, this->df[i].x2);
    }
}
