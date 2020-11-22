#include <iostream>
#include <string>
#include <vector>
#include "Material.h"
using namespace std;
std::string inputName = "../input/input.txt";

void print_structure(const Structure& structure) {
    std::cout << '\n';
    std::cout << "Num of Nodes :: " << structure.nodeNum << '\n';
    std::cout << "Num of Beams :: " << structure.beamNum << '\n';
    std::cout << "concentrated Force :: " << structure.cfNum << '\n';
    std::cout << "distributed Force :: " << structure.dfNum << '\n';
    std::cout << "totalLen :: " << structure.totalLen << '\n';
    std::cout << "how many meshes :: " << structure.meshNum << "\n\n";

    std::cout << "----- nodes -----" << '\n';
    for(int i = 0; i < structure.nodeNum; i++) {
        printf("%2d %9.3lf\n", i, structure.node[i].x);
    }
    std::cout << "----- beams -----" << '\n';
    for(int i = 0; i < structure.beamNum; i++) {
        printf("%2d %9.3lf %9.3lf %9.3lf| ", i, structure.beam[i].E, structure.beam[i].I, structure.beam[i].len);
        printf("from :: %2d %9.3lf | ", i, structure.beam[i].from->x);
        printf("to   :: %2d %9.3lf\n", i + 1, structure.beam[i].to->x);
    }

    std::cout << "----- cLoad -----" << '\n';
    for(int i = 0; i < structure.cfNum; i++) {
        printf("%2d %9.3lf %9.3lf\n", i, structure.cf[i].q, structure.cf[i].x);
    }

    std::cout << "----- dLoad -----" << '\n';
    for(int i = 0; i < structure.dfNum; i++) {
        printf("%2d %9.3lf %9.3lf %9.3lf %9.3lf\n", i, structure.df[i].q1, structure.df[i].q2, structure.df[i].x1, structure.df[i].x2);
    }

    std::cout << "----- moment -----" << '\n';
    for (int i = 0; i < structure.mNum; i++) {
        printf("%2d %9.3lf %9.3lf\n", i, structure.moment[i].m, structure.moment[i].x);
    }

    std::cout << "----- reaction -----" << '\n';
    for (int i = 0; i < structure.reaction.size(); i++) {
        printf("%2d %9.3lf %9.3lf\n", i, structure.reaction[i].q, structure.reaction[i].x);
    }
}

vector<Mesh> structure2mesh(const Structure& strt) {
    double ITER = ((double) strt.totalLen / strt.meshNum);
    double pos = 0;
    vector<Mesh>* ans = new vector<Mesh>;
    for (size_t i = 0; i < strt.meshNum; i++) {
        Mesh tmpMesh;
        tmpMesh.x = pos;
        for(auto x : strt.cf) {
            tmpMesh.sfd += x.shearForce(pos);
        }
        for(auto x : strt.df) {
            tmpMesh.sfd += x.shearForce(pos);
        }
        for(auto x : strt.reaction) {
            tmpMesh.sfd += x.shearForce(pos);
        }
        tmpMesh.sfd *= -1;
        for (auto x : *ans) {
            tmpMesh.bmd += x.sfd * ITER;
        }
        for (auto x : strt.moment) {
            tmpMesh.bmd += x.bodyMoment(pos);
        }
        ans->push_back(tmpMesh);
        pos += ITER;
    }
    std::cout << "calc end" << '\n';
    return *ans;
}

int main() {
    Structure st;

    FILE* fin;
    fin = fopen("../input/input.txt", "r");
    fscanf(fin, "%d %d %d %d %lf %d", &st.nodeNum, &st.cfNum, &st.dfNum, &st.mNum, &st.totalLen, &st.meshNum);
    st.beamNum = st.nodeNum - 1;

    int idx = 0; int tmp;
    for(int i = 0; i < st.nodeNum; i++) {
        Node nodeForInput; 
        fscanf(fin, "%d %lf", &tmp, &nodeForInput.x);
        st.node.push_back(nodeForInput);
    }
    for(int i = 0; i < st.beamNum; i++) {
        Beam beamForInput;
        fscanf(fin, "%d %lf %lf", &tmp, &beamForInput.E, &beamForInput.I);
        st.beam.push_back(beamForInput);
    }
    for(int i = 0; i < st.cfNum; i++) {
        cLoad cfForInput;
        fscanf(fin, "%d %lf %lf", &tmp, &cfForInput.q, &cfForInput.x);
        st.cf.push_back(cfForInput);
    }
    for(int i = 0; i < st.dfNum; i++) {
        dLoad dfForInput;
        fscanf(fin, "%d %lf %lf %lf %lf", &idx, &dfForInput.q1, &dfForInput.q2, &dfForInput.x1, &dfForInput.x2);
        st.df.push_back(dfForInput);
    }
    for(size_t i = 0; i < st.mNum; i++)
    {
        Moment momentForInput;
        fscanf(fin, "%d %lf %lf", &idx, &momentForInput.m, &momentForInput.x);
        st.moment.push_back(momentForInput);
    }
    for (size_t i = 0; i < st.beamNum; i++) {
        st.beam[i].from = &st.node[i];
        st.beam[i].to = &st.node[i+1];
        st.beam[i].len = st.node[i+1].x - st.node[i].x;
    }    
    fclose(fin);
    print_structure(st);
    
    st.moment3();
    const auto& meshes = structure2mesh(st);

    FILE* fp;
    fp = fopen("../output/data.txt", "w");
    for (const auto& mesh: meshes) {
        fprintf(fp, "%lf %lf %lf\n", mesh.x, mesh.sfd, mesh.bmd);
    }
    fclose(fp);

    return 0;
}
