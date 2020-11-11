#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "Material.h"
#include "Matrix.h"
#include "Vector.h"
std::string inputName = "../input/input.txt";

int main() {
    Structure s;
    std::ifstream ifs(inputName);
    ifs >> s; ifs.close();
    s.moment3();
    std::vector<Mesh> meshes;
    s.calcGraph(meshes);
    s.ls();
    std::ofstream ofs("../data.txt");
    for(auto x : meshes) {
        ofs << x;
    }
    ofs.close();


    return 0;
}
