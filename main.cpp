#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include "Material.h"
#include "Matrix.h"
#include "Vector.h"
#include <emscripten/bind.h>

using namespace emscripten;

std::string toString( std::ostream& str ) {
    std::ostringstream ss;
    ss << str.rdbuf();
    return ss.str();
}


std::string inputName = "../input/input.txt";
std::string input = "3 1 1 1 20 10000\
\n1 0\
\n2 8\
\n3 20\
\n1 10000 1\
\n2 10000 1\
\n1 50000 4\
\n1 3000 3000 8 20\
\n1 3000 2";

std::string calculator(std::string input_string) {
    std::istringstream str(input_string);
    Structure s;
    str >> s;
    std::cout << "input succeess!\n";
    s.moment3();
    std::vector<Mesh> meshes;
    s.calcGraph(meshes);
    s.ls();

    std::ostringstream os;
    for(auto x : meshes) {
        os << x;
    }
    return os.str();
}

EMSCRIPTEN_BINDINGS(my_module) {
    function("calculator", &calculator);
}

int main() {
    std::cout << calculator(input);
    // std::istringstream str(input);     


    // Structure s;
    // // std::ifstream ifs(inputName);
    // // ifs >> s; ifs.close();
    // str >> s;
    // std::cout << "input succeess!\n";
    // s.moment3();
    // std::vector<Mesh> meshes;
    // s.calcGraph(meshes);
    // s.ls();
    // // std::ofstream ofs("../data.txt");
    // // for(auto x : meshes) {
    // //     ofs << x;
    // // }
    // // ofs.close();
    // ostringstream os;
    // for(auto x : meshes) {
    //     os << x;
    // }
    // cout << os.str();


    return 0;
}
