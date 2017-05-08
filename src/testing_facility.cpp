#include <yaml-cpp/yaml.h>
#include <iostream>

int main( int argc , char* argv[] ) {

    YAML::Node input = YAML::LoadFile(argv[1]);

    std::string paths_file;
    std::string complex_file;

    paths_file = input[0]["paths_file"].as<std::string>();
    complex_file = input[0]["complex_file"].as<std::string>();

    std::cout << paths_file << "\n";
    std::cout << complex_file << "\n";
}
