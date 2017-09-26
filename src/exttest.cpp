#include <iostream>
#include <scomplex/exterior_prod.hpp>

using namespace std;

int main() {

    ext_vector a({1,0,0});
    ext_vector b({0,1,0});

    ext_vector c = b^a;

    cout << a.to_str();
    cout << b.to_str();
    cout << c.to_str();
    cout << c.norm();

    return 0;
}
