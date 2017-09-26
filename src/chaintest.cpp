#include <iostream>
#include <scomplex/chains.hpp>

using namespace std;

using namespace gsimp;
int main() {
    vector<double> vec({1,2,3});
    chain my_chain(1,vec);

    cout << my_chain.dimension() << ": ";
    for (auto p : my_chain.get_dense()) cout << p << " ";
    cout << "\n";


    chain other_chain = my_chain+my_chain;

    cout << other_chain.dimension() << ": ";
    for (auto p : other_chain.get_dense()) cout << p << " ";
    cout << "\n";

    cout << (my_chain^my_chain) << "\n";

    my_chain.to_sparse();

}
