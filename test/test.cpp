#include <Eigen/Core>
#include <iostream>
using namespace Eigen;
using namespace std;

int main(int argc, char* argv[]) {
//    ArrayXXd a(2, 2);
//    a << 1, 2, 3, 4;
//    ArrayXXd b(2, 2);
//    b << 1, 2, 3, 4;

    int a = 1;
    while (true) {
        int b = 10;
        a++;
        std::cout << a << std::endl;
        if (a == b) break;
    }

//    if (a.isApprox(b))
//        cout << "equal" << endl;
//    else
//        cout << "not equal" << endl;

    return 0;
}