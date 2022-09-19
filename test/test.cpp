#include <math.h>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <ostream>

#include "const.h"
#include "field.h"
#include "utlis.h"

class Classy {
   private:
    int arraysize;
    // std::unique_ptr<int[]> myarray;
    int *myarray;

   public:
    Classy(int parraysize);
    void printarray();
};

Classy::Classy(int parraysize)
    : arraysize{parraysize}, myarray{new int[arraysize]} {
    for (int i = 0; i < arraysize; i++) {
        myarray[i] = i * i * 2;
    }
}

void Classy::printarray() {
    for (int i = 0; i < arraysize; i++) {
        std::cout << myarray[i] << std::endl;
    }
}

int main(int argc, char *argv[]) { return 0; }
