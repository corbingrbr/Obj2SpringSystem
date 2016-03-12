#include <iostream>
#include "ObjFrame.h"
#include <cstdlib>
#include <assert.h>

using namespace std;

int main(int argc, char **argv) {

    if (argc == 3) {
       
        ObjFrame objf(atoi(argv[2]));
        objf.framerize(argv[1]);
        
        cout << objf.f2Ndx(.5) << endl;
        




    } else {
        cout << "specify .obj and res" << endl;
    }

    return 0;
}
