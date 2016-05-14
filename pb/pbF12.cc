#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>

#include "F1F209Wrapper.hh"

using namespace std;

int main(int argc, char *argv[])
{
    int in = argc > 1 ? atoi(argv[1]) : 1;
    int Z = argc > 2 ? atoi(argv[2]) : 1;
    int A = argc > 3 ? atoi(argv[3]) : 1;
    double F1, F2;
    double R = 0;
    
    char infname[50], outfname[50];
    sprintf(infname, "WQ2list%i%i%s.dat", Z, A, argc > 4 ? argv[4] : "");
    sprintf(outfname, "PBF12%i%i%s.dat", Z, A, argc > 4 ? argv[4] : "");

    ifstream infile(infname);
    ofstream outfile(outfname);

    double W2, Q2;
    F1F209Wrapper *model = new F1F209Wrapper();
    int i = 0;

    while (infile >> W2 >> Q2) {
        i++;
        if (in) {
            model->GetF1F2IN09(Z, A, Q2, W2, F1, F2, R);
        }
        else {
            model->GetF1F2QE09(Z, A, Q2, W2, F1, F2);
        }
        if (i % 10000 == 0) {
            printf("%i %i %i %f %f %f %f\n", i, Z, A, Q2, W2, F1, F2);
        }
        outfile << F1 << " " << F2 << " " << R << "\n";
    }

    infile.close();
    outfile.close();

    return 0;
}
