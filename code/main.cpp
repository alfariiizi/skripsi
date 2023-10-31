#include "NumCpp.hpp"
#include "Potential.h"

#include <cstdlib>
#include <iostream>
#include <bits/stdc++.h>
using namespace std;

int main()
{
    time_t start, end;
    time(&start);
    ios_base::sync_with_stdio(false);

    const int nz = 200'000;
    const int nr = 100;
    const int batch = nz / nr;

    const int totalData = 180'000;

    assert(nr * nz == totalData);

    auto rho = nc::fromfile<double>("../build/uji_hasil_gaussian/rho_gaussian_silinder_neumann81.csv", ',');

    rho.reshape(nz, nr);
    Solver solver(nz, nr);
    solver.rho_ndArr = rho;
    solver.setextents(0.025, 0.0, 0.050, 0.035);
    solver.setParam(100'000, 0.01);

    solver.writeSolveGS("phi_gaussian_pred_GS.csv", nr, nr);

    time(&end);
    double time_taken = double(end-start);
    cout << "Time taken by program is: " << fixed << time_taken << setprecision(5);
    cout << "sec" << endl;

    return EXIT_SUCCESS;
}