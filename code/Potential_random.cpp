#include "Potential.h"
#include "Field.h"

#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

Solver::Solver(int nz, int nr)
    : nz{nz},
      nr{nr},
      phi(nz, nr),
      rho(nz, nr),
      object_id(nz, nr)
{
}

void Solver::setextents(double zmax, double zmin, double rmax, double rmin)
{
    this->x0[0] = zmin;
    this->x0[1] = rmin;

    this->xm[0] = zmax;
    this->xm[1] = rmax;

    this->dh[0] = (zmax - zmin) / 100;
    this->dh[1] = (rmax - rmin) / 100;

    this->L[0] = zmax - zmin;
    this->L[1] = rmax - rmin;
}

void Solver::setParam(int iter, double tol)
{
    this->max_solver_it = iter;
    this->tolerance = tol;
}


void Solver::writeSolveGS(const std::string& filename, int reshape_rows, int reshape_cols) {
    std::ofstream out(filename);

    double idz = 1 / this->dh[0];
    double idr = 1 / this->dh[1];
    double idz2 = idz * idz;
    double idr2 = idr * idr;

    double L2 = 0;
    bool converged = false;

    double crz = 0.5 / (idz2 + idr2);
    for (unsigned it = 0; it < this->max_solver_it; it++)
    {
        for (int i = 0; i < this->nz; i++)
        {
            for (int j = 0; j < this->nr; j++)
            {
                if (i == 0 || i % reshape_rows == 0) 
                {
                    continue;
                }
                else if (i % reshape_rows == 99) 
                {

                    phi[i][j] = phi[i-1][j];
                }
                else if (j == 0) 
                {
                    continue; 

                }
                else if (j == (this->nr) - 1) 
                {
                    continue; 
                }
                else // selain di syarat batas
                {
                    double crj = 0.5 / ((this->x0[1] + j * this->dh[1]) * this->dh[1]);
                    double phi_baru = crz * ((this->rho_ndArr(i, j) / Const::EPS_0) + (idz2 * (this->phi[i + 1][j] + this->phi[i - 1][j])) + (this->phi[i][j + 1] * (idr2 + crj)) + (this->phi[i][j - 1] * (idr2 - crj)));
                    // lanjutkan dengan SOR
                    phi[i][j] += 1.4 * (phi_baru - this->phi[i][j]);
                }

            }
        }
        if (it % 100 == 0)
        {
            double sum = 0;
            for (int i = 1; i < reshape_rows - 1; i++)
            {
                for (int j = 1; j < reshape_cols - 1; j++)
                {
                    double crj = 0.5 / ((this->x0[1] + j * this->dh[1]) * this->dh[1]);
                    double R = -2 * this->phi[i][j] * (idr2 + idz2) + ((this->rho_ndArr(i, j) / Const::EPS_0)
                            + idz2 * (this->phi[i + 1][j]
                            + this->phi[i - 1][j]) + this->phi[i][j + 1] * (idr2 + crj)
                            + this->phi[i][j - 1] * (idr2 - crj));
                    sum += R * R;
                }
            }
            L2 = sqrt(sum / (reshape_cols * reshape_rows));
            if (L2 < tolerance)
            {
                converged = true;
                break;
            }
        }
    }
    if (!converged)
    {
        cerr << "Gauss seidel standar gagal konvergen, L2=" << L2 << endl;
    }
    for(int i=0; i<nz; i++)
    {
        for(int j=0; j<nr; j++)
        {
            if(j==this->nr-1)
            {
                out<<phi[i][j]<<'\n';
            }
            else
            {
                out<<phi[i][j]<<',';
            }
        }
    }
    out.close();
}