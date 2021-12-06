#include <iostream>
#include <fstream>
#include <math.h>
#include "mgmres.h"

void poisson(std::string name, int eps1, int eps2, int nx, int ny, int V1, int V2, int V3, int V4, double p1, double p2){
    
    double delta = 0.1;
    int N = (nx + 1) * (ny + 1);
    int k = -1;
    double* a = new double[5*N];
    int* ja = new int[5*N];
    int* ia = new int[N+1];
    double* b = new double[N];
    double* epsl = new double[N];
    double* V = new double[N];
    int nz_num;
    int itr_max = 500;
    int mr = 500;
    double tol_abs = 1e-8;
    double tol_rel = 1e-8;
    int i, j;

    double sigma = (delta*nx)/10;

    for(int i=0; i<N; i++){
        V[i] = 0;
        b[i] = 0;
        ia[i] = -1;
    }
    ia[N] = -1;

    for(int l=0; l<N; l++){
        j = l/(nx+1);
        i = l - j * (nx + 1);

        if (i <= nx/2)
            epsl[l] = eps1;
        else
            epsl[l] = eps2;
    }

    for(int l=0; l<N; l++){

        int brzeg = 0;
        double vb = 0;
        j = l/(nx+1);
        i = l - j * (nx + 1);

        if(i == 0){
            brzeg = 1;
            vb = V1;
        }

        if(j == ny){
            brzeg = 1;
            vb = V2;
        }

        if(i == nx){
            brzeg = 1;
            vb = V3;
        }

        if(j == 0){
            brzeg = 1;
            vb = V4;
        }

        if(V1 == 0 && V2 == 0 && V3 == 0 && V4 == 0){
            b[l] = - (exp( -(pow(delta*i - 0.25*delta*nx, 2) / pow(sigma, 2)) - (pow(delta*j - 0.5*delta*ny, 2) / pow(sigma, 2))) +
          (-exp(-(pow(delta*i - 0.75*delta*nx, 2)/pow(sigma, 2)) -(pow(delta*j - 0.5*delta*ny, 2) / pow(sigma, 2)))));
        } else {
            b[l] = - (p1 + p2);
        }

        if (brzeg == 1){
            b[l] = vb;
        }

        ia[l] = -1;

        if(l-nx-1 >= 0 && brzeg == 0){
            k++;
            if(ia[l] < 0)
                ia[l] = k;
            a[k] = epsl[l]/(delta*delta);
            ja[k] = l - nx - 1;
        }

        if(l-1 >= 0 && brzeg == 0){
            k++;
            if(ia[l] < 0)
                ia[l] = k;
            a[k] = epsl[l]/(delta*delta);
            ja[k] = l-1;
        }

        k++;
        if(ia[l] < 0)
            ia[l] = k;
        
        if(brzeg == 0){
            a[k] = -(2*epsl[l] + epsl[l+1] + epsl[l+nx+1]) / (delta*delta);
        } else {
            a[k] = 1;
        }
        ja[k] = l;

        if(l<N && brzeg == 0){
            k++;
            a[k] = epsl[l+1]/(delta*delta);
            ja[k] = l+1;
        }

        if(l < N-nx-1 && brzeg == 0){
            k++;
            a[k] = epsl[l+nx+1] / (delta*delta);
            ja[k] = l + nx + 1;
        } 
    }

    nz_num = k+1;
    ia[N] = nz_num;

    pmgmres_ilu_cr(N, nz_num, ia, ja, a, V, b, itr_max, mr, tol_abs, tol_rel);

    if(nx == 4 && ny == 4){
        std::fstream file;
        file.open(name, std::ios::out);
        if(file.good()){
            for(int l=0; l<5*N; l++){
                file << l << "\t" << a[l] << std::endl;
            }
            file.close();
        }

        std::fstream file2;
        file2.open("matrix_b.txt", std::ios::out);
        if(file2.good()){
            for(int l=0; l<N; l++){
                file2 << l << "\t" << b[l] << std::endl;
            }
            file2.close();
        }
    } else {
        std::fstream file;
        file.open(name, std::ios::out);
        if(file.good()){
            for(int l=0; l<N; l++){
                double j = l/(nx+1);
                double i = l - j * (nx + 1);
                file << i << " " << j << " " << V[l] << std::endl;
            }
            file.close();
        }
    }

    delete [] a;
    delete [] ia;
    delete [] ja;
    delete [] V;
    delete [] b;
    delete [] epsl;
}

int main(){

    poisson("matrix_a.txt", 1, 1, 4, 4, 10, -10, 10, -10, 0, 0);
    poisson("map1_a.txt", 1, 1, 50, 50, 10, -10, 10, -10, 0, 0);
    poisson("map1_b.txt", 1, 1, 100, 100, 10, -10, 10, -10, 0, 0);
    poisson("map1_c.txt", 1, 1, 200, 200, 10, -10, 10, -10, 0, 0);
    poisson("map2_a.txt", 1, 1, 100, 100, 0, 0, 0, 0, 0, 0);
    poisson("map2_b.txt", 1, 2, 100, 100, 0, 0, 0, 0, 0, 0);
    poisson("map2_c.txt", 1, 10, 100, 100, 0, 0, 0, 0, 0, 0);

    return 0;
}