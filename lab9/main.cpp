#include <iostream>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

double d2t(gsl_vector* T, int l, int nx, double delta) {
    return (( gsl_vector_get(T, l+1) - 2 * gsl_vector_get(T, l) + gsl_vector_get(T, l-1)) / pow(delta, 2))
            + (( gsl_vector_get(T, l+nx+1) - 2 * gsl_vector_get(T, l) + gsl_vector_get(T, l-nx-1)) / pow(delta, 2));
}

void thermoDiffusion() {
    
    //1
    int nx = 40;
    int ny = 40;
    int N = (nx + 1) * (ny + 1);
    double delta = 1;
    double deltaT = 1;
    int tA = 40;
    int tB = 0;
    int tC = 30;
    int tD = 0;
    double kB = 0.1;
    double kD = 0.6;
    int IT_MAX = 2000;

    //2
    gsl_matrix* A = gsl_matrix_calloc(N, N);
    gsl_matrix* B = gsl_matrix_calloc(N, N);
    gsl_vector* C = gsl_vector_calloc(N);

    int l;
    for(int i=1; i<nx; i++) {
        for(int j=1; j<ny; j++) {
            l = i + j * (nx + 1);
            gsl_matrix_set(A, l, l-nx-1, deltaT / (2*pow(delta, 2)));
            gsl_matrix_set(A, l, l-1, deltaT / (2*pow(delta, 2)));
            gsl_matrix_set(A, l, l+1, deltaT / (2*pow(delta, 2)));
            gsl_matrix_set(A, l, l+nx+1, deltaT / (2*pow(delta, 2)));
            gsl_matrix_set(A, l, l, -(2*deltaT)/pow(delta,2) - 1);
            
            gsl_matrix_set(B, l, l-nx-1, -deltaT / (2*pow(delta, 2)));
            gsl_matrix_set(B, l, l-1, -deltaT / (2*pow(delta, 2)));
            gsl_matrix_set(B, l, l+1, -deltaT / (2*pow(delta, 2)));
            gsl_matrix_set(B, l, l+nx+1, -deltaT / (2*pow(delta, 2)));
            gsl_matrix_set(B, l, l, (2*deltaT)/pow(delta,2) - 1);
        }
    }

    for(int i=0; i<=nx; i+=nx) {
        for(int j=0; j<=ny; j++) {
            l = i + j * (nx + 1);
            gsl_matrix_set(A, l, l, 1); 
            gsl_matrix_set(B, l, l, 1);
            gsl_vector_set(C, l, 0);
        }
    }

    for(int i=1; i<nx; i++) { 
        l = i + ny * (nx + 1);
        gsl_matrix_set(A, l, l-nx-1, -1.0 / (kB*delta)); 
        gsl_matrix_set(A, l, l, 1 + 1.0 / (kB*delta)); 
        gsl_vector_set(C, l, tB);
        for(int k=0; k<N; k++) {
            gsl_matrix_set(B, l, k, 0);
        }
    }

    for(int i=1; i<nx; i++) {
        l = i + 0 * (nx + 1);
        gsl_matrix_set(A, l, l+nx+1, -1.0 / (kD*delta)); 
        gsl_matrix_set(A, l, l, 1 + 1.0 / (kD*delta)); 
        gsl_vector_set(C, l, tD);
        for(int k=0; k<N; k++) {
            gsl_matrix_set(B, l, k, 0);
        }
    }

    //3
    gsl_vector* T = gsl_vector_calloc(N);

    for(int j=0; j<=ny; j++) { 
        l = 0 + j * (nx + 1);
        gsl_vector_set(T, l, tA); 
        l = nx + j * (nx + 1);
        gsl_vector_set(T, l, tC); 
    }

    for(int i=1; i<nx; i++) {
        for(int j=0; j<=ny; j++) {
            l = i + j * (nx + 1);
            gsl_vector_set(T, l, 0);
        }
    }

    //4
    gsl_permutation* p = gsl_permutation_calloc(N);
    int signum = 0;
    
    gsl_linalg_LU_decomp(A, p, &signum);

    //5, 6
    gsl_vector *d = gsl_vector_calloc(N);

    for(int it=0; it<=IT_MAX; it++) {
        //d = B*T
        gsl_blas_dgemv(CblasNoTrans, 1, B, T, 0, d);
        //d += C
        gsl_blas_daxpy(1, C, d);

        gsl_linalg_LU_solve(A, p, d, T);

        //7, 8
        if(it == 100 || it == 200 || it == 500 || it == 1000 || it == 2000) {
            FILE* f1 = fopen(("T(x,y)_" + std::to_string(it) + ".txt").c_str(), "w");
            FILE* f2 = fopen(("d2T(x,y)_" + std::to_string(it) + ".txt").c_str(), "w");
            for(int i=1; i<nx; i++) {
                for( int j=1; j<ny; j++) {
                    fprintf(f1, "%g %g %g\n", i*delta, j*delta, gsl_vector_get(T, i + j * (nx + 1)));
                    fprintf(f2, "%g %g %g\n", i*delta, j*delta, d2t(T, i + j * (nx + 1), nx, delta));
                }
                fprintf(f1, "\n");
                fprintf(f2, "\n");
            }
            fclose(f1);
            fclose(f2);
        }
    }
}

int main() {

    thermoDiffusion();

    return 0;
}