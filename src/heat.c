//heat.c
/*
Autor: Emmanuel Alejandro Larralde Ortiz | emmanuel.larralde@cimat.mx
Descripción:
    Resuelve la ecuación de calor en una barra unidimensional.
Compilación:
    gcc src/heat.c include/matrices/matrices.c -o output/matrices.o
    chmod +x output/matrices.o
Uso:
    ./output/heat.o <Q> <K> <phi_0> <phi_n> <L> <n>
*/

#include <stdlib.h>
#include <stdio.h>
#include "../include/matrices/matrices.h"


int main(int argc, char **argv){
    if(argc < 7){
        printf("Error. Incorrect usage\n");
        printf("./output/heat.o <Q> <K> <phi_0> <phi_n> <L> <n>\n");
        return 0;
    }

    double Q, K, phi_0, phi_n, L;
    int n;
    Q = atof(argv[1]);
    K = atof(argv[2]);
    phi_0 = atof(argv[3]);
    phi_n = atof(argv[4]);
    L = atof(argv[5]);
    n = atoi(argv[6]);

    int nxs = n-1; //Number of unknowns.
    double dx, bi;
    dx = L/((double) n); //size of the finite element.
    bi = (Q * dx*dx)/K; //Qdx^2/K

    double *mat, *b;
    //Knowns vector init.
    b = (double *) malloc(nxs * sizeof(double));
    for(int i=0; i<nxs; ++i)
        *(b + i) = bi;
    *b += phi_0;
    *(b + nxs - 1) += phi_n;

    //Matrix initialization.
    mat = (double *) malloc(nxs * nxs * sizeof(double));
    for(int i=0; i<nxs; ++i)
        *(mat + i*nxs + i) = 2; //Diagonal elements.
    for(int i=0; i<nxs-1; ++i)
        *(mat + i*nxs + i + 1) = -1; //Above diagonal.
    for(int i=1; i<nxs; ++i)
        *(mat + i*nxs + i - 1) = -1; //Below diagonal.

    //Cholesky decomposition
    int nxs2 = nxs*nxs;
    double *l;
    l = (double *) malloc(nxs2 * sizeof(double));
    int ok;
    ok = cholesky(mat, l, nxs);
    if(ok == -1){ //Method safeguard
        printf("Couldn't perform LL' decomposition.\n");
        return 0;
    }

    //LL' = mat Check
    double *lt;
    lt = transpose(l, nxs, nxs);
    double *mat2;
    mat2 = matmul(l, lt, nxs, nxs, nxs);
    double error;
    error = norm(mat, mat2, nxs2);
    printf("Error = norm(flatten(mat), flatten(LL')): %lf\n", error);

    //Solving
    //LL'x = b
    //L(L'x) = b
    double *ltx;
    ltx = solve_l(l, b, nxs);

    //Solving L'x = ltx
    double *x;
    x = solve_u(lt, ltx, nxs);
    printf("x: ");
    print_matrix(x, 1, nxs);

    //Mx = b check
    double *b2;
    b2 = matmul(mat2, x, nxs, nxs, 1);
    error = norm(b, b2, nxs);
    printf("Error = norm(b, b'): %lf\n", error);
    
    printf("Answer:\n");
    printf("%d,%lf\n", n, *(x + (nxs>>1)));

    free(x);
    free(b);
    free(mat);
    free(l);
    free(lt);
    free(ltx);
    free(mat2);
    free(b2);

    return 0;
}
