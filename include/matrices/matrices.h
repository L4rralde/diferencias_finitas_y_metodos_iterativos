//matrices.h
/*
Autor: Emmanuel Alejandro Larralde Ortiz    | emmanuel.larralde@cimat.mx
Descripcion:
    Header de una librería con funciones varias para manejar matrices y vectores, 
    y resolver sistemas de ecuaciones.
*/
#ifndef MATRICES_H
#define MATRICES_H

double *vec_from_txt(char *, int *);
double *mat_from_txt(char *, int *, int *);
void print_matrix(double *, int, int);

double *matmul(double *, double *, int, int, int);
double *solve_d(double *, double *, int);
double *solve_u(double *, double *, int);
double *solve_l(double *, double *, int);

double norm(double *, double *, int);

int lu_crout(double *, double *, double *, int);
int cholesky(double *, double *, int);
int cholesky_ldl(double *, double *, double *, int);

double *transpose(double *, int, int);

double *gauss_seidel(double *, double *, int);
double *jacobi(double *, double *, int);

#endif