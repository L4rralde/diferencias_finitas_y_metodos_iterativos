//jacobi.c
/*
Autor: Emmanuel Alejandro Larralde Ortiz | emmanuel.larralde@cimat.mx
Descripción:
    Resuelve un sistema de ecuaciones Ax = b usando el método de Gauss-Seidel

Compilación:
    gcc src/jacobi.c include/matrices/matrices.c -o output/jacobi.o
    chmod +x output/jacobi.o

Uso:
    ./output/jacobi.o <path_matriz> <path_vector>
*/

#include <stdlib.h>
#include <stdio.h>
#include "../include/matrices/matrices.h"

int main(int argc, char **argv){
    double *matrix, *vector;
    int m, n, len;

    //Lee matriz A
    matrix = mat_from_txt(argv[1], &m, &n);
    if(matrix == NULL){
        printf("No se pudo leer el archivo de matrices\n");
        return 0;
    }
    //Lee vector b
    vector = vec_from_txt(argv[2], &len);
    if(matrix == NULL){
        printf("No se pudo leer el archivo de vectores\n");
        return 0;
    }
    if(m != n || n != len){
        printf("No se puede sistema de ecuaciones con diferentes dimensiones\n");
        return 0;
    }

    //Resuelve Ax = b con Jacobi
    double *x;
    x = jacobi(matrix, vector, n);
    printf("Solucion: ");
    print_matrix(x, 1, n);

    //Comprobación
    double *b2;
    b2 = matmul(matrix, x, n, n, 1);
    printf("b': ");
    print_matrix(b2, 1, n);
    printf("b (original): ");
    print_matrix(vector, 1, n);
    double dist;
    dist = norm(vector, b2, n);
    printf("Error = norm(b - b'): %lf\n", dist);

    //Liberación de memoria dinámica
    free(b2);
    free(x);
    free(matrix);
    free(vector);

    return 0;
}
