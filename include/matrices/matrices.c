//matrices.c
/*
Autor: Emmanuel Alejandro Larralde Ortiz    | emmanuel.larralde@cimat.mx
Descripcion:
    Libreria con funciones varias para manejar matrices y vectores, 
    y resolver sistemas de ecuaciones.
*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "matrices.h"


//Lee un vector almacenado en un archivo txt.
double * vec_from_txt(char *fname, int *len){
    //fname: ruta del archivo, len: puntero para almacenar longitud del vector.
    FILE *file = fopen(fname, "r");
    if (file == NULL)
        return NULL;
    fscanf(file, "%d", len);

    double * ptr = (double *) malloc(sizeof(double) * (*len));
    for(int i=0; i<(*len); ++i)
        fscanf(file, "%lf", ptr + i);

    fclose(file);
    return ptr;
}


//Lee una matriz almacenada en un archivo de texto
double * mat_from_txt(char *fname, int *m, int *n){
    //fname: ruta del archivo.
    //m: puntero al número de renglones.
    //n: puntero al número de columnas.
    FILE *file = fopen(fname, "r");
    if (file == NULL)
        return NULL;

    fscanf(file, "%d %d", m, n);

    int len = (*m) * (*n);
    double * ptr = (double *) malloc(sizeof(double) * len);
    for(int i=0; i<len; ++i)
        fscanf(file, "%lf", ptr + i);

    fclose(file);
    return ptr;
}


//Imprime el contenido de una matriz.
void print_matrix(double *mat, int m, int n){
    //mat: puntero a la matriz, m: número de renglones, n: número de columnas.
    for(int i=0; i<m; ++i){
        for(int j=0; j<n; ++j)
            printf("%lf ", *(mat + i*n + j));
        printf("\n");
    }
}

//Genera una nueva matriz (memoria dinámica) con el mismo contenido que una matriz dada mat
double *copy_matrix(double *mat, int m, int n){
    int len = m * n;
    double *copy;
    copy = (double *) malloc(len * sizeof(double));
    for(int i=0; i<len; ++i)
        *(copy +i) = *(mat +i);
    return copy;
}

//Copia los elementos de una matriz mat a una matrix copy.
void copy_matrix_no_alloc(double *mat, double *copy, int m, int n){
    int len = m*n;
    for(int i=0; i<len; ++i)
        *(copy + i) = *(mat + i);
}

//Multiplica dos matrices, retorna un puntero a la matriz resultante.
double *matmul(double * mat_a, double * mat_b, int m, int n, int o){
    //mat_a: puntero a la matriz de la izquierda (en la multiplicación).
    //mat_b: puntero a la matriz de la derecha (en la multiplicación).
    //m: número de renglones de la primera matriz.
    //n: número de columnas de la primera y de renglones de la segunda.
    //o: número de columnas de la segunda matriz.
    int max = m * o;
    double *mat_res = (double *) malloc(sizeof(double) * max);
    
    for(int i=0; i < max; ++i)
        *(mat_res + 1) = 0.0;
    
    for(int i=0; i<m; ++i)
        for(int j=0; j<o; ++j)
            for(int k=0; k<n; ++k)
                *(mat_res + i*o + j) += (*(mat_a + i*n + k))*(*(mat_b + k*o + j));
    
    return mat_res;
}


//Resuelve un sistema de ecuaciones de la forma Dx = b. Con D diagonal.
double * solve_d(double *mat, double *b, int n){
    //mat: puntero a la matriz diagonal D.
    //b: puntero al vector b.
    //n: número de renglones y columnas de la matriz, y longitud del vector.
    double *x = (double *) malloc(sizeof(double) *  n);
    for(int i=0; i<n; ++i)
        *(x + i) = (*(b + i))/(*(mat + i*(n + 1)));
    return x;
}


//Resuelve un sistema de ecuaciones de la forma Ux = b. Con U triangular superior.
double * solve_u(double *matrix, double *vector, int n){
    //martrix: puntero a la matriz triangular superior.
    //vector: puntero al vector b.
    //n: número de renglones y columnas de la matriz, y longitud del vector.
    double *x;
    x = (double *) malloc(sizeof(double) * n);
    for(int i=0; i<n; ++i)
        *(x + i) = 0.0;

    double sum;
    for(int i=n-1; i>=0; --i){
        sum = 0.0;
        for(int j=i+1; j<n; ++j)
            sum += *(matrix + i*n + j) * *(x + j);
        *(x + i) = (*(vector + i) - sum)/(*(matrix + i*(n + 1)));
    }
    return x;
}


//Resuelve un sistema de ecuaciones de la forma Lx = b. Con L triangular inferior.
double * solve_l(double *matrix, double *vector, int n){
    //martrix: puntero a la matriz triangular inferior.
    //vector: puntero al vector b.
    //n: número de renglones y columnas de la matriz, y longitud del vector.
    double *x;
    x = (double *) malloc(sizeof(double) * n);
    for(int i=0; i<n; ++i)
        *(x + i) = 0.0;

    double sum;
    for(int i=0; i<n; ++i){
        sum = 0.0;
        for(int j=0; j<i; ++j)
            sum += *(matrix + i*n +j) * *(x + j);
        *(x + i) = (*(vector + i) - sum)/(*(matrix + i*(n + 1)));
    }
    return x;
}

//Calcula el cuadrado de la distancia entre dos vectores.
double norm(double *vec_a, double *vec_b, int n){
    //vec_a: puntero de un vector.
    //vec_b: puntero del otro vector.
    //n: longitud de los punteros
    double sum = 0.0;
    double diff;
    for(int i=0; i<n; ++i){
        diff = *(vec_a + i) - *(vec_b + i);
        sum += diff * diff;
    }

    return sum;
}

//Factoriza una matriz cuadrada A nxn en dos matrices L y U por el método de Crout.
int lu_crout(double *A, double *L, double *U, int n){
    //A: matriz nxn a fatorizar.
    //L: Matriz triangular inferior resultado de la factorización.
    //U: Matriz triangular superior resultado de la factorización.
    //n: dimensión de la matriz cuadrada.
    int max = n*n;
    double acc;
    *U = 1; //U[0][0] = 1 
    *L = *A; //L[0][0] = A[0][0]
    for(int i=1; i<n; ++i){
        //Resuelve L[0:i-1][0:i-1]U[0:i-1][i] = A[0:i-1][i]
        for(int ui=0; ui<i; ++ui){
            if(*(L + ui*n + ui) == 0)
                return -1;
            acc = 0;
            for(int k=0; k<ui; ++k)
                acc += (*(L + ui*n + k)) * (*(U + k*n +i));
            *(U + ui*n + i) = (*(A + ui*n + i) - acc)/(*(L + ui*n + ui));
        }
        //U[i][i] = 1;
        *(U + i*n + i) = 1.0;
        //Resuelve U[1:i][1:i]^TL[i][1:i]^T = A[i][1:i]^T
        for(int j=0; j<=i; ++j){
            acc = 0;
            for(int k=0; k<j; ++k)
                acc += (*(L + i*n + k)) * (*(U + k*n +j));
            *(L +i*n + j) = *(A +i*n + j) - acc;
        }
    }
    return 0;
}

//Factoriza una matriz simétrica cuadrada nxn definida positiva en LL' = A.
int cholesky(double *mat, double *l, int n){
    //mat: matriz a factorizar.
    //l: matriz triangular inerior talque ll' = mat
    //n: dimensión de la matriz cuadrada
    double acc, aux;
    for(int j=0; j<n; ++j){ //Por cada columna
        //Cálculo de L[j][j]
        acc = 0;
        for(int k=0; k<j; k++)
            acc += (*(l + j*n +k)) * (*(l + j*n +k));
        aux = *(mat + j*n + j) - acc;
        if(aux <= 0) //Nan safeguard.
            return -1;
        *(l + j*n + j) = sqrt(aux); //L[j][j]
        for(int i=j+1; i<n; ++i){ //Para cada elemento de la columna debajo de la diagonal.
            //Cálculo de L[i][j], i>j
            acc = 0;
            for(int k=0; k<j; ++k)
                acc += (*(l + j*n + k)) * (*(l + i*n + k));
            *(l + i*n + j) = (*(mat + j*n + i) - acc)/(*(l + j*n + j));
        }
    }
    return 0;
}

//Factoriza una matriz simétrica cuadrada nxn definida positiva mat en mat = ldl'
int cholesky_ldl(double *mat, double *l, double *d, int n){
    //mat: matriz a factorizar.
    //l: matriz triangular inerior resultado de la factorización.
    //d: matriz diagonal resultado de la factorización.
    //n: dimensión de la matriz.
    double acc, aux;
    for(int j=0; j<n; ++j){ //Por cada columna
        //Cálculo de D[j][j]
        acc = 0;
        for(int k=0; k<j; ++k)
            acc += (*(d + k*n +k)) * (*(l + j*n +k))*(*(l + j*n +k));
        aux = *(mat + j*n + j) - acc;
        if(aux == 0) //Nan safeguard
            return -1;
        *(d + j*n + j) = aux; //D[j][j]
        *(l + j*n + j) = 1; //L[j][j]
        for(int i=j+1; i<n; ++i){ //Para cada elemento de la columna y debajo de la diagonal.
            //Calculo de L[i][j], i>j
            acc = 0;
            for(int k=0; k<j; ++k)
                acc += (*(d + k*n + k)) * (*(l + i*n + k)) * (*(l + j*n + k));
            *(l +i*n + j) = (*(mat + i*n +j) - acc)/aux;
        }
    }
    return 0;
}

//Calcula la transpuesta de una matriz mxn
double *transpose(double *mat, int m, int n){
    //mat: matriz.
    //m: numero de renglones
    //n: numero de columnas
    double *mat2;
    mat2 = (double *) malloc(m * n * sizeof(double));
    for(int i=0; i<m; ++i)
        for(int j=0; j<n; ++j)
            *(mat2 + j*n + i) = *(mat + i*n +j); //mat2[j][i] = mat[i][j]
    return mat2;
}

//Realiza el metodo de Gauss-Seidel usando los valores de xs mas recientes 
//o de la estimacion anterior
double *_gauss_seidel(double *mat, double *b, int n, int no_stale){
    double *x;
    x = (double *) calloc(n, sizeof(double));
    //Safeguard: non-0s check
    for(int i=0; i<n; ++i){
        if(*(mat + i*n + i) == 0){
            printf("EXCEPTION: found zeros in diagonal\n");
            return x;
        }
    }
    double *prev_x, *xptr, acc, error;
    prev_x = (double *) calloc(n, sizeof(double));
    //if no_stale==1:
    //  Usa los valores más recientes de x
    //else:
    //  Usa los valores de la iteración anterior.
    xptr = (no_stale == 1)? x: prev_x;
    for(int rep=0; rep<2000; rep++){
        copy_matrix_no_alloc(x, prev_x, 1, n);
        for(int i=0; i<n; ++i){//Para cada elemento de x.
            acc = 0;
            for(int j=0; j<n; ++j){//Para cada elemento de x...
                if(i==j) //Excepto xi (o sea, el mismo elemento de x)
                    continue;
                acc += (*(mat + i*n + j)) * (*(xptr + j)); //mat[i][j]*mat[j]
            }
            //x[i] = (b[i] - acc)/(mat[i][i])
            *(x + i) = (*(b + i) - acc)/(*(mat + i*n +i));
        }
        //Calcula la distancia entre la iteración actual y la anterior
        error = norm(x, prev_x, n); 
        if(error <= 0.0001) //Si no cambia mucho, termina.
            break;
    }
    if(error > 0.0001)
        printf("WARNING: Couldn't find a solution\n");

    free(prev_x);
    return x;
}

//Obtiene la solucion de Ax=b usando el método de Gauss-Seidel
double *gauss_seidel(double *mat, double *b, int n){
    return _gauss_seidel(mat, b, n, 1);
}

//Obtiene la solucion de Ax=b usando el método de Jacobi
double *jacobi(double *mat, double *b, int n){
    return _gauss_seidel(mat, b, n, 0);
}