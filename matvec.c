#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#define NREPS 10000

void printMatrix(double *matrix, int m, int n, int rank, char matrixName);

/*
 * Multiplicación de una matriz banda por un vector
 *  w = A*v, con A matriz cuadrada de dimensión N y ancho de banda b
 *  Algoritmo orientado a filas
 */
void matvec(int N, int b, double *A, double *v, double *w) {
  int i, j, li, ls;

  for (i = 0; i < N; i++) {
    w[i] = 0.0;
    li = i - b < 0 ? 0 : i - b;         /* limite inferior */
    ls = i + b > N - 1 ? N - 1 : i + b; /* limite superior */
    for (j = li; j <= ls; j++) {
      w[i] += A[i * N + j] * v[j];
      /*printf("rank = %d -- i = %d, j = %d, N = %d, iGNj = %d, li = %d, ls %d,
       * "*/
      /*"w[%d] = %f, A[%d] = %f, v[%d] = %f\n",*/
      /*-1, i, j, N, i * N + j, li, ls, i, w[i], i * N + j,*/
      /*A[i * N + j], j, v[j]);*/
    }
  }
}

int main(int argc, char **argv) {
  int i, j, k, N = 50, b = 4;
  double *A, *v, *w;

  /* Extracción de argumentos */
  if (argc > 1) { /* El usuario ha indicado el valor de N */
    if ((N = atoi(argv[1])) < 0)
      N = 50;
  }
  if (argc > 2) { /* El usuario ha indicado el valor de b */
    if ((b = atoi(argv[2])) < 0)
      b = 1;
  }
  if (b >= N) { /* Valor de b incorrecto */
    printf("Error: ancho de banda excesivo, N=%d, b=%d\n", N, b);
    exit(1);
  }

  /* Reserva de memoria */
  A = (double *)calloc(N * N, sizeof(double));
  v = (double *)calloc(N, sizeof(double));
  w = (double *)calloc(N, sizeof(double));

  /* Inicializar datos */
  for (i = 0; i < N; i++) {
    A[i * N + i] = 2 * b;
    /*printf("i * N + i = %d, i = %d, N = %d, i = %d\n", i * N + i, i, N, i);*/
  }
  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++) {
      if (i != j && abs(i - j) <= b) {
        A[i * N + j] = -1.0;
      }
    }
  }
  // printMatrix(A, N, N, 0, 'A');

  for (i = 0; i < N; i++)
    v[i] = 1.0;
  // printMatrix(v, 1, N, 0, 'v');

  clock_t start = clock();
  /* Multiplicación de matrices */
  for (k = 0; k < NREPS; k++)
    matvec(N, b, A, v, w);
  clock_t end = clock();
  float seconds = (float)(end - start) / CLOCKS_PER_SEC;

  /* Imprimir solución */
  if (N < 100)
    for (i = 0; i < N; i++)
      printf("w[%d] = %g\n", i, w[i]);

  printf("Tiempo transcurrido: %7.5f\n", seconds);
  free(A);
  free(v);
  free(w);

  return 0;
}

void printMatrix(double *matrix, int m, int n, int rank, char matrixName) {
  int i, j;

  printf("Process--> %d\n", rank);
  /*printf("Matrix--> %c\n", matrixName);*/
  /*printf("m = %d || n = %d \n", m, n);*/

  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      printf("%7.3f ", matrix[i * n + j]);
    }
    printf(" \n");
  }
  printf("------------------------\n");
}
