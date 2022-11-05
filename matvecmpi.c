#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#define NREPS 1000000

void printMatrix(double *matrix, int m, int n, int rank, char matrixName);

/*
 * Multiplicación de una matriz banda por un vector
 *  w = A*v, con A matriz cuadrada de dimensión N y ancho de banda b
 *  Algoritmo orientado a filas
 */
void matvec(int N, int b, double *A, double *v, double *w) {
  int i, j, li, ls, rank, above, below, size, iG, n, jaux;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  n = N / size;

  if (rank == 0)
    above = MPI_PROC_NULL;
  else
    above = rank - 1;
  if (rank == size - 1)
    below = MPI_PROC_NULL;
  else
    below = rank + 1;

  MPI_Sendrecv(&v[b], b, MPI_DOUBLE, above, 0, &v[0], b, MPI_DOUBLE, above, 0,
               MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  MPI_Sendrecv(&v[n], b, MPI_DOUBLE, below, 0, &v[b + n], b, MPI_DOUBLE, below,
               0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  /*printMatrix(v, 1, n + 2 * b, rank, 'v');*/
  for (i = 0; i < n; i++) {
    iG = rank * n + i;
    w[i] = 0.0;
    li = iG - b < 0 ? 0 : iG - b;         /* limite inferior */
    ls = iG + b > N - 1 ? N - 1 : iG + b; /* limite superior */
    for (j = li; j <= ls; j++) {
      jaux = j;
      if (rank == 0) {
        jaux = j + b;
      }
      w[i] += A[iG * N + j] * v[jaux];
      /*printf("rank = %d -- i = %d, j = %d, N = %d, iGNj = %d, li = %d, ls %d,
       * "*/
      /*"ls - li = %d, "*/
      /*"w[%d] = %f, A[%d] = %5.3f, v[%d] = %5.3f\n",*/
      /*rank, iG, j, N, i * N + j, li, ls, ls - li, i, w[i], iG * N + j,*/
      /*A[iG * N + j], jaux, v[jaux]);*/
    }
  }
}

int main(int argc, char **argv) {
  int i, j, k, N = 50, b = 4, rank, size, n, iG;
  double *A, *v, *w, *Aux, time;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

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

  n = N / size;
  /*printf("N = %d, n = %d, b = %d, size = %d, rank = %d\n", N, n, b, size,
   * rank);*/
  /* Reserva de memoria */
  A = (double *)calloc(N * N, sizeof(double));
  v = (double *)calloc(N + 2 * b, sizeof(double));
  /*printf("n + 2 * b = %d\n", n + 2 * b);*/
  w = (double *)calloc(n, sizeof(double));
  Aux = (double *)calloc(N, sizeof(double));

  /* Inicializar datos */
  for (i = 0; i < n; i++) {
    iG = rank * n + i;
    A[iG * N + iG] = 2 * b;
  }

  for (i = 0; i < n; i++) {
    iG = rank * n + i;
    for (j = 0; j < N; j++) {
      if (iG != j && abs(iG - j) <= b) {
        A[iG * N + j] = -1.0;
      }
    }
  }
  /*printMatrix(A, N, N, rank, 'A');*/

  for (i = 0; i < N + 2 * b; i++)
    v[i] = 1.0;

  /* Multiplicación de matrices */
  MPI_Barrier(MPI_COMM_WORLD);
  time = MPI_Wtime();
  for (k = 0; k < NREPS; k++)
    matvec(N, b, A, v, w);
  MPI_Barrier(MPI_COMM_WORLD);
  time = MPI_Wtime() - time;

  /*printMatrix(w, 1, n, rank, n);*/
  MPI_Gather(&w[0], n, MPI_DOUBLE, Aux, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  if (rank == 0) {
    if (N < 100)
      for (i = 0; i < N; i++) {
        printf("Aux[%d] = %g\n", i, Aux[i]);
      }
    free(Aux);
    printf("Tiempo de ejecución: %7.5f\n", time);
  }

  /*printMatrix(w, 1, n, rank, 'w');*/
  /*if (N < 100) {*/
  /*for (i = 0; i < n; i++)*/
  /*printf("rank = %d -- w[%d] = %g\n", rank, i, w[i]);*/
  /*}*/
  MPI_Finalize();
  free(A);
  free(v);
  free(w);

  return 0;
}

void printMatrix(double *matrix, int m, int n, int rank, char matrixName) {
  int i, j;

  printf("Process--> %d || Matrix --> %c\n", rank, matrixName);

  for (i = 0; i < m; i++) {
    printf("r%d -- ", rank);
    for (j = 0; j < n; j++) {
      printf("%7.3f ", matrix[i * n + j]);
    }
    printf(" \n");
  }
  printf("------------------------\n");
}
