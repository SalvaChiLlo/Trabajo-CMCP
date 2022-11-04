#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#define NREPS 10000

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
    w[i] = 0.0;
    iG = rank * n + i;
    li = iG - b < 0 ? 0 : iG - b;         /* limite inferior */
    ls = iG + b > N - 1 ? N - 1 : iG + b; /* limite superior */
    for (j = li; j <= ls; j++) {
      jaux = j;
      if (rank == 0) {
        jaux = j + b;
      }
      w[i] += A[iG * N + j] * v[jaux];
      /*printf("rank = %d -- i = %d, iG = %d, j = %d, N = %d, iGNj = %d, li =
       * "*/
      /*"%d, ls %d, "*/
      /*"w[%d] = %f, A[%d] = %f, v[%d] = %f\n",*/
      /*rank, i, iG, j, N, iG * N + j, li, ls, iG, w[i], iG * N + j,*/
      /*A[iG * N + j], j, v[j]);*/
    }
  }
}

int main(int argc, char **argv) {
  int i, j, k, N = 50, b = 4, rank, size, n, iG;
  double *A, *v, *w, *Aux;

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
  printf("n = %d\n", n);
  /* Reserva de memoria */
  A = (double *)calloc(N * N, sizeof(double));
  v = (double *)calloc(n + 2 * b, sizeof(double));
  w = (double *)calloc(n, sizeof(double));
  Aux = (double *)calloc(N*N, sizeof(double));

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

  for (i = 0; i < n; i++)
    v[b + i] = 1.0;

  /* Multiplicación de matrices */
  for (k = 0; k < NREPS; k++)
    matvec(N, b, A, v, w);

  printf("w[2] = %f\n", w[2]);
  MPI_Gather(w, n, MPI_DOUBLE, Aux, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  if (rank == 0)
    if (N < 100)
      for (i = 0; i < N; i++)
        printf("Aux[%d] = %g\n", i, Aux[i]);

  /*if (N < 100) {*/
  /*printf("PROCESO --> %d\n", rank);*/
  /*for (i = 0; i < n+99; i++)*/
  /*printf("rank = %d -- w[%d] = %g\n", rank, i, w[i]);*/
  /*}*/
  free(A);
  free(v);
  free(w);
  free(Aux);
  MPI_Finalize();

  return 0;
}

void printMatrix(double *matrix, int m, int n, int rank, char matrixName) {
  int i, j;

  printf("Process--> %d\n", rank);

  for (i = 0; i < m; i++) {
    printf("r%d -- ", rank);
    for (j = 0; j < n; j++) {
      printf("%7.3f ", matrix[i * n + j]);
    }
    printf(" \n");
  }
  printf("------------------------\n");
}
