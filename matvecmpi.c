#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#define NREPS 10000

/*
 * Multiplicación de una matriz banda por un vector
 *  w = A*v, con A matriz cuadrada de dimensión N y ancho de banda b
 *  Algoritmo orientado a filas
 */
void matvec(int N, int b, double *A, double *v, double *w) {
  int i, j, li, ls, rank, above, below, size, iG;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

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

  MPI_Sendrecv(&v[N], b, MPI_DOUBLE, below, 0, &v[b + N], b, MPI_DOUBLE, below,
               0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  for (i = 0; i < N; i++) {
    iG = N * rank + i;
    w[i] = 0.0;
    li = iG - b < 0 ? 0 : iG - b;                       /* limite inferior */
    ls = iG + b > N * size - 1 ? N * size - 1 : iG + b; /* limite superior */
    for (j = li; j <= ls; j++) {
      w[i] += A[i * N * size + j] * v[j];
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
  /* Reserva de memoria */
  A = (double *)calloc(N * N, sizeof(double));
  v = (double *)calloc(n + 2 * b, sizeof(double));
  w = (double *)calloc(n, sizeof(double));
  Aux = (double *)calloc(N, sizeof(double));

  /* Inicializar datos */
  for (i = 0; i < n; i++) {
    iG = n * rank + i;
    A[i * N + iG] = 2 * b;
  }

  /*  for (i = 0; i < n; i++) {*/
  /*iG = n * rank + i;*/
  /*li = iG - b < 0 ? 0 : iG - b;*/
  /*ls = iG + b > N - 1 ? N - 1 : iG + b;*/
  /*for (j = li; j <= ls; j++) {*/
  /*A[i * N + j] = -1.0;*/
  /*}*/
  /*}*/
  for (i = 0; i < n; i++) {
    iG = n * rank + i;
    for (j = 0; j < N; j++) {
      if (iG != j && abs(iG - j) <= b)
        A[i * N + j] = -1.0;
    }
  }

  for (i = 0; i < n; i++)
    v[b + i] = 1.0;

  /* Multiplicación de matrices */
  for (k = 0; k < NREPS; k++)
    matvec(n, b, A, v, w);

  /*MPI_Gather(w, n, MPI_DOUBLE, Aux, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);*/
  /*if (rank == 0)*/
  /*[> Imprimir solución <]*/
  /*if (N < 100)*/
  /*[>printf("PROCESO --> %d\n", rank);<]*/
  /*for (i = 0; i < N; i++)*/
  /*printf("w[%d] = %g\n", i, Aux[i]);*/

  if (N < 100) {
    printf("PROCESO --> %d\n", rank);
    for (i = 0; i < n; i++)
      printf("w[%d] = %g\n", i, w[i]);
  }
  free(A);
  free(v);
  free(w);
  MPI_Finalize();

  return 0;
}
