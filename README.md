# Proyecto CMCP
Conceptos y Métodos de la Computación Paralela
![Figura 1](Figura1.png)

## Producto matriz-vector con matriz banda
Para este proyecto se proporciona la versión secuencial **matvec.c** del producto matriz-vector *w = Av* de una matriz banda cuadrada *A ∈ R^N×N* con ancho de banda b. La estructura banda de la matriz se refiere a que el elemento *a_ij* es cero si *|i − j| > b*.

En una implementación real, únicamente se almacenaría la parte de la matriz correspondiente a la banda (elementos no nulos). Sin embargo, para simplificar, el programa *matvec.c* utiliza un almacenamiento denso por filas, es decir, se reserva espacio también para los ceros.

El producto matriz-vector es la operación básica utilizada en muchos métodos iterativos para diferentes problemas de álgebra lineal. En estos métodos se repite la operación de producto matriz-vector muchas veces. Para simplificar, en nuestro caso vamos a analizar las prestaciones de un bucle simple que ejecute el producto matriz-vector (secuencial o paralelo) un número determinado de veces, por ejemplo 1000 o 10000.

**Versión paralela básica** Habría que desarrollar una versión paralela *matvecmpi.c*, en la que la distribución de la matriz se realice por bloques de filas contiguas, siendo el número de filas locales a cada proceso, *n*, aproximadamente el mismo pero no necesariamente igual en todos los procesos. Para simplificar la paralelización, se puede asumir el caso sencillo en que el ancho de banda es siempre menor que el tamaño local de la matriz, *n*.

En la Figura 1 (izquierda) se muestra el esquema que debe seguir el código paralelo. Básicamente, antes de realizar el cálculo de su porción local del vector *w*, cada proceso debe recibir *b* elementos del vector *v* del vecino de arriba y *b* elementos del vecino de abajo. Asimismo, cada proceso deberá enviar *b* elementos de su vector *v* al proceso de arriba y *b* elementos al proceso de abajo. Esta comunicación debe realizarse de forma que se evite la secuencialización y/o el interbloqueo, por ejemplo, utilizando primitivas de comunicación combinada (*MPI_Sendrecv*). Habitualmente, la implementación paralela se haría de forma que cada proceso almacena *n* filas completas de la matriz y *n* elementos del vector *w*. En el caso del vector *v*, se reservaría memoria para *n + 2b* elementos, es decir, los *n* elementos locales más *b* recibidos de cada uno de los vecinos.
En el estudio de prestaciones, se debería considerar tanto el tamaño de la matriz *N* como el ancho de banda *b*, ya que este último parámetro influye tanto en el coste aritmético como en la longitud de los mensajes.

### Versión MPI con solapamiento
Además de la versión paralela básica descrita antes, el objetivo de este proyecto es hacer una versión que intente mejorar las prestaciones paralelas solapando la comunicación con parte de la computación. Esto es posible porque hay parte del cálculo (correspondiente a las columnas que coinciden con el rango local de filas) que no depende de los datos de otros procesos. Para ello, hay que utilizar primitivas de comunicación no bloqueantes. El algoritmo quedaría así:
```pseudocode
matvec_parallel(N,n,b,A,v,w):
  start send to neighbour above
  start send to neighbour below
  start receive to neighbour above
  start receive to neighbour below
  compute partially the local vector w, using only local columns of A and local elements of v
  finish all pending communication operations (MPI_Waitall)
  complete computation of w for remaining columns of A
```
Realizar un estudio de prestaciones de la versión básica y de la versión mejorada, para ver si en este caso se nota una mejora sustancial o no.