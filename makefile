ORIGINAL_OBJS	= matvec.o
ORIGINAL_OUT	= matvec.exe
CC	 = gcc
MPI_OBJS	= matvecmpi.o
MPI_OUT	= matvecmpi.exe
MPI_CC	 = mpicc
FLAGS	 = -g -c -Wall -O3
LFLAGS	 = -Wall -O3

original: $(ORIGINAL_OBJS)
	$(CC) -g $(ORIGINAL_OBJS) -o $(ORIGINAL_OUT) $(LFLAGS)

matvec.o: matvec.c
	$(CC) $(FLAGS) matvec.c 

mpi: $(MPI_OBJS)
	$(MPI_CC) -g $(MPI_OBJS) -o $(MPI_OUT) $(LFLAGS)

matvecmpi.o: matvecmpi.c
	$(MPI_CC) $(FLAGS) matvecmpi.c 


clean:
	rm -f $(ORIGINAL_OBJS) $(ORIGINAL_OUT)