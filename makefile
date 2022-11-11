ORIGINAL_OBJS	= matvec.o
ORIGINAL_OUT	= matvec.exe
CC	 = gcc
MPI_OBJS	= matvecmpi.o
MPI_OUT	= matvecmpi.exe
MPIS_OBJS	= matvecmpis.o
MPIS_OUT	= matvecmpis.exe
MPI_CC	 = mpicc
FLAGS	 = -g -c -Wall -O3
LFLAGS	 = -Wall -O3

original: $(ORIGINAL_OBJS)
	$(CC) -g $(ORIGINAL_OBJS) -o $(ORIGINAL_OUT) $(LFLAGS)

matvec.o: matvec.c
	$(CC) $(FLAGS) matvec.c 

mpi: $(MPI_OBJS)
	$(MPI_CC) -g $(MPI_OBJS) -o $(MPI_OUT) $(LFLAGS)

mpis: $(MPIS_OBJS)
	$(MPI_CC) -g $(MPIS_OBJS) -o $(MPIS_OUT) $(LFLAGS)

matvecmpi.o: matvecmpi.c
	$(MPI_CC) $(FLAGS) matvecmpi.c 

matvecmpis.o: matvecmpis.c
	$(MPI_CC) $(FLAGS) matvecmpis.c 


clean:
	rm -f $(ORIGINAL_OBJS) $(ORIGINAL_OUT)