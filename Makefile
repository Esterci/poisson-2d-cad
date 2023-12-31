CC   = gcc
MPICC=mpicc
CFLAGS=-g -O3
LDFLAGS=-lm
OBJS=*.o
BINS= serial_stencil 

# Habilite essa opcao caso tenha compilacao com MPI
USE_MPI = 0

ifeq ($(USE_MPI),1)
	CC=mpicc
	CFLAGS+=-DHAVE_MPI
endif

# Habilite essa opcao caso tenha compilacao com OMP
USE_OMP = 0

ifeq ($(USE_OMP),1)
	CC=mpicc
	CFLAGS+=-DHAVE_OMP -fopenmp
endif

# Habilite essa opcao caso tenha compilacao com Híbrido
USE_HB = 1

ifeq ($(USE_HB),1)
	CC=mpicc
	CFLAGS+=-DHAVE_MPI -fopenmp
endif

all: $(BINS)

%.o: %.c 
	$(CC) $(CFLAGS) $< -c -o $@

serial_stencil: serial_stencil.o save_array.o appctx.o appctx.h Makefile
	$(CC) $(CFLAGS) -o $@ serial_stencil.o save_array.o appctx.o $(LDFLAGS)

mpi_stencil: mpi_stencil.o save_array.o appctx.o appctx.h Makefile
	$(CC) $(CFLAGS) -o $@ mpi_stencil.o save_array.o appctx.o $(LDFLAGS)

omp_stencil: omp_stencil.o save_array.o appctx.o appctx.h Makefile
	$(CC) $(CFLAGS) -o $@ omp_stencil.o save_array.o appctx.o $(LDFLAGS)

hb_stencil: hb_stencil.o save_array.o appctx.o appctx.h Makefile
	$(CC) $(CFLAGS) -o $@ hb_stencil.o save_array.o appctx.o $(LDFLAGS)



clean:
	rm -f $(BINS) $(OBJS)
	rm -f output*.bmp *.txt *.bin *.vti 
