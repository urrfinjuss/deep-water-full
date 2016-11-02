FSW     = fsw_dw.x

#List of sources
C_FSW   = main.c memory.c

# Compilers -Ofast -ffast-math -mfpmath=387
CC      = gcc
LINK    = gcc
OPT     = -march=native -Ofast -ffast-math -mfpmath=387 -std=gnu99 -Wall -Wpointer-arith -Wcast-align -fno-schedule-insns -fschedule-insns2 -fstrict-aliasing -funroll-loops -fprefetch-loop-arrays

#-----------------------------
#generic

LIB_MPI         =
LIB_FFT         = -L/home/orange/fftw3libd/lib  -lfftw3l_threads -lfftw3l -lm -lpthread
INC_MPI         =
INC_FFT         = -I/home/orange/fftw3libd/include
LIB_ADD         =

#-----------------------------

OBJ_FSW         = $(C_FSW:.c=.o) $(F_FSW:.f=.o)
LIB_FSW         = $(LIB_MPI) $(LIB_FFT) $(LIB_ADD)
INC_FSW         = $(INC_MPI) $(INC_FFT)

#-----------------------------

default: fsw

.f.o:
	$(FC) $(FFLAGS) -c $<


fsw:
	$(CC) $(OPT) $(DEF_FSW) $(INC_FSW) -c $(C_FSW)
	$(LINK) $(OPT) $(OBJ_FSW) $(LIB_FSW) -o $(FSW)
	cp fsw_dw.x ./debug/fsw_dw.x

hostname:
	@echo $(HOSTNAME) $(INC_FFT)

clean:
	@echo "cleaning ..."
	rm -f *~ *.o
	rm -f ./debug/fsw_dw.x

#-------------------------------


