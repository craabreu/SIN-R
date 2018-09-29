# Compilers and their basic options:
FORT = gfortran

BASIC_OPTS = -march=native -m64 -fopenmp -cpp -fmax-errors=1
BASIC_OPTS += -Wall -Wno-maybe-uninitialized
FAST_OPTS = -Ofast
DEBUG_OPTS = -g -Og -fcheck=all -Ddebug
ifeq ($(DEBUG), 1)
  FLAGS = $(BASIC_OPTS) $(DEBUG_OPTS)
else
  FLAGS = $(BASIC_OPTS) $(FAST_OPTS)
endif

SRCS = $(wildcard *.f90)
OBJS = $(SRCS:%.f90=%.o)
LIBFILE = libmscommon.a

all: lib

lib: $(LIBFILE)

$(LIBFILE): $(OBJS)
	ar -cr $(LIBFILE) $(OBJS)

%.o: %.f90 mGlobal.o
	$(FORT) $(FLAGS) -c -o $@ $<

mGlobal.o: mGlobal.f90
	$(FORT) $(FLAGS) -c -o $@ $<

clean:
	rm -f *.o
	rm -f *.mod
	rm -f $(LIBFILE)
	rm -rf html/

