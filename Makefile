# Compilers and their basic options:
FORT ?= gfortran
EXECUTABLE ?= nvt

BASIC_OPTS = -march=native -m64 -fopenmp -cpp -fmax-errors=1
BASIC_OPTS += -Wall -Wno-maybe-uninitialized
FAST_OPTS = -Ofast
DEBUG_OPTS = -g -Og -fcheck=all -Ddebug
ifeq ($(DEBUG), 1)
  FLAGS = $(BASIC_OPTS) $(DEBUG_OPTS)
else
  FLAGS = $(BASIC_OPTS) $(FAST_OPTS)
endif

ifdef EMDEE
  LN_OPTS = -I$(EMDEE)/include $(LN_SO_OPT) -Wl,-rpath,'$(EMDEE)/lib'
endif

#FORT = gfortran
#FLAGS = -fno-range-check -march=native -ffast-math -funroll-loops -fstrict-aliasing -O3 -Wunused -cpp -fopenmp
COMMON = common
SRCS = $(wildcard $(COMMON)/*.f90)
LIB = mscommon
LIBFILE = $(COMMON)/lib$(LIB).a

.PHONY: all lj lj_coul

all: lj_coul

lj: $(EXECUTABLE)_lj

lj_coul: $(EXECUTABLE)_lj_coul

$(EXECUTABLE)_lj_coul: $(EXECUTABLE).f90 $(LIBFILE) $(SRCS)
	$(FORT) $(FLAGS) $(LN_OPTS) -Dcoul -J$(COMMON) -L$(COMMON) -o $@ $< -l$(LIB) -lemdee

$(EXECUTABLE)_lj: $(EXECUTABLE).f90 $(LIBFILE) $(SRCS)
	$(FORT) $(FLAGS) -J$(COMMON) -L$(COMMON) -o $@ $< -l$(LIB) -lemdee

$(LIBFILE): $(SRCS)
	make -C $(COMMON) lib

clean:
	rm -f $(EXECUTABLE)_lj
	rm -f $(EXECUTABLE)_lj_coul

clean-all:
	rm -f $(EXECUTABLE)_lj
	rm -f $(EXECUTABLE)_lj_coul
	make -C $(COMMON) clean

