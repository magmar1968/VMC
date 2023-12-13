BINDIR := ../bin
SRCDIR := src
OBJDIR := obj

SRCEXT := .F90
OBJEXT := .o

FC      := gfortran
FFLAGS  := 
LDFLAGS := 
LIB    :=  
.SUFFIXES: .f90 .o .F90

#-----------------------------------------
OTHER :=  src/prec src/string src/iofile \
		src/random \
		src/arrays_utility \
		src/HS_puregas \
		src/VMC_params \
		src/VMC_print \
		src/VMC \

PROGS := src/twf_param_opt 

ALL_FILES  = $(OTHER) $(PROGS)
ALL_SOURCE = $(ALL_FILES:%=%.F90)
ALL_OBJS   := $(subst $(SRCDIR)/,$(OBJDIR)/,$(ALL_FILES:%=%.o))



.PHONY: all
all: $(PROGS)

.PHONY: remake
remake: clean
remake: all

.PHONY: debug
debug: clean
debug: FFLAGS += -g -W -Wall -fcheck=all #-fsanitize=address,null,undefined
debug: LDFLAGS+= -g -W -Wall -fcheck=all #-fsanitize=address,null,undefined
debug: all

.PHONY:parallel
parallel: clean
parallel: FFLAGS += -fopenmp
parallel: LDFLAGS+= -fopenmp
parallel: all

.PHONY: fast
fast: clean
fast: FFLAGS += -O3
fast: all

.PHONY: profiling
profiling: clean
profiling: FFLAGS += -pg
profiling: LDFLAGS+= -pg
profiling: all


$(PROGS): $(ALL_OBJS)
	$(FC) $(LDFLAGS) $(ALL_OBJS) -o $(notdir $@).x 

$(OBJDIR)/%.o: $(SRCDIR)/%$(SRCEXT)
	$(FC) $(FFLAGS) -c -o $@ $<


.PHONY: clean
clean:
	-rm -f obj/*.o *.mod *.x