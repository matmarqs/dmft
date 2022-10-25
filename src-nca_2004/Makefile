PLATFORM = PC
RELEASE = 1
include Makefile.$(PLATFORM)

executable = $(DIR)/nca
generate = $(DIR)/generate
fobjects =
cobjects = $(DIR)/Main.o
objects = $(fobjects) $(cobjects)

base_headers = Common.h average.h function.h mesh.h

util_headers = assert.h  complex.h  util.h zeroin.h parser.h  timer.h

headers = $(base_headers) $(util_headers)

OFLAGS = $(OPT) -DNO_ARG_CHECK
GFLAGS = -g -DCTMA_DEBUG
ADFLAGSCC = `if [ x$(RELEASE) != x ]; then echo $(OFLAGS); else echo $(GFLAGS); fi`
ADFLAGS = $(OPT)

FFLAGS  = $(ADFLAGS)   $(PFFLAGS)
CFLAGS  = $(ADFLAGSCC) $(PFFLAGSCC)
LDFLAGS = $(ADFLAGSCC) $(PFFLAGSCC) $(LDPFLAGS) $(LIBPATH)
cFLAGS  = $(ADFLAGSCC) $(PFFLAGSCC)

.SUFFIXES : .cc
.cc.o:  $(base_headers)
	$(C++) $(CFLAGS) -c $<

.SUFFIXES : .c
.c.o:  $(base_headers)
	$(cc) $(cFLAGS) -c $<

.SUFFIXES : .f
.f.o:
	$(F77) $(FFLAGS) -c $<

all : $(executable) $(generate)

$(executable) : $(objects)
	$(C++) $(LDFLAGS) -o $@ $(objects) $(LIBS)

$(DIR)/generate : generate.cc
	$(C++) $(LDFLAGS) -o $(generate) generate.cc

$(DIR)/Main.o : Main.cc $(headers)
	$(C++) $(CFLAGS) -c -o $@ Main.cc

$(DIR)/mrqmin.o : mrqmin.c
	$(cc) $(CFLAGS) -c -o $@ mrqmin.c

$(DIR)/erfc.o : erfc.f
	$(F77) $(FFLAGS) -c -o $@ $?
clean :
	- rm -f $(objects) $(executable) $(generate)

