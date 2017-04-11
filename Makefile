PROF =
WARN = -Wall -Wstrict-aliasing=2
CFLAGS += -std=c++11
DFLAGS = -g -DS_DEBUG #-ftrapv
OFLAGS = 

TARGET = test_dir

CFLAGS = $(WARN) $(DFLAGS) $(OFLAGS) $(PROF) -pipe $(GUI_FLAGS)
IFLAGS = -I$(HOME)/include -I../common
LFLAGS = $(CFLAGS) -L$(HOME)/lib
LIBS = -lm -lstdc++ -lgsl -lgslcblas

DOC_INSTALL = cd doc/html; installdox -l biome.tag@/home/martin/doc/Biome/html;
DOC_SOURCE = *.cc *.h

DEP_SRC = *.cc
DEP_OPTS = -I../common

SOURCE = *.cc *.h

OBJECTS = test_drift.o network_io.o 


all : $(TARGET)

include $(HOME)/lib/make/Makefile.common

include Makefile.dep

MOC = $(QTDIR)/bin/moc
MAKE_PARAMS=$(HOME)/bin/makeParams

install : $(TARGET)
	cp $(TARGET) $(HOME)/bin

install_strip : install
	strip $(HOME)/bin/$(TARGET)

release :
	$(MAKE) DFLAGS="" OFLAGS="-O3 $(ARCH)" && strip $(TARGET)

PBO: clean
	$(MAKE) PROF=-fprofile-generate $(PBO_TARGET)
	$(MAKE) benchmark
	$(MAKE) clean
	$(MAKE) PROF=-fprofile-use $(PBO_TARGET)

new_release : version release

clean :
	rm -f $(OBJECTS) 

all_clean : clean
	rm -f $(TARGET)

benchmark: $(TARGET)
	time ./$(TARGET) $(BENCH_ARGS)
