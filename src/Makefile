OBJFILES = main.o libbpfind.o
CPP = g++
CPPFLAGS = -Wall -static
FORTFLAGS = -lgfortran -lquadmath

OBJDIR = obj
BIN    = bin
SRC    = src
metallic : $(OBJFILES)
	$(CPP) -o $(BIN)/$@ $^ $(CPPFLAGS) $(FORTFLAGS)

$(OBJDIR)/main.o: $(SRC)/main.cpp $(SRC)/atom.h $(SRC)/metalbind.h $(SRC)/metaldefs.h \
      $(SRC)/molecule.h $(SRC)/residue.h $(SRC)/struct.h  $(SRC)/sysparams.h $(SRC)/tree.h $(SRC)/util.h
	$(CPP) -c -o  $(FLAGS) main.cpp
$(OBJDIR)/libbpfind.o: $(SRC)/libbpfind.f
	gfortran -c -o  $(FLAGS) libbpfind.f


.PHONY: clean

clean:
	rm -f *.o

