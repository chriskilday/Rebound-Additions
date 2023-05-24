export OPENGL=0
export MPI=1
export CC=mpicc

# Path to /rebound/src/Makefile.defs
include ../src/Makefile.defs

all: librebound
	@echo ""
	@echo "Compiling problem file ..."
	$(CC) -I../src/ -Wl,-rpath,./ $(OPT) $(PREDEF) polar_binary.c output.c -L. -lrebound $(LIB) -g -o pb.exe
	@echo ""
	@echo "REBOUND compiled successfully."

librebound:
	@echo "Compiling shared library librebound.so ..."
	$(MAKE) -C ../src/
	@-rm -f librebound.so
	@ln -s ../src/librebound.so .

clean:
	@echo "Cleaning up shared library librebound.so ..."
	@-rm -f librebound.so
	$(MAKE) -C ../src/ clean
	@echo "Cleaning up local directory ..."
	@-rm -vf pb.exe

