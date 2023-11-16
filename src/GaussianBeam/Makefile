
GCC=gcc

MPICC=mpicc

FLAGS = -ld_classic -lm #if not macOS, remove -ld_classic

PAR =  atopkge.o errpkge.o getpars.o subcalls.o syscalls.o

CORE = alloc.o complex.o franuni.o fft.o sinc.o

migration : migration.o $(CORE) $(PAR)
	$(MPICC) -o migration migration.o $(CORE) $(PAR) $(FLAGS) 
forward : forward.o $(CORE) $(PAR)
	$(GCC) -o forward forward.o $(CORE) $(PAR) $(FLAGS)
makehead : makehead.o $(CORE) $(PAR)
	$(GCC) -o makehead makehead.o $(CORE) $(PAR) $(FLAGS)
lsm : lsm.o $(CORE) $(PAR)
	$(GCC) -o lsm lsm.o $(CORE) $(PAR) $(FLAGS)
subtract : subtract.o $(CORE) $(PAR)
	$(GCC) -o subtract subtract.o $(CORE) $(PAR) $(FLAGS)
generate : generate.o $(CORE) $(PAR)
	$(GCC) -o generate generate.o $(CORE) $(PAR) $(FLAGS)
migrationGBS : migrationGBS.o $(CORE) $(PAR)
	$(MPICC) -o migrationGBS migrationGBS.o $(CORE) $(PAR) $(FLAGS)

forward.o : forward.c
	$(GCC) -c forward.c
migration.o : migration.c
	$(MPICC) -c migration.c
migrationGBS.o : migrationGBS.c
	$(MPICC) -c migrationGBS.c
migrationCBM.o : migrationCBM.c
	$(MPICC) -c migrationCBM.c
# core
alloc.o : alloc.c
	$(GCC) -c alloc.c
complex.o : complex.c
	$(GCC) -c complex.c
franuni.o : franuni.c
	$(GCC) -c franuni.c
fft.o : fft.c
	$(GCC) -c fft.c
sinc.o : sinc.c
	$(GCC) -c sinc.c

# par
atopkge.o : atopkge.c
	$(GCC) -c atopkge.c
errpkge.o : errpkge.c
	$(GCC) -c errpkge.c
getpars.o : getpars.c
	$(GCC) -c getpars.c
subcalls.o : subcalls.c
	$(GCC) -c subcalls.c
syscalls.o : syscalls.c
	$(GCC) -c syscalls.c

# head
makehead.o : makehead.c
	$(GCC) -c makehead.c

# lsm
subtract.o : subtract.c
	$(GCC) -c subtract.c

lsm.o : lsm.c
	$(GCC) -c lsm.c

# tools
generate.o : generate.c
	$(GCC) -c generate.c


clean:
	rm forward migration makehead subtract lsm generate migrationGBS *.o

all: forward migration makehead subtract lsm generate migrationGBS

install: forward migration makehead subtract lsm generate migrationGBS
	cp forward migration makehead subtract lsm generate migrationGBS ../../bin
