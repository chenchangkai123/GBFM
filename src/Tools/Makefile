CC=gcc
MPICC=mpicc

FLAG= -g

PAR =  atopkge.o errpkge.o getpars.o subcalls.o syscalls.o

CORE = alloc.o complex.o franuni.o fft.o sinc.o split.o SegyRead.o

segycut :  segycut.o $(CORE) $(PAR) 
	$(CC) $(FLAG) -o segycut  segycut.o $(CORE) $(PAR) 

makemodel : makemodel.o $(CORE) $(PAR) 
	$(CC) $(FLAG) -o makemodel  makemodel.o $(CORE) $(PAR)
	
segycut.o:  segycut.c
	$(CC) $(FLAG) -c   segycut.c

makemodel.o: makemodel.c
	$(CC) $(FLAG) -c   makemodel.c

atopkge.o: atopkge.c
	$(CC) $(FLAG) -c atopkge.c

errpkge.o : errpkge.c
	$(CC) $(FLAG) -c  errpkge.c

subcalls.o : subcalls.c
	$(CC) $(FLAG) -c subcalls.c

syscalls.o : syscalls.c
	$(CC) $(FLAG) -c  syscalls.c

getpars.o : getpars.c
	$(CC) $(FLAG) -c  getpars.c


alloc.o: alloc.c
	$(CC) $(FLAG) -c alloc.c

complex.o: complex.c
	$(CC) $(FLAG) -c complex.c


fft.o: fft.c
	$(CC) $(FLAG) -c fft.c

franuni.o: franuni.c
	$(CC) $(FLAG) -c franuni.c


sinc.o: sinc.c
	$(CC) $(FLAG) -c sinc.c

split.o: split.c
	$(CC) $(FLAG) -c split.c

SegyRead.o: SegyRead.c
	$(CC) $(FLAG) -c SegyRead.c

clean:
	rm -f *.o  segycut makemodel

all: segycut makemodel

install: segycut makemodel
	cp segycut makemodel ../../bin
