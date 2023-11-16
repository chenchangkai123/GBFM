
all:
	cd src/GaussianBeam && $(MAKE)
	cd src/Tools && $(MAKE)

install:
	cd src/GaussianBeam && $(MAKE) install
	cd src/Tools && $(MAKE) install

clean:
	cd src/GaussianBeam && $(MAKE) clean
	cd src/Tools && $(MAKE) clean
	cd bin && rm *
	


