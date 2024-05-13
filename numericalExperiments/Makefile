include make.inc

default: experiments.exe timeHSEQR.exe subModules

subModules:
	+$(MAKE) -C hqr

experiments.exe: experiments.o
	$(CXX) $(CXXFLAGS) $(INCLUDES) $^ -o $@

timeHSEQR.exe: timeHSEQR.o refFiles/dhseqr.o
	$(FC) $(FFLAGS) $^ $(OPTBLAS) -o $@

timeHSEQR.o: timeHSEQR.c
	$(CC) -c $(CFLAGS) $^ -o $@

experiments.o: experiments.cpp
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) $^ -o $@

clean:
	+$(MAKE) -C hqr clean
	rm *exe *o
