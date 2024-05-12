include make.inc

experiments.exe: experiments.o
	$(CXX) $(CXXFLAGS) $(INCLUDES) $^ -o $@

experiments.o: experiments.cpp
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) $^ -o $@

clean:
	rm *exe *o
