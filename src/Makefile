
FC = g++
FCLAGS = -O3 -Wall -shared -std=c++11 -undefined dynamic_lookup 
PBFLAGS = $$(python3.9 -m pybind11 --includes)
CFFLAGS = $$(python3.9-config --extension-suffix)


SOURCE = functions.cpp pybind.cpp

all: abp

abp: $(SOURCE)
	$(FC) $(FCLAGS) $(PBFLAGS) $(SOURCE) -o $@$(CFFLAGS)

clean:
	-rm -f *.txt *.png abp.*