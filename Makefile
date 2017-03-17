all: lapmod

lapmod: lapmod.cpp
	$(CXX) -std=c++14 -O2 lapmod.cpp -o lapmod

.PHONY: clean

clean:
	@rm lapmod
