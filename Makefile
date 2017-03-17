all: lapmod

lapmod: lapmod.cpp
	$(CXX) -std=c++14 -Ofast lapmod.cpp -o lapmod

.PHONY: clean

clean:
	@rm lapmod
