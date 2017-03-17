all: lapmod

lapmod: lapmod.cpp
	$(CXX) -std=c++14 -march=native -Ofast -fno-rtti -Wall -Wextra -Werror -pedantic -Wshadow -Wmissing-include-dirs -Winvalid-pch -Wformat=2 lapmod.cpp -o lapmod

.PHONY: clean

clean:
	-@$(RM) lapmod 2> /dev/null
