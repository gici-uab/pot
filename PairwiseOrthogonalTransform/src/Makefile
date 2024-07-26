
#CXXFLAGS=-O0 -g3 -Wall
CXXFLAGS=-O3 -pedantic -Wall -Wstrict-aliasing=1 -Wstrict-null-sentinel -Wno-non-template-friend -Woverloaded-virtual -Wsign-promo -static
# -fno-strict-aliasing -static

all: example

clean:
	rm -f example

example: half.c example.cpp half.h pot.h
	$(CXX) $(CXXFLAGS) half.c example.cpp -o example

