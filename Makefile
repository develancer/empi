# University of Warsaw, Department of Biomedical Physics
# See LICENCE for details.
CXXFLAGS = -std=c++11 -fopenmp -Iinc -Wall -Wextra -Wpedantic -O2
LDFLAGS = -lfftw3

.PHONY : default

default : empi

clean :
	@rm -frv empi obj

empi : src/empi.cpp obj/classes.o obj/conf.o obj/gabor.o obj/io.o obj/mmp.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

obj/%.o : src/%.cpp | obj
	$(CXX) $(CXXFLAGS) -c -o $@ $^

obj :
	mkdir -pv obj
