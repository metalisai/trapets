CXX=g++
CXXFLAGS=-std=c++11 -Wall -Wextra -ggdb
LDFLAGS=
OBJ=main.o

trapets: $(OBJ)
	$(CXX) $(LDFLAGS) -o $@ $^

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@
