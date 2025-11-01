CXX=g++
CXXFLAGS=-std=c++17 -Wall -Wextra -ggdb
LDFLAGS=
OBJ=main.o

trapets: $(OBJ)
	$(CXX) $(LDFLAGS) -o $@ $^

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@
