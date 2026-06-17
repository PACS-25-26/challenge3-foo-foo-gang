# change this with your ROOT
PACS_ROOT = /home/alessandro-romani/Desktop/PACS/pacs-examples/Examples

# Eigen library dependencies
EIGEN_DIR= $(PACS_ROOT)/eigen3		 # pacs-examples directory
EIGEN_SYS_DIR = /usr/include/eigen3  # local one

CXX = mpicxx # compiler for parallel code -Wall -Wextra 
CXXFLAGS = -std=c++20 -O3 -pedantic -fopenmp \
		   -Iinclude \
		   -I$(EIGEN_DIR) \
		   -I$(EIGEN_SYS_DIR)

LDLIBS = -lmuparser

# get all files *.cpp in /src folder
SRCS=$(wildcard src/*.cpp)
# get the corresponding object file
OBJS = $(SRCS:.cpp=.o)
# get all headers in the /include folder
HEADERS=$(wildcard include/*.hpp)

# create executable
exe_sources=$(filter src/main%.cpp, $(SRCS))
EXEC=$(notdir $(exe_sources:.cpp=))

# object definitions
.PHONY= all clean distclean

.DEFAULT_GOAL = all

all: $(EXEC)

# compilation rules
src/%.o: src/%.cpp $(HEADERS)
		$(CXX) $(CXXFLAGS) -c $< -o $@

$(EXEC): $(OBJS)
		$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS) $(LDLIBS)

# cleaning rules
clean:
	$(RM) -f $(EXEC) $(OBJS)

# remove also vtk solution file
distclean: clean
	$(RM) *.out *.vtk *~
