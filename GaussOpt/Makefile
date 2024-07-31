#PROJECT_ROOT = $(dir $(abspath $(lastword $(MAKEFILE_LIST))))

OBJS = Gauss.o Solve.o Init.o
HEADER = Gauss.h

CXXFLAGS = -o3 --fast-math

all:	a.out

fdebug: CXXFLAGS = -g -lm -o3 -Wall -Wextra -pedantic -Werror
fdebug: a.out

debug: CXXFLAGS = -g -lm -pg -Wall -Wextra -pedantic -Werror
debug : a.out

a.out:	$(OBJS) 
	$(CXX) -pg -o $@ $^

%.o:	%.cpp $(HEADER)
	$(CXX) -c $(CFLAGS) $(CXXFLAGS) $(CPPFLAGS) -o $@ $<

#%.o:	$(PROJECT_ROOT)%.c
#	$(CC) -c $(CFLAGS) $(CPPFLAGS) -o $@ $<

clean:
	rm -fr a.out $(OBJS)
