#PROJECT_ROOT = $(dir $(abspath $(lastword $(MAKEFILE_LIST))))

OBJS = Eigenvalues.o Solve.o Init.o #ThreadComp.o
HEADER = Eigenvalues.h
CXXFLAGS = -O3 --fast-math


all: a.out

dbgprint: CXXFLAGS += -D__DEBUG__
dbgprint: a.out

fdebug: CXXFLAGS = -g -lm -Wall -Wextra -pedantic -Werror -D__DEBUG__
fdebug: a.out

debug: CXXFLAGS = -g -lm -pg -Wall -Wextra -pedantic -Werror -D__DEBUG__
debug : a.out

a.out:	$(OBJS) 
	$(CXX) -o $(CXXFLAGS) -o $@ $^

%.o:	%.cpp $(HEADER)
	$(CXX) -c $(CXXFLAGS) -o $@ $<

#%.o:	$(PROJECT_ROOT)%.c
#	$(CC) -c $(CFLAGS) $(CPPFLAGS) -o $@ $<

clean:
	rm -fr a.out $(OBJS)
