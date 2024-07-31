#PROJECT_ROOT = $(dir $(abspath $(lastword $(MAKEFILE_LIST))))

OBJS = Gauss.o Solve.o Init.o #ThreadComp.o
HEADER = Gauss.h
CXXFLAGS = -O3 -pthread --fast-math


all: a.out

RESO: CXXFLAGS += -D__RESO__
RESO: a.out

dbgprint: CXXFLAGS += -D__DEBUG__
dbgprint: a.out

dbgcomp: CXXFLAGS += -g -lm -pg -Wall -Wextra -pedantic -Werror
dbgcomp: a.out

fulldbg: CXXFLAGS += -g -lm -pg -Wall -Wextra -pedantic -Werror -D__DEBUG__
fulldbg: a.out

a.out:	$(OBJS) 
	$(CXX) -o $(CXXFLAGS) -o $@ $^

%.o:	%.cpp $(HEADER)
	$(CXX) -c $(CXXFLAGS) -o $@ $<

#%.o:	$(PROJECT_ROOT)%.c
#	$(CC) -c $(CFLAGS) $(CPPFLAGS) -o $@ $<

clean:
	rm -fr a.out $(OBJS)
