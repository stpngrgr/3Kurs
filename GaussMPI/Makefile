#PROJECT_ROOT = $(dir $(abspath $(lastword $(MAKEFILE_LIST))))

MPI = mpicxx

OBJS = Gauss.o Solve.o Init.o #ThreadComp.o
HEADER = Gauss.h
CXXFLAGS = -O3 -pthread --fast-math


all: a.out


fullres: CXXFLAGS += -D__FULLRES__
fullres: a.out

reso: CXXFLAGS += -D__RESO__
reso: a.out

dbgprint: CXXFLAGS += -D__DEBUG__
dbgprint: a.out

fdebug: CXXFLAGS = -g -lm -Wall -Wextra -pedantic -Werror -D__DEBUG__
fdebug: a.out

debug: CXXFLAGS = -g -lm -pg -Wall -Wextra -pedantic -D__DEBUG__
debug : a.out

warning: CXXFLAGS = -g -lm -pg -Wall -Wextra -pedantic -Werror -D__DEBUG__
warning : a.out

a.out:	$(OBJS) 
	$(MPI) -o $(CXXFLAGS) -o $@ $^

%.o:	%.cpp $(HEADER)
	$(MPI) -c $(CXXFLAGS) -o $@ $<

#%.o:	$(PROJECT_ROOT)%.c
#	$(CC) -c $(CFLAGS) $(CPPFLAGS) -o $@ $<

clean:
	rm -fr a.out $(OBJS)
