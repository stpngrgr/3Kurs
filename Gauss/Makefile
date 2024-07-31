#PROJECT_ROOT = $(dir $(abspath $(lastword $(MAKEFILE_LIST))))

OBJS = Gauss.o Solve.o Init.o
HEADER = Gauss.h

CXXFLAGS += -g -lm -o3 -Wall -Wextra -pedantic -Werror

all:	a.out

opt: CFLAGS +=
opt : a.out

a.out:	$(OBJS) 
	$(CXX) -pg -o $@ $^

%.o:	%.cpp $(HEADER)
	$(CXX) -c $(CFLAGS) $(CXXFLAGS) $(CPPFLAGS) -o $@ $<

#%.o:	$(PROJECT_ROOT)%.c
#	$(CC) -c $(CFLAGS) $(CPPFLAGS) -o $@ $<

clean:
	rm -fr a.out $(OBJS)
