
OBJS = Gauss.o Solve.o Init.o
HEADER = Gauss.h

CXXFLAGS = -O3 --fast-math

all: a.out

dbgprint: CXXFLAGS += -D__DEBUG__
dbgprint: a.out

prof: CXXFLAGS += -pg
prof: a.out



fdebug: CXXFLAGS = -g -lm -Wall -Wextra -pedantic -Werror
fdebug: a.out

debug: CXXFLAGS = -g -lm -pg -Wall -Wextra -pedantic -Werror
debug : a.out

a.out:	$(OBJS) 
	$(CXX) $(CXXFLAGS) -o $@ $^

%.o:	%.cpp $(HEADER)
	$(CXX) -c $(CFLAGS) $(CXXFLAGS) $(CPPFLAGS) -o $@ $<

#%.o:	$(PROJECT_ROOT)%.c
#	$(CC) -c $(CFLAGS) $(CPPFLAGS) -o $@ $<

clean:
	rm -fr a.out $(OBJS)
