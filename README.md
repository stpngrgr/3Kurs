# 3Kurs
Code I've written during my third year in university.

I was not in a habit of adding comments back then, so the code can be a bit hard to follow. I've been meaning to come back and review it to make it a little more readable.

The folders whose names start with "Gauss" are programms for solving a system of linear equations.
- "Gauss" is a straigtforward solution.
- "GaussOpt" is slightly optimised by unrolling the loops to make better use of instruction pipelining.
- "GaussThreads" is a multithreaded implementation.
- "GaussMPI" is an MPI implementation.
- "Eigenvalues" is a programm for approximating eigenvalues of a matrix.
- "Approximation" has two programms both for consstructing polynomial interpolation, one for interpolation of unary functions and a multithreaded one for binary functions. Both programms are visualised in Qt5 and so should be built using qmake.
