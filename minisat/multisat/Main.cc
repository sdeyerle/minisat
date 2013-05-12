/*****************************************************************************************[Main.cc]
Copyright (c) 2003-2006, Niklas Een, Niklas Sorensson
Copyright (c) 2007,      Niklas Sorensson

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
associated documentation files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge, publish, distribute,
sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or
substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
**************************************************************************************************/

#include <errno.h>
#include <zlib.h>
#include "minisat/utils/System.h"
#include "minisat/utils/ParseUtils.h"
#include "minisat/utils/Options.h"
#include "minisat/core/Dimacs.h"
#include "minisat/multisat/SolverGroup.h"
#include "minisat/multisat/SimpSolver.h"

using namespace Minisat;

//Note: All the graceful exit stuff has been removed because the ManySAT guys never bothered
//to actually make it work and I'm sick of the warnings

//=================================================================================================
// Main:

int main(int argc, char** argv)
{
    try {
        setUsageHelp("USAGE: %s [options] <input-file> <result-output-file>\n\n  where input may be either in plain or gzipped DIMACS.\n");
        setX86FPUPrecision();
        
        // Extra options:
        //
        IntOption    verb    ("MAIN", "verb",   "Verbosity level (0=silent, 1=some, 2=more).", 1, IntRange(0, 2));
        IntOption    nthreads("MAIN", "nthreads",   "Number of concurrent threads", 4, IntRange(1, 32));
        BoolOption   pre     ("MAIN", "pre",    "Completely turn on/off any preprocessing.", true);
        BoolOption   solve   ("MAIN", "solve",  "Completely turn on/off solving after preprocessing.", true);
        StringOption dimacs  ("MAIN", "dimacs", "If given, stop after preprocessing and write the result to this file.");
        IntOption    cpu_lim ("MAIN", "cpu-lim","Limit on CPU time allowed in seconds.\n", 0, IntRange(0, INT32_MAX));
        IntOption    mem_lim ("MAIN", "mem-lim","Limit on memory usage in megabytes.\n", 0, IntRange(0, INT32_MAX));
        BoolOption   strictp ("MAIN", "strict", "Validate DIMACS header during parsing.", false);

        parseOptions(argc, argv, true);


	SolverGroup group(nthreads);  //Initialize a parallel solver with nthreads
	gzFile in = (argc == 1) ? gzdopen(0, "rb") : gzopen(argv[1], "rb");
        if (in == NULL)
            printf("ERROR! Could not open file: %s\n", argc == 1 ? "<stdin>" : argv[1]), exit(1);
       
        group.parse_DIMACS(in);
        gzclose(in);

	for(int i=0; i<nthreads; i++) {
		group.solvers[i]->verbosity = verb;
	}

	group.eliminate_parallel(true);

	if( !group.solvers[0]->okay() ) {
		printf("-- UNSAT --  Solved by simplification\n");
		exit(0);
	}

	int winning_thread;
	lbool ret = group.solve_parallel(&winning_thread);
	printf("Winning thread is: %d\n", winning_thread);
	printf(ret == l_True ? "SATISFIABLE\n" : ret == l_False ? "UNSATISFIABLE\n" : "INDETERMINATE\n");
	if(winning_thread > 0)
		group.solvers[winning_thread]->printStats();

        
 //       SimpSolver  S;
 //       double      initial_time = cpuTime();

 //       if (!pre) S.eliminate(true);

 //       S.verbosity = verb;
 //       
 //       solver = &S;
 //       // Use signal handlers that forcibly quit until the solver will be able to respond to
 //       // interrupts:
 //       sigTerm(SIGINT_exit);

 //       // Try to set resource limits:
 //       if (cpu_lim != 0) limitTime(cpu_lim);
 //       if (mem_lim != 0) limitMemory(mem_lim);

 //       if (argc == 1)
 //           printf("Reading from standard input... Use '--help' for help.\n");

 //       gzFile in = (argc == 1) ? gzdopen(0, "rb") : gzopen(argv[1], "rb");
 //       if (in == NULL)
 //           printf("ERROR! Could not open file: %s\n", argc == 1 ? "<stdin>" : argv[1]), exit(1);
 //       
 //       if (S.verbosity > 0){
 //           printf("============================[ Problem Statistics ]=============================\n");
 //           printf("|                                                                             |\n"); }
 //       
 //       parse_DIMACS(in, S, (bool)strictp);
 //       gzclose(in);
 //       FILE* res = (argc >= 3) ? fopen(argv[2], "wb") : NULL;

 //       if (S.verbosity > 0){
 //           printf("|  Number of variables:  %12d                                         |\n", S.nVars());
 //           printf("|  Number of clauses:    %12d                                         |\n", S.nClauses()); }
 //       
 //       double parsed_time = cpuTime();
 //       if (S.verbosity > 0)
 //           printf("|  Parse time:           %12.2f s                                       |\n", parsed_time - initial_time);

 //       // Change to signal-handlers that will only notify the solver and allow it to terminate
 //       // voluntarily:
 //       sigTerm(SIGINT_interrupt);

 //       S.eliminate(true);
 //       double simplified_time = cpuTime();
 //       if (S.verbosity > 0){
 //           printf("|  Simplification time:  %12.2f s                                       |\n", simplified_time - parsed_time);
 //           printf("|                                                                             |\n"); }


 //       //The solver is now ready to run, so let's duplicate it once for each thread.
 //       SimpSolver solver;
 //       SimpSolver solver2 = SimpSolver(solver);

 //       if (!S.okay()){
 //           if (res != NULL) fprintf(res, "UNSAT\n"), fclose(res);
 //           if (S.verbosity > 0){
 //               printf("===============================================================================\n");
 //               printf("Solved by simplification\n");
 //               S.printStats();
 //               printf("\n"); }
 //           printf("UNSATISFIABLE\n");
 //           exit(20);
 //       }

 //       lbool ret = l_Undef;

 //       if (solve){
 //           vec<Lit> dummy;
 //           ret = S.solveLimited(dummy);
 //       }else if (S.verbosity > 0)
 //           printf("===============================================================================\n");

 //       if (dimacs && ret == l_Undef)
 //           S.toDimacs((const char*)dimacs);

 //       if (S.verbosity > 0){
 //           S.printStats();
 //           printf("\n"); }
 //       printf(ret == l_True ? "SATISFIABLE\n" : ret == l_False ? "UNSATISFIABLE\n" : "INDETERMINATE\n");
 //       if (res != NULL){
 //           if (ret == l_True){
 //               fprintf(res, "SAT\n");
 //               for (int i = 0; i < S.nVars(); i++)
 //                   if (S.model[i] != l_Undef)
 //                       fprintf(res, "%s%s%d", (i==0)?"":" ", (S.model[i]==l_True)?"":"-", i+1);
 //               fprintf(res, " 0\n");
 //           }else if (ret == l_False)
 //               fprintf(res, "UNSAT\n");
 //           else
 //               fprintf(res, "INDET\n");
 //           fclose(res);
 //       }

 //#ifdef NDEBUG
 //       exit(ret == l_True ? 10 : ret == l_False ? 20 : 0);     // (faster than "return", which will invoke the destructor for 'Solver')
 //#else
 //       return (ret == l_True ? 10 : ret == l_False ? 20 : 0);
 //#endif
    } catch (OutOfMemoryException&){
        printf("===============================================================================\n");
        printf("INDETERMINATE\n");
        exit(0);
    }
}
