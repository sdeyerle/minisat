
#include "SolverGroup.h"
#include "SimpSolver.h"

//Placeholder for real solver thread
void *solver_thread(void *) {
	return NULL;
}

namespace Minisat {
	SolverGroup::SolverGroup(int nthreads) {
		this->nthreads = nthreads;

		//Allocate all the Minisat Solvers into an array
		solvers = new SimpSolver* [nthreads];
		for(int i=0; i<nthreads; i++) {
			solvers[i] = new SimpSolver();
		}
		pthread_mutex_init(&lock, NULL);
		pthread_cond_init(&signal_complete, NULL);
		pthread_rwlock_init(&exportedClauseLock, NULL);
	}

	void SolverGroup::parse_DIMACS(gzFile input_stream) {
		StreamBuffer in(input_stream);

		vec<Lit> lits;
		int vars    = 0;
		int clauses = 0;
		int cnt     = 0;
		for (;;){
			skipWhitespace(in);
			if (*in == EOF) break;
			else if (*in == 'p'){
				if (eagerMatch(in, "p cnf")){
					vars    = parseInt(in);
					clauses = parseInt(in);
				}else{
					printf("PARSE ERROR! Unexpected char: %c\n", *in), exit(3);
				}
			} else if (*in == 'c' || *in == 'p')
				skipLine(in);
			else{
				cnt++;
				this->readClause(in, lits);
				for(int t = 0; t < nthreads;t++)
					solvers[t]->addClause(lits); }
		}
		if (vars != solvers[0]->nVars())
			fprintf(stderr, "WARNING! DIMACS header mismatch: wrong number of variables.\n");
		if (cnt  != clauses)
			fprintf(stderr, "WARNING! DIMACS header mismatch: wrong number of clauses.\n");


	}

	void SolverGroup::readClause(StreamBuffer in, vec<Lit>& lits) {
		int     parsed_lit, var;
		lits.clear();
		for (;;){
			parsed_lit = parseInt(in);
			if (parsed_lit == 0) break;
			var = abs(parsed_lit)-1;
			int initvar = var;
			for(int t = 0; t < nthreads; t++){
				var = initvar;	
				while (var >= solvers[t]->nVars()) solvers[t]->newVar();
			}
			lits.push( (parsed_lit > 0) ? mkLit(var) : ~mkLit(var) );
		}
	}

	SolverGroup::~SolverGroup() {
		for(int i=0; i<nthreads; i++) {
			free(solvers[i]);
		}
		free(solvers);
	}
}


