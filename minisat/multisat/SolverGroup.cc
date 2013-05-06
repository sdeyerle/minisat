
#include "SolverGroup.h"
#include "SimpSolver.h"
#include <signal.h>
#include <omp.h>
#include <math.h>
#include <vector>

namespace Minisat {
	//Placeholder for real solver thread
	void *solver_thread(void *in) {
		struct sg_thread_status *status = (struct sg_thread_status*) in;
		SimpSolver *solver = status->solver;
		status->trail = status->solver->getTrail();
		
		printf("About to begin solving thread #%d\n", status->thread_id);

		if(*status->group_done) {
			printf("Thread #%d aborting. 1 \n", status->thread_id);
			pthread_exit(NULL);
		}


		lbool ret = solver->solveLimited( *status->assumps );
		status->done = true;
		status->result = ret;

		if(*status->group_done) {
			printf("Thread #%d aborting. 2 \n", status->thread_id);
			pthread_exit(NULL);
		}

		//Signal the main thread that we are done
		pthread_mutex_lock(status->lock);
		
		pthread_cond_signal(status->signal_complete);

		pthread_mutex_unlock(status->lock);
		
		printf("Finished solving thread #%d\n", status->thread_id);

		pthread_exit(NULL);
	}

	SolverGroup::SolverGroup(int nthreads) {
		this->nthreads = nthreads;

		//Allocate all the Minisat Solvers into an array
		solvers = new SimpSolver* [nthreads];
		thread_status = new struct sg_thread_status* [nthreads];
		threads = new pthread_t [nthreads];

		for(int i=0; i<nthreads; i++) {
			solvers[i] = new SimpSolver();
			thread_status[i] = new struct sg_thread_status;
		}
		//Create joinable threads for portability
		pthread_attr_init(&attr);
		pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

		pthread_mutex_init(&lock, NULL);
		pthread_cond_init(&signal_complete, NULL);
		pthread_rwlock_init(&exportedClauseLock, NULL);

		omp_set_num_threads(nthreads);
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
					printf("Clauses: %d Vars: %d\n", clauses, vars);
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

	void SolverGroup::readClause(StreamBuffer &in, vec<Lit>& lits) {
		int     parsed_lit, var;
		lits.clear();
		for (;;){
			parsed_lit = parseInt(in);
			if (parsed_lit == 0) break;
			var = abs(parsed_lit)-1;
			int initvar = var;
			for(int t = 0; t < nthreads; t++){
				var = initvar;	
				while (var >= solvers[t]->nVars()) solvers[t]->newVar(l_False);
			}
			lits.push( (parsed_lit > 0) ? mkLit(var) : ~mkLit(var) );
		}
	}

	//The input number of paths which is rounded up to 2^n
	vec< vec<Lit>* > *SolverGroup::createGuidingPaths(int num) {
		int num_guiding_paths = pow(2, ceil(log2(num)));
		printf("Num Guiding Paths: %d\n", num_guiding_paths);

		int num_vars = log2(num_guiding_paths); //The number of variables in the guiding paths

		typedef vec<Lit> simp_clause;
		vec<Lit> guiding_path;
		vec<simp_clause*> *guiding_path_queue = new vec< vec<Lit>* >;

		printf("Picking guiding path vars \n");
		for(int i=0; i<num_vars; i++) {
			Lit pickedLit =  solvers[0]->pickGuidingPathLit();
			guiding_path.push( pickedLit );
			printf("%d ", var(pickedLit) );
		}
		printf("\n");
		for(int i=0; i<num_guiding_paths; i++) {
			printf("Guiding Path %d:", i);
			vec<Lit> *new_gp = new vec<Lit>;
			for(int j=0; j<num_vars; j++) {
				Lit temp;
				if( (i>>j) % 2 ) {
					temp = guiding_path[j];
				}
				else {
					temp = ~guiding_path[j];
				}
				new_gp->push(temp);
			}
			for(int j=0; j<num_vars; j++) {
				printf("%s%d ", sign( (*new_gp)[j] ) ? "-" : "", var( (*new_gp)[j] ) );
			}
			printf("\n");
			guiding_path_queue->push(new_gp);
		}
		return guiding_path_queue;

	}

	lbool SolverGroup::solve_parallel(int *winning_thread) {
		bool group_done = false;

		enum stat {
			INDETERM,
			SAT,
			UNSAT
		};

		vec< vec<Lit>* > *guiding_path_queue = createGuidingPaths(nthreads);
		int num_paths = guiding_path_queue->size();
		int *guiding_path_status = new int[num_paths];
	
		for(int i=0; i<num_paths; i++) {
			guiding_path_status[i] = INDETERM;
		}

		//Lock the primary mutex for the pthread_cond_t
		pthread_mutex_lock(&lock);

		for(int i=0; i<nthreads; i++) {
			thread_status[i]->assumps = guiding_path_queue->last();
			guiding_path_queue->pop();
			thread_status[i]->signal_complete = &signal_complete;
			thread_status[i]->lock = &lock;
			thread_status[i]->thread_id = i;
			thread_status[i]->done = false;
			thread_status[i]->solver = solvers[i];
			thread_status[i]->group_done = &group_done;
			int ret = pthread_create(&threads[i], &attr, solver_thread, thread_status[i]); 
			printf("Created thread %d with return value %d\n", i, ret);
		}

		while(1) {
			pthread_cond_wait(&signal_complete, &lock);
			for(int i=0; i<nthreads; i++) {
				if(thread_status[i]->done == true && thread_status[i]->result == l_True) {
					printf("Thread:%d was SAT\n", i);	
					group_done = true; //Tell all threads to exit
					
					for(int j=0; j<nthreads; j++) {
						pthread_kill(threads[i], 9);
					}
					
					pthread_mutex_unlock(&lock);
					guiding_path_status[i] = SAT;
					*winning_thread = i;
					return thread_status[i]->result; //Return which thread finished first
				}
				else if(thread_status[i]->done == true && thread_status[i]->result == l_False) {
					printf("Thread:%d was UNSAT\n", i);	
					guiding_path_status[i] = UNSAT;
				}
			}
			printf("Thread 0 trail: ");
			for(int i=0; i< thread_status[0]->trail->size(); i++) {
				printf("%s%d ", sign( (*thread_status[0]->trail)[i] ) ? "-" : "", var( (*thread_status[0]->trail)[i] ));
			}
			printf("\n");
			printf("Thread 0 trail decisions: ");
			for(int i=0; i<  solvers[i]->trail_lim.size(); i++) {
				printf("%d ",  solvers[i]->trail_lim[i] ); 
			}
			printf("\n");

			int num_unsat = 0;
			for(int i=0; i<nthreads; i++) {
				if(guiding_path_status[i] == UNSAT)
					num_unsat++;
			}
			if(num_unsat == num_paths) {
				return l_False;
			}
			num_unsat = 0;
		}
	}
	
	bool SolverGroup::eliminate_parallel(bool in) {
#pragma omp parallel for
		for (int i=0; i<nthreads; i++) {
			solvers[i]->eliminate(in);
			printf("Simplified thread #%d\n", i);
		}
		
		return true;
	}

	SolverGroup::~SolverGroup() {
		for(int i=0; i<nthreads; i++) {
			free(solvers[i]);
			free(thread_status[i]);
		}
		free(solvers);
		free(thread_status);
	}
}


