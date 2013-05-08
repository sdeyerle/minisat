
#include "SolverGroup.h"
#include "SimpSolver.h"
#include "minisat/utils/System.h"
#include <signal.h>
#include <omp.h>
#include <math.h>
#include <vector>


namespace Minisat {

	//Debugger helper function
	void printPath( char *prefix, vec<Lit> *in) {
		printf("%s: ", prefix);
		for(int j=0; j<in->size(); j++) {
			printf("%s%d ", sign( (*in)[j] ) ? "-" : "", var( (*in)[j] ) );
		}
		printf("\n");
	}

	void *solver_thread(void *in) { //Mode is MODE_RAND or MODE_GP
		struct sg_thread_status *status = (struct sg_thread_status*) in;
		SimpSolver *solver = status->solver;
	
		printf("Starting a solver -- ");
		printPath((char *)"Assumptions", status->assumps);

		if(status->mode == MODE_RAND) {
			solver->setRandomPolarity(true);	
			solver->setRandomFreq( 0.1 ); //Set random freq to 10% after done
		} //Presently will be irreversible 

		double start_time = cpuTime();

		printf("About to begin solving thread #%d\n", status->thread_id);
		
		lbool ret = solver->solveLimited( *status->assumps );
		//vec<Lit> dummy;
		//lbool ret = solver->solveLimited( dummy );
		status->done = true;
		status->result = ret;

		//Signal the main thread that we are done
		pthread_mutex_lock(status->lock);
		
		pthread_cond_signal(status->signal_complete);

		pthread_mutex_unlock(status->lock);
		
		printf("Finished solving thread #%d.  Time: %f s\n", status->thread_id, cpuTime() - start_time);

		pthread_exit(NULL);
	}

	SolverGroup::SolverGroup(int nthreads) {
		this->nthreads = nthreads;

		//Allocate all the Minisat Solvers into an array
		solvers = new SimpSolver* [nthreads];
		thread_status = new struct sg_thread_status* [nthreads];
		threads = new pthread_t [nthreads];
		thread_guiding_path = new int[nthreads];
		guiding_path_status = new int[num_guiding_paths];
		guiding_path_queue = NULL;

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
		double start_time = cpuTime();
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

		printf("Parsing time: %f s\n", cpuTime() - start_time);

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

		if(nthreads==1) {
			num_guiding_paths = 0;
			vec< vec<Lit> *> *guiding_path_queue = new vec< vec<Lit>* >;
			vec<Lit> *dummy = new vec<Lit>;
			guiding_path_queue->push(dummy);
			return guiding_path_queue;
		}
		
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
				if( (i>>j) % 2 ) { //Create every combination of inversions
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
		guiding_path_queue = createGuidingPaths(nthreads);
		num_guiding_paths = guiding_path_queue->size();
				
		//Initialize thread/guiding path relation to indeterminate at start
		for(int i=0; i<nthreads; i++) {
			thread_guiding_path[i] = -1;;
		}
		for(int i=0; i<num_guiding_paths; i++) {
			guiding_path_status[i] = INDETERM;
		}

		//Lock the primary mutex for the pthread_cond_t
		pthread_mutex_lock(&lock);
		
		current_guiding_path = 0;
		for(int i=0; i<nthreads; i++) {
			thread_status[i]->guiding_path_num = current_guiding_path;
			thread_guiding_path[i] = current_guiding_path;
			thread_status[i]->assumps = (*guiding_path_queue)[current_guiding_path++];
			thread_status[i]->signal_complete = &signal_complete;
			thread_status[i]->lock = &lock;
			thread_status[i]->thread_id = i;
			thread_status[i]->done = false;
			thread_status[i]->exit_now = false;
			thread_status[i]->mode = MODE_GP;
			thread_status[i]->solver = solvers[i];
			int ret = pthread_create(&threads[i], &attr, solver_thread, thread_status[i]); 
			printf("Created thread %d with return value %d\n", i, ret);
		}

		while(1) {
			pthread_cond_wait(&signal_complete, &lock);
			printf("Signal caught.....\n");
			lbool status;
			bool done = processCompleteSolvers(status, winning_thread);
			if(done)
				return status;
		}
	}

	int SolverGroup::findOldThread() {
		for(int i=0; i<num_guiding_paths; i++) {
			if(guiding_path_status[i] == INDETERM)
				return i;
		}
		return -1;
	}

	bool SolverGroup::processCompleteSolvers(lbool &status, int *winning_thread) {
		printf("processorCompleteSolvers() called\n");
		for(int i=0; i<nthreads; i++) {
			if(thread_status[i]->done == true && thread_status[i]->result == l_True) {
				printf("Thread:%d was SAT\n", i);	

				for(int j=0; j<nthreads; j++) {
					thread_status[i]->exit_now = true;	
				}
				for(int j=0; j<nthreads; j++) {
					pthread_join(threads[i], 0);
				}

				pthread_mutex_unlock(&lock);
				guiding_path_status[i] = SAT;
				*winning_thread = i;
				status = thread_status[i]->result; //Return which thread finished first
				return true;
			}
			else if(thread_status[i]->done == true) {
				guiding_path_status[ thread_status[i]->guiding_path_num ] = UNSAT;

				for(int j=0; j<nthreads; j++) {
					if( thread_status[j]->guiding_path_num ==  thread_status[i]->guiding_path_num ) {
						thread_status[j]->exit_now = true;
					}
				}

				printf("Guiding path %d is complete\n", thread_status[i]->guiding_path_num);

				if(current_guiding_path != num_guiding_paths) { //Still working on guiding path case
					printf("Launching guiding_path %d on thread %d\n", current_guiding_path, i);
					thread_guiding_path[i] = current_guiding_path;
					thread_status[i]->guiding_path_num = current_guiding_path;
					thread_status[i]->assumps = (*guiding_path_queue)[current_guiding_path++];
					thread_status[i]->done = false;
					pthread_create(&threads[i], &attr, solver_thread, thread_status[i]);
				}
				else { //Random solver case
					int oldest_gp = findOldThread();

					if(oldest_gp < 0) //If findOldThread returns 0, every path has already been solved
						break;

					printf("Launching random solver for guiding_path %d on thread %d\n", oldest_gp, i);
					thread_guiding_path[i] = oldest_gp;
					thread_status[i]->guiding_path_num = oldest_gp;
					thread_status[i]->assumps = (*guiding_path_queue)[oldest_gp];
					thread_status[i]->done = false;
					thread_status[i]->exit_now = false;
					thread_status[i]->mode = MODE_RAND;
					pthread_create(&threads[i], &attr, solver_thread, thread_status[i]);
				}
			}
		}

		int num_unsat = 0;
		for(int i=0; i<nthreads; i++) {
			if(guiding_path_status[i] == UNSAT)
				num_unsat++;
		}
		if(num_unsat == num_guiding_paths) {
			*winning_thread = -1;
			status = l_False;
			return true;
		}
		num_unsat = 0;
		return false;
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
		free(thread_guiding_path);
	}
}


