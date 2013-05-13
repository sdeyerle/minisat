
#include "SolverGroup.h"
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

        IntOption    nthreads("PARALLEL", "nthreads",   "Number of concurrent threads", 4, IntRange(1, 32));
        IntOption    share_lim("PARALLEL", "share_lim",   "Maximum learnt clause size that is shared to other threads", 8, IntRange(1, 2000));
	
	bool var_stat_compare (const struct VarCount &i, const struct VarCount &j) {
		int num_neg = i.num_neg;
		int num_pos = i.num_pos;
		double controversy_i =  ((double)num_neg + (double)num_pos)  / ( abs(num_neg - num_pos)/10.0 + 1);
		num_neg = j.num_neg;
		num_pos = j.num_pos;
		double controversy_j = ((double)num_neg + (double)num_pos) / ( abs(num_neg - num_pos)/10.0 + 1);
		return controversy_i < controversy_j;
	}

	void *solver_thread(void *in) { //Mode is MODE_RAND or MODE_GP
		struct sg_thread_status *status = (struct sg_thread_status*) in;
		SimpSolver *solver = status->solver;
	
		printf("Starting solver %d -- ", status->thread_id);
		printPath((char *)"Assumptions", status->assumps);

		if(status->mode == MODE_RAND) {
			solver->setRandomPolarity(true);	
			solver->setRandomFreq( 0.1 ); //Set random freq to 10% after done
			printf("Now settings solver %d to RAND mode\n", status->thread_id);
		} //Presently will be irreversible 

		double start_time = cpuTime();

		//Catch if this is already set and exit early
		solver->exit_now = &status->exit_now;
		solver->thread_id = status->thread_id;

		printf("About to begin solving thread #%d\n", status->thread_id);
		
		lbool ret = solver->solveLimited( *status->assumps );
		//vec<Lit> dummy;
		//lbool ret = solver->solveLimited( dummy );
	
		printf("Thread %d -- Result: %d -- Attempting to lock mutex on thread\n", status->thread_id, toInt(ret));

		//Signal the main thread that we are done
		pthread_mutex_lock(status->lock);
	
		//This must happen after lock is acquired to prevent deadlock
		status->done = true;
		status->result = ret;
	
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
		shared_clause_iters = new unsigned int[nthreads];
		shared_unit_iters = new unsigned int[nthreads];
		guiding_path_queue = NULL;

		
		for(int i=0; i<nthreads; i++) {
			solvers[i] = new SimpSolver();
			thread_status[i] = new struct sg_thread_status;
			solvers[i]->group = this;
			solvers[i]->exportSizeLim = share_lim;
			shared_clause_iters[i] = 0;
			shared_unit_iters[i] = 0;
		}
		//Create joinable threads for portability
		pthread_attr_init(&attr);
		pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

		pthread_mutex_init(&lock, NULL);
		pthread_cond_init(&signal_complete, NULL);
		pthread_rwlock_init(&sharedClauseLock, NULL);
		pthread_rwlock_init(&sharedUnitLock, NULL);

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
					
					struct VarCount dummy;
					dummy.num_pos = 0;
					dummy.num_neg = 0;
					for(int i=0; i<vars; i++) {
						dummy.var = i;
						var_stats.push_back(dummy);
					}
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

		//printVarStats();
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
				while (var >= solvers[t]->nVars()) solvers[t]->newVar(); //Default polarity set here
			}
			lits.push( (parsed_lit > 0) ? mkLit(var) : ~mkLit(var) );
			if(parsed_lit > 0) {
				var_stats[var].num_pos++;
			} else {
				var_stats[var].num_neg++;
			}
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
		std::sort(var_stats.begin(), var_stats.end(), var_stat_compare );
		for(int i=0; i<num_vars; i++) {
			Lit pickedLit = mkLit(var_stats.back().var);
			var_stats.pop_back();
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
		guiding_path_status = new int[num_guiding_paths];
				
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

				pthread_mutex_unlock(&lock);

				for(int j=0; j<nthreads; j++) {
					thread_status[i]->exit_now = true;	
				}
				for(int j=0; j<nthreads; j++) {
					pthread_join(threads[i], 0);
					printf("Thread %d has been joined\n", j);

				}

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
						printf("Waiting for thread %d to exit\n", j);
						pthread_mutex_unlock(&lock);
						pthread_join(threads[j], 0);
						pthread_mutex_lock(&lock);
						printf("Thread %d has been joined\n", j);
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
					printf("thread %d -- exit_now? %s\n", i, thread_status[i]->exit_now ? "true" : "false");
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

	void SolverGroup::exportSharedClause( int thread_id, vec<Lit> &in ) {
		pthread_rwlock_wrlock(&sharedClauseLock);

		vec<Lit> *tmp = new vec<Lit>;
		in.copyTo(*tmp);
		struct sharedClause clause_tmp;
		clause_tmp.clause = tmp;
		clause_tmp.thread_id = thread_id;

		shared_clauses.push_back(clause_tmp);
		//printf("Exported Clause Size: %d\n", (int) shared_clauses.size());

		pthread_rwlock_unlock(&sharedClauseLock);
	}

	bool SolverGroup::getNextSharedUnit( int thread_id, Lit &out ) 
	{
		if( shared_units.size() == 0)
			return false;
		if( shared_unit_iters[thread_id] == shared_units.size()-1 )
			return false;

		pthread_rwlock_rdlock(&sharedUnitLock);

		out = shared_units[ shared_unit_iters[thread_id]++ ]; 

		pthread_rwlock_unlock(&sharedUnitLock);

		printf("Thread: %d -- Importing unit: %s%d\n", thread_id, sign(out) ? "" : "-", var(out) ); 
		return true;
	}

	bool SolverGroup::getNextSharedClause( int thread_id, vec<Lit> &out ) 
	{
		if( shared_clauses.size() == 0)
			return false;
		if( shared_clause_iters[thread_id] == shared_clauses.size()-1 )
			return false;

		pthread_rwlock_rdlock(&sharedClauseLock);
		
		do {
			shared_clause_iters[thread_id]++;
			if( shared_clause_iters[thread_id] == shared_clauses.size()-1 ) {
				pthread_rwlock_unlock(&sharedClauseLock);
				return false;
			}
		} while(  shared_clauses[ shared_clause_iters[thread_id] ].thread_id == thread_id );

		//Copy the shared clause to the output
		(*shared_clauses[ shared_clause_iters[thread_id]++ ].clause).copyTo(out); 

		pthread_rwlock_unlock(&sharedClauseLock);

		//printf("Thread %d -- Importing Shared Clause", thread_id);
		//printPath( (char *)"", &out);
		return true;
	}

	void SolverGroup::exportSharedUnit( Lit &in ) {
		pthread_rwlock_wrlock(&sharedUnitLock);

		shared_units.push_back(in);

		pthread_rwlock_unlock(&sharedUnitLock);
		printf("Exported Unit: %s%d\n",  sign(in) ? "" : "-", var(in) ); 

	}
	
	SolverGroup::~SolverGroup()
	{
		for(int i=0; i<nthreads; i++) {
			free(solvers[i]);
			free(thread_status[i]);
		}
		free(solvers);
		free(thread_status);
		free(thread_guiding_path);
	}

	void SolverGroup::printVarStats() {
		std::sort(var_stats.begin(), var_stats.end(), var_stat_compare );
		for(std::vector<VarCount>::iterator i = var_stats.begin(); i != var_stats.end(); i++) {
			int num_neg = i->num_neg;
			int num_pos = i->num_pos;
			double controversy =  ((double)num_neg + (double)num_pos) / ( abs(num_neg - num_pos)/10.0 + 1);
			printf("%5d -- pos: %4d - neg: %4d - controversy: %8.2f\n", i->var, i->num_pos, i->num_neg, controversy);
		}
	};

	void SolverGroup::resetIters(int thread_id) 
	{
		shared_clause_iters[thread_id] = 0;
		shared_unit_iters[thread_id] = 0;
	}
}



