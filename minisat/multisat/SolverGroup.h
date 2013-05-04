#ifndef Minisat_SolverGroup_H
#define Minisat_SolverGroup_H

#include <vector>
#include "SimpSolver.h"
#include <pthread.h>

namespace Minisat {

struct thread_status {
	vec<Lit>	assumps;
	vec<Lit>	trail;
	bool		done;
	pthread_cond_t	*signal_complete;
	pthread_mutex_t *lock;
	int 		thread_id;
};

class SolverGroup {
public:
	SolverGroup(int nthreads);
	SolverGroup();
	~SolverGroup();

	void parse_DIMACS(gzFile in);
	int solve_parallel();
	void exportClause(vec<Lit> clause);
private:
	void readClause(StreamBuffer in, vec<Lit>& lits);

public:
	SimpSolver ** solvers;
private:
	vec<SharedClause> exportClauses;
	pthread_attr_t attr;
	pthread_t *threads;
	pthread_cond_t signal_complete; //Thread completion
	pthread_mutex_t lock;  //Main mutex for threads
	pthread_rwlock_t exportedClauseLock;
	int nthreads;


}; //class SolverGroup

} //namespace Minisat

#endif //Minisat_SolverGroup_h
