#ifndef Minisat_SolverGroup_H
#define Minisat_SolverGroup_H

#include <vector>
#include "SimpSolver.h"
#include <pthread.h>

namespace Minisat {

struct sg_thread_status {
	vec<Lit>	*assumps;
	vec<Lit>	*trail;
	bool		done;
	lbool		result;
	pthread_cond_t	*signal_complete;
	pthread_mutex_t *lock;
	int 		thread_id;
	SimpSolver 	*solver;
	bool		*group_done;
};

class SolverGroup {
public:
	SolverGroup(int nthreads);
	SolverGroup();
	~SolverGroup();

	void parse_DIMACS(gzFile in);
	lbool solve_parallel(int *winning_thread);
	bool eliminate_parallel(bool in);
	void exportClause(vec<Lit> clause);
	vec< vec<Lit>* > *createGuidingPaths(int num);
private:
	void readClause(StreamBuffer &in, vec<Lit>& lits);

public:
	SimpSolver ** solvers;
private:
	vec<SharedClause> exportClauses;
	pthread_attr_t attr;
	struct sg_thread_status **thread_status;
	pthread_t *threads;
	pthread_cond_t signal_complete; //Thread completion
	pthread_mutex_t lock;  //Main mutex for threads
	pthread_rwlock_t exportedClauseLock;
	int nthreads;


}; //class SolverGroup

} //namespace Minisat

#endif //Minisat_SolverGroup_h
