#ifndef Minisat_SolverGroup_H
#define Minisat_SolverGroup_H

#include <vector>
#include "SimpSolver.h"
#include <pthread.h>
#include <vector>
#include <algorithm>
#include <math.h>

namespace Minisat {

class SimpSolver;

class SolverGroup;

struct sg_thread_status {
	SolverGroup     *group;
	vec<Lit>	*assumps;
	bool		done;
	lbool		result;
	pthread_cond_t	*signal_complete;
	pthread_mutex_t *lock;
	int 		thread_id;
	SimpSolver 	*solver;
	bool		exit_now;
	int 		mode;
	int		guiding_path_num;
};

enum stat {
	INDETERM,
	SAT,
	UNSAT
};

enum solver_mode {
	MODE_GP,
	MODE_RAND
};

struct VarCount {
	int var;
	int num_pos;
	int num_neg;
};

struct sharedClause {
	int thread_id;
	vec<Lit> *clause;
};

class SolverGroup {
public:
	SolverGroup(int nthreads);
	SolverGroup();
	~SolverGroup();

	void parse_DIMACS(gzFile in);
	lbool solve_parallel(int *winning_thread);
	bool eliminate_parallel(bool in);
	vec< vec<Lit>* > *createGuidingPaths(int num);
	bool processCompleteSolvers( lbool &status, int *winning_thread );
	int  findOldThread();
	void exportSharedClause( int thread_id, vec<Lit> &clause );
	void printVarStats(); 
	bool getNextSharedClause(int thread_id, vec<Lit> &clause);
	bool getNextSharedUnit(int thread_id, Lit &out);
	void exportSharedUnit( Lit &in );
	void resetIters( int thread_id);

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
	pthread_rwlock_t sharedClauseLock;
	pthread_rwlock_t sharedUnitLock;
	int nthreads;
	int current_guiding_path;
	int *thread_guiding_path; //Array containing relation of which thread is
				  //working on which guiding path
	int *guiding_path_status;
	int num_guiding_paths;
	vec< vec<Lit>* > *guiding_path_queue;
	std::vector< VarCount > var_stats;
	std::vector< struct sharedClause > shared_clauses;
	std::vector< Lit > shared_units;
	unsigned int *shared_clause_iters;
	unsigned int *shared_unit_iters;


}; //class SolverGroup

} //namespace Minisat

#endif //Minisat_SolverGroup_h
