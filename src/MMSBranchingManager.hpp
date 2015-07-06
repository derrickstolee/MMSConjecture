/*
 * MMSBranchingManager.hpp
 *
 * The MMSBranchingManager class is a TreeSearch implementation of a branching
 * procedure to verify the MMS Conjecture for a pair (n,k) using linear programming.
 *
 *  Created on: Jul 19, 2012
 *      Author: stolee
 */

//#define CPLEX
// #define GLPK
// #define COINOR
// #define GUROBI
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stack>
#include "Set.hpp"
#include "SearchManager.hpp"

#ifdef CPLEX
extern "C"
{
#include "ilcplex/cplex.h"
}
#endif

#ifdef GLPK
extern "C"
{
#include <glpk.h>
}
#endif

#ifdef COINOR

#endif

#ifdef GUROBI

#endif

class MMSNode: public SearchNode
{
	public:
		MMSNode( LONG_T label );
		virtual ~MMSNode();
		int branch_rank;
};

class MMSBranchingManager: public SearchManager
{
	protected:
		typedef SearchManager super;

		bool write_proof;

		/**
		 * n and k are the parameters of the conjecture we are testing.
		 * c and r are such that n = ck + r (r in 0...k-1)
		 */
		int n, k, c, r;
		int nchoosek;

		/**
		 * gvals[m] stores g(m,k).
		 *
		 * Initialized using Baranyai.
		 */
		int* gvals;

		/**
		 * Store {n-1 choose k-1}, or whatever target we are looking for.
		 */
		int mms_target;

		/**
		 * Store { (c-1)k - 1 choose k-1 } where n = ck+r.
		 */
		int Baranyai_value;

		bool strong_mode;

		/* used for pruning AFTER a solution was found, avoids duplicates, hopefully */
		bool found_solution;
		int found_solution_at;

		char* markers;
		int initMarkers();
		void freeMarkers();
		void clearMarkers();

		/**
		 * The left- and right- neighbors define the Hasse diagram of
		 * the left-shift poset.
		 * Given as _ranks_, not sets.
		 * Index is always k*rank, even if not all previous ranks had k neighbors.
		 */
		int* left_neighbors;
		int* right_neighbors;
		int* left_degree;
		int* right_degree;
		bool* contains0; /* does the set with this rank contain 0? */
		int initNeighbors();
		void freeNeighbors();

		int getMeetRank( int num_sets, int* set_ranks );
		int getJoinRank( int num_sets, int* set_ranks );

		/**
		 * The left- and right-shift values stored in a table.
		 */
		int* left_shift;
		int* right_shift;
		int initLeftShift();
		int initRightShift();
		void freeLeftShift();
		void freeRightShift();
		/**
		 * The left shift function calculates how many subsets of the given size are
		 * to the left of the given set (including that set). It is a recursive formula, using
		 * the formula on lower sizes.
		 */
		int getLeftShift( int* set, int size, int rank );

		/**
		 * The right shift function calculates how many subsets of the given size are to the right of
		 * the given set (including that set). This is only used for sets of size k and simply reverses
		 * the set {1,...,n} and asks for the left shift.
		 */
		int getRightShift( int* set, int rank );

		/**
		 * These tables will store the number of certain values to the right or left of a given set THAT IS IN C^*
		 *
		 * When a set changes from C^* to C^+ or C^-, the values are recomputed.
		 * This is part of the addPositiveSet and addNegativeSet methods
		 */
		bool use_breadth_first_search_right;
		bool use_inclusion_exclusion_right;
		bool use_incremental_method_right;
		bool use_breadth_first_search_left;
		bool use_inclusion_exclusion_left;
		bool use_incremental_method_left;

		int max_depth_left_c_star;
		int cur_depth_left_c_star;
		int** left_c_star;
		int* left_c_star_updated_to; // stores the number of the POSITIVE sets the current level is updated to.

		bool use_right_c_star;
		int max_depth_right_c_star;
		int cur_depth_right_c_star;
		int** right_c_star;
		int* right_c_star_updated_to; // stores the number of the NEGATIVE sets the current level is updated to.
		int initCStarValues();
		void freeCStarValues();
		void clearCStarValues();
		bool recalculating_c_star_values;
		void snapshotCStarValues();
		void rollbackCStarValues();
		std::stack<int> stack_left_c_star_depth;
		std::stack<int> stack_right_c_star_depth;

		void recalculateAllCStarValues( int mode );

		void recalculateRightCStarBFS();
		void recalculateRightCStarInclusionExclusion();
		void recalculateRightCStarIterative( int* new_set );

		void recalculateLeftCStarBFS();
		void recalculateLeftCStarInclusionExclusion();
		void recalculateLeftCStarIterative( int* new_set );

		int getLeftCStarValue( int rank );
		int getRightCStarValue( int rank );

		/* used in BFS method */
		int recalculateLeftCStarValueRecurse( int rank );

		/* used in BFS Method */
		int recalculateRightCStarValueRecurse( int rank );

		int num_forced_positive;
		void addNewPositiveSets( int number );
		void snapshotNumPos();
		void rollbackNumPos();
		void clearNumPos();
		std::stack<int> stackNumPos;

		/***
		 * Store the list of sets that generate C^+ and C^-
		 */
		int num_genpos;
		int max_genpos;
		int* genpos;
		int num_genneg;
		int max_genneg;
		int* genneg;
		void initGenSets();
		void freeGenSets();
		void clearGenSets();
		void snapshotGenSets();
		void rollbackGenSets();
		std::stack<int> stack_num_genpos;
		std::stack<int> stack_num_genneg;

		/**
		 * The C^+, C^- and C^* sets contain the current sets of positively-selected
		 * k-sets, negatively-selected k-sets, and currently-open k-sets.
		 *
		 * These take different values within the set_labels array.
		 *
		 * 0 - C^*
		 * 1 - C^+
		 * 2 - C^-
		 */
#define LABEL_C_STAR 0
#define LABEL_C_POS 1
#define LABEL_C_NEG 2
		char* set_labels;
		int min_c_star_rank;
		int max_c_star_rank;
		int num_in_c_star;
		int initSetLabels();
		void freeSetLabels();
		void snapshotSetLabels();
		void rollbackSetLabels();
		void clearSetLabels();
		std::stack<int> stack_max_c_star;
		std::stack<int> stack_min_c_star;
		std::stack<int> stack_set_labels;
		std::stack<int> stack_set_label_sizes;
		std::stack<int> stack_num_in_c_star;

#ifdef CPLEX
		/**
		 * The LP instance in terms of CPLEX
		 */
		int num_constraints;
		double* solution;
		int* isolution;
		int num_positive;
		int num_sure_positive;
		int num_integer_positive;
		CPXENVptr env;
		CPXLPptr lpx;
#endif

		int sol_pos_sets;

#ifdef GLPK
		/**
		 * The LP instance in terms of GLPK
		 */
		// TODO: add GLPK THINGS
		glp_prob *lp;
		glp_smcp* parm;
		int num_rows;
		double* solution;
		std::stack<int> stack_lp_rows;
#endif

#ifdef COINOR
		/**
		 * The LP instance in terms of COINOR
		 */
		// TODO: add COINOR THINGS
#endif

#ifdef GUROBI
		/**
		 * The LP instance in terms of GUROBI
		 */
		// TODO: add GUROBI THINGS
#endif

		void initLP();
		void freeLP();
		void clearLP();
		void addPositiveConstraint( int rank );
		void addNegativeConstraint( int rank );
		bool isLPFeasible();
		bool isLPOptimalACounterexample();
		void printLPOptimal();
		void snapshotLP();
		void rollbackLP();
		// TODO: does this stack type work with the LP implementation?
		std::stack<int> stackLPconstraints;
		std::stack<int> stackLPconstraints_snapshots;

		/**
		 * Handle all of the memory things.
		 */
		void initAll();
		void freeAll();
		void snapshot();
		void rollback();

		/**
		 * Select a set to branch upon.
		 */
#define BRANCH_MIDDLE 1
#define BRANCH_BALANCED 2
#define BRANCH_OPTHALF 4
#define BRANCH_LEFTHALF 8
#define BRANCH_LMAXPOS 16
#define BRANCH_LMINNEG 32
		char branching_rule;
		int broke_optimal;
		int getBranchSetMiddle();
		int getBranchSetBalanced();
		int getBranchSetLeftHalf();
		double left_divisor;

		// LP-based methods:

		int getBranchSetOptHalf();
		int getBranchSetMaxLPositive();
		int getBranchSetMinLNegative();

		/**
		 * given the current situation, determine if some sets _must_ be negative in a counterexample.
		 */
#define PROPAGATE_NONE 0
#define PROPAGATE_GAC 1
#define PROPAGATE_BRANCH_SGAC 2
#define PROPAGATE_SGAC 4
		char propagation_mode;
		bool propagate_positive;
		bool propagateGAC();
		bool propagateBranchSGAC();
		bool propagateSGAC();
		bool propagatePositive();

		bool lp_disabled;

		bool use_branchless;
		bool use_stochastic;
		bool done_branchless_or_stochastic;

	public:
		/**
		 * Default constructor
		 */
		MMSBranchingManager( int n, int k, bool write_proof );

		/**
		 * Destructor
		 */
		virtual ~MMSBranchingManager();

		/**
		 * pushNext -- deepen the search to the next child
		 * 	of the current node.
		 *
		 * @return the label for the new node. -1 if none.
		 */
		virtual LONG_T pushNext();

		/**
		 * pushTo -- deepen the search to the specified child
		 * 	of the current node.
		 *
		 * @param child the specified label for the new node
		 * @return the label for the new node. -1 if none, or failed.
		 */
		virtual LONG_T pushTo( LONG_T child );

		/**
		 * pop -- remove the current node and move up the tree.
		 *
		 * @return the label of the node after the pop.
		 * 		This return value is used for validation purposes
		 * 		to check proper implementation of push*() and pop().
		 */
		virtual LONG_T pop();

		/**
		 * prune -- Perform a check to see if this node should be pruned.
		 *
		 * @return 0 if no prune should occur, 1 if prune should occur.
		 */
		virtual int prune();

		/**
		 * isSolution -- Perform a check to see if a solution exists
		 * 		at this point.
		 *
		 *  By "Solution" we of course mean a "counterexample" to the conjecture.
		 *  Test the LP feasible point for a counterexample.
		 *
		 * @return 0 if no solution is found, 1 if a solution is found.
		 */
		int isSolution();

		/**
		 * writeSolution -- create a buffer that contains a
		 * 	description of the solution.
		 *
		 * Write the counterexample, as well as the number of non-negative k-sets and how that violates expectation
		 *
		 * @return a string of data allocated with malloc().
		 * 	It will be deleted using free().
		 */
		virtual char* writeSolution();

		/**
		 * writeStatistics -- create a buffer that contains a
		 * 	description of the solution.
		 *
		 * Statistics take the following format in each line:
		 *
		 * T [TYPE] [ID] [VALUE]
		 *
		 * @return a string of data allocated with malloc().
		 * 	It will be deleted using free().
		 */
		virtual char* writeStatistics();

		/**
		 * ClearAll will restart the search from scratch, but won't reallocate.
		 * It just clears the data to zero.
		 */
		void clearAll();

		/**
		 * setBranchingRule
		 */
		void setBranchingRule( char rule );
		void setLeftDivisor( double divisor ); // used for lefthalf rule

		/**
		 * setPropagationMode
		 */
		void setPropagationMode( char mode );

		/**
		 * setPropagationPositive
		 */
		void setPropagationPositive( bool proppos );

		/**
		 * setProofMode
		 *
		 * Should we write the proof to output?
		 */
		void setProofMode( bool mode );

		/**
		 * setTarget
		 *
		 * Change the target number of non-negative k-sums to avoid.
		 *
		 * Useful for computing g(n,k) when it is strictly less than (n-1 choose k-1).
		 */
		void setTarget( int target );

		/**
		 * setStrongMode
		 *
		 * Change whether or not to test for strong pair (n,k).
		 */
		void setStrongMode( bool mode );

		/**
		 * setGValue
		 *
		 * Store a previously-computed value of g(m,k) for some other m < n.
		 */
		void setGValue( int m, int value );

		/**
		 * useRightCStar()
		 *
		 * Should we use right c-star values or not?
		 */
		void useRightCStar( bool useRCS );

		/**
		 * useIncrementalMethod()
		 *
		 * Set us up for the Incremental Method (and appropriately changes the other settings).
		 */
		void useIncrementalMethodLeft();
		void useIncrementalMethodRight();

		/**
		 * useInclusionExclusion()
		 *
		 * Set us up for the Inclusion/Exclusion Method (and appropriately changes the other settings).
		 */
		void useInclusionExclusionLeft();
		void useInclusionExclusionRight();

		/**
		 * useBFS()
		 *
		 * Set us up for the BFS Method (and appropriately changes the other settings).
		 */
		void useBFSLeft();
		void useBFSRight();

		/**
		 * disableLP ONLY FOR TESTING PURPOSES!
		 */
		void disableLP();

		/**
		 * Propagate.
		 *
		 * Return FALSE if infeasible!
		 */
		bool propagate();

		/* public for tests */
		void branchPositiveSet( int rank );
		void branchNegativeSet( int rank );
		/**
		 * 0 for both, < 0 for only left, > 0 for only right
		 */
		int addPositiveSet( int rank );
		int addNegativeSet( int rank );
		int addMinimalPositiveSet( int rank );
		int addMaximalNegativeSet( int rank );

		int getBranchSet();

		/** test for the correct left-shift values */
		bool testAgainst( MMSBranchingManager* obj1, MMSBranchingManager* obj2 );

		/**
		 * doBranchlessSearch
		 *
		 * Do ALL negative propagations, then do ALL positive propagations.
		 * REPEAT.
		 *
 	     * returns FALSE if it reaches a stable point with propagations to make.
		 */
		bool doBranchlessSearch();
		void useBranchlessSearch(bool branchless);

		/**
		 * doStochasticSearch
		 *
		 * Do ALL negative propagations, then randomly select sets to check if they make the LP infeasible
		 * when negative. Make these sets positive, propagate negatives, and repeat!
		 *
		 */
		bool doStochasticSearch();
		void useStochasticSearch(bool stochastic);
};

