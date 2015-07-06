/***********************************************************

 Copyright Derrick Stolee 2012.

 This file is part of SearchLib.

 SearchLib is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 SearchLib is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with SearchLib.  If not, see <http://www.gnu.org/licenses/>.

 *************************************************************/

/*
 * MMSImplementationCPLEX.cpp
 *
 * This is the implementation of the LP objects for CPLEX mode.
 *
 *  Created on: Jul 24, 2012
 *      Author: stolee
 */

#define CPLEX

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stack>
#include <math.h>
#include "Set.hpp"
#include "translation.hpp"
#include "MMSBranchingManager.hpp"

extern "C"
{
#include "ilcplex/cplex.h"
}
typedef unsigned long long int ULI;

void MMSBranchingManager::initLP()
{
	/* initialize the LP */
	int status = 0;
	this->env = CPXopenCPLEX(&status);
	lpx = CPXcreateprob(this->env, &status, "MMS Conjecture: Branch-and-Cut Method");
	CPXchgobjsen(this->env, this->lpx, CPX_MIN);

	/**
	 * We have n unbounded real variables: x_1, ..., x_n.
	 */
	double *obj = new double[n]; // objective values
	double *lb = new double[n]; // lower bounds
	double *ub = new double[n]; // upper bounds

	for ( int i = 0; i < n; i++ ) // CPLEX is 0-indexed
	{
		obj[i] = 0.0;
		lb[i] = -10000.0;
		ub[i] = 10000.0;
	}
	obj[0] = 1.0; /* minimize x_1 */
	CPXnewcols(this->env, this->lpx, this->n, obj, lb, ub, NULL, NULL);
	// have xctype passed as NULL; everything continuous, so not MIP
	delete obj;
	obj = 0;
	delete lb;
	lb = 0;
	delete ub;
	ub = 0;

	this->solution = new double[this->n];
	this->isolution = new int[this->n];
	bzero(this->solution, this->n * sizeof(double));
	bzero(this->solution, this->n * sizeof(int));

	/**
	 * Now, add the constraints
	 * 0: x_1 + ... + x_n = 0
	 * 1: x_1 - x_2 >= 0
	 * 2: x_2 - x_3 >= 0
	 * ...
	 * n-1: x_{n-1} - x_n >= 0
	 */
	char* sense = new char[this->n]; // senses
	int *rmatbeg = new int[this->n]; // where do rows begin in sparse representation?
	int *rmatind = new int[3 * this->n - 2]; // what columns are used?
	double *rmatval = new double[3 * this->n - 2]; //  what values are used?

	sense[0] = 'E';
	rmatbeg[0] = 0;
	for ( int i = 0; i < this->n; i++ )
	{
		rmatind[i] = i;
		rmatval[i] = 1.0;
	}
	rmatbeg[1] = this->n;
	for ( int i = 1; i < this->n; i++ )
	{
		sense[i] = 'G';

		if ( i > 1 )
		{
			rmatbeg[i] = rmatbeg[i - 1] + 2;
		}

		rmatind[rmatbeg[i] + 0] = i - 1;
		rmatind[rmatbeg[i] + 1] = i;

		rmatval[rmatbeg[i] + 0] = 1.0;
		rmatval[rmatbeg[i] + 1] = -1.0;
	}
	status = CPXaddrows(this->env, this->lpx, 0, this->n, 3 * this->n - 2, NULL, sense, rmatbeg, rmatind, rmatval, NULL,
	                    NULL);
	delete sense;
	delete rmatbeg;
	delete rmatind;
	delete rmatval;

	this->num_constraints = this->n;
	this->stackLPconstraints.push(this->num_constraints);
	this->stackLPconstraints_snapshots.push(this->stackLPconstraints.size());

	//status=CPXsetintparam(env,CPX_PARAM_SCRIND,CPX_ON);
	//status=CPXsetintparam(env,CPX_PARAM_SIMDISPLAY,2);
	status = CPXlpopt(this->env, this->lpx);
	double objval = 0;
	status = CPXgetobjval(this->env, this->lpx, &objval);

	status = CPXlpopt(this->env, this->lpx);
}
void MMSBranchingManager::freeLP()
{
	/* free the LP */
	int status = CPXdelrows(this->env, this->lpx, 0, this->num_constraints);
//	delete this->lpx;
//	delete this->env;
	this->env = 0;
	this->lpx = 0;

	this->num_constraints = 0;

	delete this->solution;
	delete this->isolution;

	while ( this->stackLPconstraints.size() > 0 )
	{
		this->stackLPconstraints.pop();
	}
	while ( this->stackLPconstraints_snapshots.size() > 0 )
	{
		this->stackLPconstraints_snapshots.pop();
	}
}

void MMSBranchingManager::clearLP()
{
	/* clear the LP in a silly way... */
	this->freeLP();
	this->initLP();
}

void MMSBranchingManager::addPositiveConstraint( int rank )
{
	/* Add constraints to the LP */
	/**
	 * Now, add the constraint
	 * 0: sum( x_i : i in A_rank ) >= 0
	 */
	char* sense = new char[1]; // senses
	int *rmatbeg = new int[1]; // where do rows begin in sparse representation?
	int *rmatind = new int[this->k]; // what columns are used?
	double *rmatval = new double[this->k]; //  what values are used?

	sense[0] = 'G';
	rmatbeg[0] = 0;
	int* A = new int[this->k];

	lexIndexToSet(this->n, this->k, rank, A);

//	printf(" adding positive constraint: { ");
	for ( int i = 0; i < this->k; i++ )
	{
//		printf("%2d ", A[i]);
		rmatind[i] = A[i];
		rmatval[i] = 1.0;
	}
//	printf("}\n");

	int status = CPXaddrows(this->env, this->lpx, 0, 1, this->k, NULL, sense, rmatbeg, rmatind, rmatval, NULL, NULL);
	delete sense;
	delete rmatbeg;
	delete rmatind;
	delete rmatval;
	delete A;

	this->stackLPconstraints.push(this->num_constraints);
	this->num_constraints = this->num_constraints + 1;
}

void MMSBranchingManager::addNegativeConstraint( int rank )
{
	/* Add constraints to the LP */
	/**
	 * Now, add the constraint
	 * 0: sum( x_i : i in A_rank ) <= -1
	 */
	char* sense = new char[1]; // senses
	double* rhs = new double[1]; // right hand side
	int *rmatbeg = new int[1]; // where do rows begin in sparse representation?
	int *rmatind = new int[this->k]; // what columns are used?
	double *rmatval = new double[this->k]; //  what values are used?

	sense[0] = 'L';
	rmatbeg[0] = 0;
	rhs[0] = -1;
	int* A = new int[this->k];

	lexIndexToSet(this->n, this->k, rank, A);

//	printf(" adding negative constraint: { ");
	for ( int i = 0; i < this->k; i++ )
	{
//		printf("%2d ", A[i]);
		rmatind[i] = A[i];
		rmatval[i] = 1.0;
	}
//	printf("}\n");

	int status = CPXaddrows(this->env, this->lpx, 0, 1, this->k, rhs, sense, rmatbeg, rmatind, rmatval, NULL, NULL);
	delete sense;
	delete rmatbeg;
	delete rmatind;
	delete rmatval;
	delete rhs;
	delete A;

	this->stackLPconstraints.push(this->num_constraints);
	this->num_constraints = this->num_constraints + 1;
}

bool MMSBranchingManager::isLPFeasible()
{
	//status=CPXsetintparam(env,CPX_PARAM_SCRIND,CPX_ON);
	//status=CPXsetintparam(env,CPX_PARAM_SIMDISPLAY,2);
	int status = CPXlpopt(this->env, this->lpx);
	int lpstat = CPXgetstat(this->env, this->lpx);

//	printf("RAN CPLEX WITH STATUS %d and LPSTAT %d...", status, lpstat);

	if ( lpstat == CPX_STAT_OPTIMAL )
	{
//		printf("OPTIMAL.\n");
		/* time to store the solution */
		status = CPXsolution(this->env, this->lpx, &lpstat, NULL, this->solution, NULL, NULL, NULL);
		return true;
	}
	else
	{
		//printf("OTHER: status %d lpstat %d\n", lpstat);
		/* what other statuses could happen during success? */
		return false;
	}
}

bool MMSBranchingManager::isLPOptimalACounterexample()
{
	/* check the solution */
	int num_pos = 0;
	int num_i_pos = 0;
	int num_weak_pos = 0;
	int* set = new int[this->k];
	double min_weak_value = -1.0;

	for ( int i = 0; i < this->n; i++ )
	{
		this->isolution[i] = round(this->solution[i] * 27720);
	}

	lexIndexToSet(this->n, this->k, 0, set);
	for ( int i = 0; i < this->nchoosek; i++ )
	{

		double value = 0.0;
		int ival = 0;

		for ( int j = 0; j < this->k; j++ )
		{
			value += this->solution[set[j]];
			ival += this->isolution[set[j]];
		}

		if ( value >= -0.0 )
		{
			num_pos++;
		}

		if ( ival >= 0 )
		{
			num_i_pos++;
		}

		if ( value >= -0.0001 )
		{
			num_weak_pos++;

			if ( value < 0.0 && value > min_weak_value )
			{
				min_weak_value = value;
			}
		}

		getLexSuccessor(this->n, this->k, set);
	}

	delete set;
	this->num_positive = num_weak_pos;
	this->num_sure_positive = num_pos;
	this->num_integer_positive = num_i_pos;

	for ( int i = 0; i < this->n; i++ )
	{
		this->isolution[i] = round(this->solution[i] * 1000);
	}

	if ( num_weak_pos < this->mms_target )
	{
		if ( num_pos < this->mms_target )
		{
			printf("-- FOUND A STRONG(?) EXAMPLE with MIN WEAK VALUE %2.15lf:\n", min_weak_value);
			return true;
		}

		printf("-- FOUND A WEAK EXAMPLE with MIN WEAK VALUE %2.15lf:\n", min_weak_value);
		this->printLPOptimal();
		printf("\t\t(%5d of those sums were surely non-negative!)\n", num_pos);
		printf("\t\t(%5d of those sums were non-negative as integers!)\n", num_pos);
		return false;
	}

	return false;
}

void MMSBranchingManager::printLPOptimal()
{
	double sum_val = 0.0;

	int isum = 0;

	for ( int i = 0; i < this->n; i++ )
	{
		isum += this->isolution[i];
	}

	if ( isum < 0 )
	{
		printf("-- This example has negative integer sum: %d\n", isum);
	}

	for ( int i = 0; i < this->n; i++ )
	{
		sum_val += this->solution[i];
		printf("\tx_{%2d} = %3.20lf\t1000x_{%2d} = %5d\n", i, this->solution[i], i, this->isolution[i]);
	}
	printf("\t...has sum %3.20lf \n\t..and %5d non-negative %d-sums (%d integer non-neg), but we wanted %5d.\n",
	       sum_val, this->num_positive, this->k, this->num_integer_positive, this->mms_target);

}

void MMSBranchingManager::snapshotLP()
{
//	for ( int i = 0; i < this->stackLPconstraints_snapshots.size(); i++ )
//	{
//		printf(" ");
//	}
//	printf("snapshot with %d constraints\n", this->num_constraints);

	/* make a snapshot of the current LP */
	this->stackLPconstraints_snapshots.push(this->stackLPconstraints.size());
}

void MMSBranchingManager::rollbackLP()
{
	/* Rollback the LP to the previous snapshot */
	int to_size = this->stackLPconstraints_snapshots.top();
	this->stackLPconstraints_snapshots.pop();

	int to_num_rows = this->num_constraints;
	while ( this->stackLPconstraints.size() > to_size )
	{
		to_num_rows = this->stackLPconstraints.top();
		this->stackLPconstraints.pop();
	}

	/* delete rows from to_num_rows to this->num_constraints */
//	printf("-- deleting rows %d through %d...", to_num_rows, this->num_constraints);
	int status = CPXdelrows(this->env, this->lpx, to_num_rows, this->num_constraints - 1);

//	if ( status == 0 )
//	{
//		printf("no error.\n");
//	}
//	else
//	{
//		printf("with an error %d!\n", status);
//	}

	this->num_constraints = to_num_rows;
//	for ( int i = 0; i < this->stackLPconstraints_snapshots.size(); i++ )
//	{
//		printf(" ");
//	}
//	printf("rollback with %d constraints\n", this->num_constraints);

}

