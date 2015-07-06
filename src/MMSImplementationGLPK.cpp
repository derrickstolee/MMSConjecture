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
 * MMSImplementationGLPK.cpp
 *
 *  Created on: Jul 24, 2012
 *      Author: stolee
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stack>
#include "Set.hpp"
#include "translation.hpp"
#include "MMSBranchingManager.hpp"

#include "glpk.h"

void MMSBranchingManager::initLP()
{
	/*  initialize the LP */
	glp_term_hook(NULL, NULL);
	glp_term_out(0);

	this->lp = glp_create_prob();
	glp_add_cols(this->lp, this->n);

	this->parm = (glp_smcp*) malloc(sizeof(glp_smcp));
	glp_init_smcp(this->parm);
	this->parm->msg_lev = GLP_MSG_OFF; //GLP_MSG_ERR
	this->parm->meth = GLP_DUALP;
	this->parm->presolve = GLP_ON;

	int* row_indices = (int*) malloc((this->n + 1) * sizeof(int));
	double* row_values = (double*) malloc((this->n + 1) * sizeof(double));
	for ( int i = 1; i <= this->n; i++ )
	{
		// every variable is unbounded
		glp_set_col_bnds(this->lp, i, GLP_FR, -100.0, 100.0);

		row_indices[i] = i;
		row_values[i] = 1.0;
	}

	/* make the sum have value 0 */
	glp_add_rows(this->lp, 1);
	glp_set_row_name(lp, 1, "sum");
	glp_set_row_bnds(lp, 1, GLP_LO, 0.0, 0.0);
	glp_set_mat_row(this->lp, 1, this->n, row_indices, row_values);

	glp_add_rows(this->lp, this->n - 1);
	for ( int i = 1; i < this->n; i++ )
	{
		char buffer[500];
		int row_is[3];
		double row_vals[3];

		/* x_i - x_{i+1} >= 0 */
		row_is[0] = 0;
		row_is[1] = i;
		row_is[2] = i + 1;

		row_vals[0] = 0.0;
		row_vals[1] = 1.0;
		row_vals[2] = -1.0;

		sprintf(buffer, "ineq%d", i);

		glp_set_row_name(lp, 1 + i, buffer);
		glp_set_mat_row(this->lp, 1 + i, 2, row_is, row_vals);
		glp_set_row_bnds(lp, 1 + i, GLP_LO, 0.0, 0.0);
	}

	this->num_rows = this->n;

	free(row_indices);
	free(row_values);

	this->stack_lp_rows.push(this->num_rows);

	// objective function
	glp_set_obj_dir(this->lp, GLP_MIN);
	glp_set_obj_coef(this->lp, 1, 1.0);

	this->solution = (double*) malloc(this->n * sizeof(double));

}

void MMSBranchingManager::freeLP()
{
	/* free the LP */
	if ( this->lp != 0 )
	{
		glp_delete_prob(this->lp);
		free(this->solution);

		free(this->parm);
	}
}

void MMSBranchingManager::clearLP()
{
	/* clear the LP */
	if ( this->stack_lp_rows.size() > 0 )
	{
		this->stack_lp_rows.pop();
	}
	this->stack_lp_rows.push(this->num_rows);

	int* row_indices = (int*) malloc((this->num_rows + 1) * sizeof(int));

	for ( int i = 1; i <= this->num_rows - this->n; i++ )
	{
		row_indices[i] = i + this->n;
	}

	if ( this->num_rows - this->n > 0 )
	{
		glp_del_rows(this->lp, this->num_rows - this->n, row_indices);
	}

	this->num_rows = this->n;

	free(row_indices);


//	this->freeLP();
//	this->initLP();
}

void MMSBranchingManager::addPositiveConstraint( int rank )
{
	/* Add constraints to the LP */
	int* set = (int*) malloc(this->k * sizeof(int));
	int* row_indices = (int*) malloc((this->k + 1) * sizeof(int));
	double* row_values = (double*) malloc((this->k + 1) * sizeof(double));

	lexIndexToSet(this->n, this->k, rank, set);

	row_indices[0] = 0;
	row_values[0] = 0.0;

	/* need to encode sum_{i \in S} x_i \geq 0 */
	for ( int i = 0; i < this->k; i++ )
	{
		row_indices[i + 1] = set[i] + 1;
		row_values[i + 1] = 1.0;
	}

	this->num_rows = this->num_rows + 1;
	glp_add_rows(this->lp, 1);
	glp_set_mat_row(this->lp, this->num_rows, this->k, row_indices, row_values);
	glp_set_row_bnds(this->lp, this->num_rows, GLP_LO, 0.0, 0.0);

	free(set);
	free(row_indices);
	free(row_values);
}

void MMSBranchingManager::addNegativeConstraint( int rank )
{
	/* Add constraints to the LP */
	int* set = (int*) malloc(this->k * sizeof(int));
	int* row_indices = (int*) malloc((this->k + 1) * sizeof(int));
	double* row_values = (double*) malloc((this->k + 1) * sizeof(double));

	lexIndexToSet(this->n, this->k, rank, set);

	row_indices[0] = 0;
	row_values[0] = 0.0;

//	printf("Adding negative constraint on the set ");
	/* need to encode sum_{i \in S} x_i \leq -1 */
	for ( int i = 0; i < this->k; i++ )
	{
//		printf("%2d ", set[i] + 1);
		row_indices[i + 1] = set[i] + 1;
		row_values[i + 1] = 1.0;
	}
//	printf("\n");
	free(set);

	this->num_rows = this->num_rows + 1;
	glp_add_rows(this->lp, 1);
	glp_set_mat_row(this->lp, this->num_rows, this->k, row_indices, row_values);
	glp_set_row_bnds(this->lp, this->num_rows, GLP_UP, -1.0, -1.0);

	free(row_indices);
	free(row_values);
}

bool MMSBranchingManager::isLPFeasible()
{
	glp_adv_basis(this->lp, 0);
	int result = glp_exact(this->lp, this->parm);
//
//	if ( result == GLP_EBADB )
//	{
//		// basis problem... what is going on?
//		printf("-- The solver had a basis problem!\n");
//		return false;
//	}
//
//	if ( result == GLP_ESING )
//	{
//		printf("--The basis was singular...\n");
//		return false;
//	}
//
//	if ( result == GLP_ENOPFS )
//	{
//		printf("-- Presolver found infeasible!\n");
//		return false;
//	}

	if ( result != 0 )
	{
		printf("--The Solver had an error! %d\n", result);
//		return false;
	}

	result = glp_get_status(this->lp);
	if ( result == GLP_INFEAS || result == GLP_NOFEAS )
	{
//		printf("-- infeasible!\n");
		return false;
	}

	return true;
}

bool MMSBranchingManager::isLPOptimalACounterexample()
{
//	printf("\t\t\tchecking for counterexample\n");
	//  check for counterexample!
//	double lowest_abs_diff = 10.0;
	for ( int i = 0; i < this->n; i++ )
	{
		this->solution[i] = glp_get_col_prim(this->lp, i + 1);

//		for ( int j = 0; j < i; j++ )
//		{
//
//			if ( abs(abs(this->solution[i]) - abs(this->solution[j])) < lowest_abs_diff )
//			{
//				lowest_abs_diff = abs(abs(this->solution[i]) - abs(this->solution[j]));
//			}
//		}
	}

	/* check all sets */
	int* set = (int*) malloc(this->k * sizeof(int));

	int num_definitely_pos = 0;
	int num_pos = 0;

	lexIndexToSet(this->n, this->k, 0, set);
	for ( int rank = 0; rank < this->nchoosek; rank++ )
	{
		double sum = 0.0;

		for ( int j = 0; j < this->k; j++ )
		{
			sum += solution[set[j]];
		}

		if ( sum >= 0.0 )
		{
			num_definitely_pos++;
		}
		if ( sum >= -0.000001 )
		{
			num_pos++;
		}

		getLexSuccessor(this->n, this->k, set);
	}

	free(set);

	this->sol_pos_sets = num_pos;
//	printf("--found %10d non-negative sets and want %10d.\n", num_pos, this->mms_target);

	if ( num_definitely_pos >= this->mms_target )
	{
		return false;
	}

	if ( num_pos >= this->mms_target )
	{
		printf("-- THE FOLLOWING OPTIMAL VECTOR HAS %d NON-NEGATIVE SETS, and %d SETS >= 0.00001 (Target: %d):\n",
		       num_definitely_pos, num_pos, this->mms_target);
		this->printLPOptimal();
		return false;
	}

	return true;
}

void MMSBranchingManager::printLPOptimal()
{
	for ( int i = 0; i < this->n; i++ )
	{
		this->solution[i] = glp_get_col_prim(this->lp, i + 1);
	}

	printf("optimal: \n");
	for ( int i = 0; i < this->n; i++ )
	{
		printf("\tx_{%2d} = %3.6lf\n", i, this->solution[i]);
	}
	printf(" with %d non-negative sets (target: %d)\n\n", this->sol_pos_sets, this->mms_target);
}

void MMSBranchingManager::snapshotLP()
{
	/* make a snapshot of the current LP */
	int num_before = this->stack_lp_rows.top();

	this->stack_lp_rows.push(this->num_rows);
}

void MMSBranchingManager::rollbackLP()
{
	/* Rollback the LP to the previous snapshot */
	int num_after = this->stack_lp_rows.top();
	this->stack_lp_rows.pop();
	int num_to_delete = this->num_rows - num_after;

	if ( num_to_delete > 0 )
	{
		int* row_indices = (int*) malloc((num_to_delete + 1) * sizeof(int));

		row_indices[0] = -1;
		for ( int i = 1; i <= num_to_delete; i++ )
		{
			row_indices[i] = num_after + i;
		}

		glp_del_rows(this->lp, num_to_delete, row_indices);

		free(row_indices);
	}

	this->num_rows = num_after;
}

