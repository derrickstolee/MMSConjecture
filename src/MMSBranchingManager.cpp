/*
 * MMSBranchingManager.cpp
 *
 *  Created on: Jul 23, 2012
 *      Author: stolee
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <stack>
#include <queue>
#include "Set.hpp"
#include "translation.hpp"
#include "shiftcalculations.hpp"
#include "MMSBranchingManager.hpp"

MMSNode::MMSNode( LONG_T label ) :
		SearchNode(label)
{
	this->branch_rank = -1;
}

MMSNode::~MMSNode()
{

}

int MMSBranchingManager::initMarkers()
{
	this->markers = (char*) malloc(this->nchoosek * sizeof(int));
	bzero(this->markers, this->nchoosek * sizeof(int));

	return this->nchoosek * sizeof(int);
}

void MMSBranchingManager::freeMarkers()
{
	if ( this->markers != 0 )
	{
		free(this->markers);
		this->markers = 0;
	}
}

void MMSBranchingManager::clearMarkers()
{
	bzero(this->markers, this->nchoosek * sizeof(int));
}

/**
 * The left- and right- neighbors define the Hasse diagram of
 * the left-shift poset.
 * Given as _ranks_, not sets.
 * Index is always k*rank, even if not all previous ranks had k neighbors.
 */
int MMSBranchingManager::initNeighbors()
{
	LONG_T allocated = 0;
	this->left_neighbors = (int*) malloc(this->k * this->nchoosek * sizeof(int));
	this->right_neighbors = (int*) malloc(this->k * this->nchoosek * sizeof(int));
	this->left_degree = (int*) malloc(this->nchoosek * sizeof(int));
	this->right_degree = (int*) malloc(this->nchoosek * sizeof(int));

	allocated += this->k * this->nchoosek * sizeof(int) * 2;
	allocated += this->nchoosek * sizeof(int) * 2;

	this->contains0 = (bool*) malloc(this->nchoosek * sizeof(bool));
	for ( int i = 0; i < this->nchoosek; i++ )
	{
		this->contains0[i] = false;
	}

	allocated += this->nchoosek * sizeof(bool);

	bzero(this->left_neighbors, this->k * this->nchoosek * sizeof(int));
	bzero(this->right_neighbors, this->k * this->nchoosek * sizeof(int));
	bzero(this->left_degree, this->nchoosek * sizeof(int));
	bzero(this->right_degree, this->nchoosek * sizeof(int));

	/* place the values by building right neighbors */
	int* set = (int*) malloc(this->k * sizeof(int));
	int* rset = (int*) malloc(this->k * sizeof(int));
	int rrank = 0;

	lexIndexToSet(this->n, this->k, 0, set);
	for ( int i = 0; i < this->nchoosek; i++ )
	{
		if ( set[0] == 0 )
		{
			this->contains0[i] = true;
		}

		bcopy(set, rset, this->k * sizeof(int));

		for ( int j = 0; j < this->k; j++ )
		{
			/* can we increase this value? */
			bool is_neighbor = false;
			if ( j == k - 1 )
			{
				if ( set[j] < this->n - 1 )
				{
					is_neighbor = true;
				}
			}
			else
			{
				/* j  < k-1 */
				/* will increasing this value push us too far? */
				if ( set[j] + 1 < set[j + 1] )
				{
					is_neighbor = true;
				}
			}

			if ( is_neighbor )
			{
				rset[j] = set[j] + 1;
				rrank = lexIndexOfSet(this->n, this->k, rset);

				this->right_neighbors[this->k * i + this->right_degree[i]] = rrank;
				(this->right_degree[i]) = (this->right_degree[i]) + 1;

				this->left_neighbors[this->k * rrank + this->left_degree[rrank]] = i;
				(this->left_degree[rrank]) = (this->left_degree[rrank]) + 1;

				/* reset */
				rset[j] = set[j];
			}
		}

		getLexSuccessor(this->n, this->k, set);
	}

	free(set);
	free(rset);

	return 0;
}

void MMSBranchingManager::freeNeighbors()
{
	if ( this->left_neighbors != 0 )
	{
		free(this->left_neighbors);
		this->left_neighbors = 0;
	}
	if ( this->right_neighbors != 0 )
	{
		free(this->right_neighbors);
		this->right_neighbors = 0;
	}
	if ( this->left_degree != 0 )
	{
		free(this->left_degree);
		this->left_degree = 0;
	}
	if ( this->right_degree != 0 )
	{
		free(this->right_degree);
		this->right_degree = 0;
	}
	if ( this->contains0 != 0 )
	{
		free(this->contains0);
		this->contains0 = 0;
	}
}

int MMSBranchingManager::getMeetRank( int num_sets, int* set_ranks )
{
	int* meetset = (int*) malloc(this->k * sizeof(int));
	int* tempset = (int*) malloc(this->k * sizeof(int));
	for ( int j = 0; j < this->k; j++ )
	{
		meetset[j] = 0;
	}

	for ( int i = 0; i < num_sets; i++ )
	{
		lexIndexToSet(this->n, this->k, set_ranks[i], tempset);

		for ( int j = 0; j < this->k; j++ )
		{
			if ( tempset[j] > meetset[j] )
			{
				meetset[j] = tempset[j];
			}
		}
	}

	int meetrank = lexIndexOfSet(this->n, this->k, meetset);
	free(meetset);
	free(tempset);

	return meetrank;
}

int MMSBranchingManager::getJoinRank( int num_sets, int* set_ranks )
{
	int* joinset = (int*) malloc(this->k * sizeof(int));
	int* tempset = (int*) malloc(this->k * sizeof(int));
	for ( int j = 0; j < this->k; j++ )
	{
		joinset[j] = this->n;
	}

	for ( int i = 0; i < num_sets; i++ )
	{
		lexIndexToSet(n, this->k, set_ranks[i], tempset);
		for ( int j = 0; j < this->k; j++ )
		{
			if ( tempset[j] < joinset[j] )
			{
				joinset[j] = tempset[j];
			}
		}
	}

	int joinrank = lexIndexOfSet(this->n, this->k, joinset);
	free(joinset);
	free(tempset);

	return joinrank;
}

int MMSBranchingManager::initLeftShift()
{
	this->left_shift = (int*) malloc(this->nchoosek * sizeof(int));
	bzero(this->left_shift, this->nchoosek * sizeof(int));

	fill_leftshift_values_form(this->n, this->k, this->left_shift);

	/* the values will be calculated as required */
	return this->nchoosek * sizeof(int);
}

int MMSBranchingManager::initRightShift()
{
	this->right_shift = (int*) malloc(this->nchoosek * sizeof(int));
	bzero(this->right_shift, this->nchoosek * sizeof(int));

	int* set = (int*) malloc(this->k * sizeof(int));
	int* revset = (int*) malloc(this->k * sizeof(int));

	lexIndexToSet(this->n, this->k, 0, set);
	for ( int rank = 0; rank < this->nchoosek; rank++ )
	{
		for ( int i = 0; i < this->k; i++ )
		{
			revset[i] = (this->n - 1) - set[(this->k - 1) - i];
		}

		int revrank = lexIndexOfSet(this->n, this->k, revset);
		this->right_shift[rank] = this->left_shift[revrank];

		getLexSuccessor(this->n, this->k, set);
	}

	free(set);
	free(revset);

	/* the values will be calculated as required */
	return this->nchoosek * sizeof(int);
}

void MMSBranchingManager::freeLeftShift()
{
	if ( this->left_shift != 0 )
	{
		free(this->left_shift);
		this->left_shift = 0;
	}
}

void MMSBranchingManager::freeRightShift()
{
	if ( this->right_shift != 0 )
	{
		free(this->right_shift);
		this->right_shift = 0;
	}
}

int MMSBranchingManager::getLeftShift( int* set, int size, int rank )
{
	if ( size == this->k && this->left_shift[rank] > 0 )
	{
		return this->left_shift[rank];
	}

	if ( size == 1 )
	{
		return set[0] + 1;
	}

	if ( size < 1 )
	{
		return 1;
	}

	/* Shouldn't I be storing these things? */
	int value = nChooseK(set[size - 1] + 1, size) - nChooseK(set[size - 1] - set[0], size);
	for ( int j = 1; j <= size - 2; j++ )
	{
		//int jrank = lexIndexOfSet(this->n, j, set); // until we actually store the smaller values, this isn't needed.
		value -= (this->getLeftShift(set, j, 0) * nChooseK(set[size - 1] - set[j], size - j));
	}

	if ( size == this->k )
	{
		this->left_shift[rank] = value;
		this->left_c_star[0][rank] = value;
	}

	return value;
}

int MMSBranchingManager::getRightShift( int* set, int rank )
{
	if ( this->right_shift[rank] > 0 )
	{
		return this->right_shift[rank];
	}

	int* revset = (int*) malloc(this->k * sizeof(int));
	for ( int i = 0; i < this->k; i++ )
	{
		revset[i] = (this->n - 1) - set[(this->k - 1) - i];
	}
	int revrank = lexIndexOfSet(this->n, this->k, revset);

	this->right_shift[rank] = this->getLeftShift(revset, this->k, revrank);
	this->right_c_star[0][rank] = this->right_shift[rank];
	free(revset);

	return this->right_shift[rank];
}

int MMSBranchingManager::initCStarValues()
{
	this->cur_depth_left_c_star = 0;
	this->max_depth_left_c_star = 100;
	this->left_c_star = (int**) malloc(this->max_depth_left_c_star * sizeof(int*));
	this->left_c_star_updated_to = (int*) malloc(this->max_depth_left_c_star * sizeof(int));
	bzero(this->left_c_star_updated_to, this->max_depth_left_c_star * sizeof(int));

	this->cur_depth_right_c_star = 0;
	this->max_depth_right_c_star = 100;
	this->right_c_star = (int**) malloc(this->max_depth_right_c_star * sizeof(int**));
	this->right_c_star_updated_to = (int*) malloc(this->max_depth_right_c_star * sizeof(int));
	bzero(this->right_c_star_updated_to, this->max_depth_right_c_star * sizeof(int));

	for ( int d = 0; d < this->max_depth_left_c_star; d++ )
	{
		this->left_c_star[d] = (int*) malloc(this->nchoosek * sizeof(int));
		this->right_c_star[d] = (int*) malloc(this->nchoosek * sizeof(int));

		if ( d == 0 )
		{
			for ( int i = 0; i < this->nchoosek; i++ )
			{
				this->left_c_star[d][i] = this->left_shift[i];
				this->right_c_star[d][i] = this->right_shift[i];
			}
		}
		else
		{
			/* to be computed later */
			bzero(this->left_c_star[d], this->nchoosek * sizeof(int));
			bzero(this->right_c_star[d], this->nchoosek * sizeof(int));
		}
	}

	this->stack_left_c_star_depth.push(0);
	this->stack_right_c_star_depth.push(0);

	//free(set);
	this->recalculating_c_star_values = false;

	return this->nchoosek * sizeof(int) * 2 * this->max_depth_left_c_star;
}

void MMSBranchingManager::freeCStarValues()
{
	if ( this->left_c_star != 0 )
	{
		for ( int d = 0; d < this->max_depth_left_c_star; d++ )
		{
			if ( this->left_c_star[d] != 0 )
			{
				free(this->left_c_star[d]);
				this->left_c_star[d] = 0;
			}
		}

		free(this->left_c_star);
		this->left_c_star = 0;
	}
	if ( this->right_c_star != 0 )
	{
		for ( int d = 0; d < this->max_depth_right_c_star; d++ )
		{
			if ( this->right_c_star[d] != 0 )
			{
				free(this->right_c_star[d]);
				this->right_c_star[d] = 0;
			}
		}
		free(this->right_c_star);
		this->right_c_star = 0;
	}

	if ( this->right_c_star_updated_to != 0 )
	{
		free(this->right_c_star_updated_to);
		this->right_c_star_updated_to = 0;
	}
	if ( this->left_c_star_updated_to != 0 )
	{
		free(this->left_c_star_updated_to);
		this->left_c_star_updated_to = 0;
	}
}

void MMSBranchingManager::clearCStarValues()
{
	/* reset the values to initial */
	for ( int i = 0; i < this->nchoosek; i++ )
	{
		this->left_c_star[0][i] = this->left_shift[i];
		this->right_c_star[0][i] = this->right_shift[i];
	}

	while ( this->stack_left_c_star_depth.size() > 0 )
	{
		this->stack_left_c_star_depth.pop();
	}
	this->stack_left_c_star_depth.push(0);

	while ( this->stack_left_c_star_depth.size() > 0 )
	{
		this->stack_left_c_star_depth.pop();
	}

	this->cur_depth_left_c_star = 0;
	this->cur_depth_right_c_star = 0;
}

void MMSBranchingManager::recalculateAllCStarValues( int mode )
{
	this->recalculating_c_star_values = true;

	if ( mode <= 0 )
	{
		/* CALCULATE TO THE LEFT */
		if ( this->cur_depth_left_c_star + 1 >= this->max_depth_left_c_star )
		{
			this->max_depth_left_c_star = this->max_depth_left_c_star + 50;
			this->left_c_star = (int**) realloc(this->left_c_star, this->max_depth_left_c_star * sizeof(int**));
			this->left_c_star_updated_to = (int*) realloc(this->left_c_star_updated_to,
			                                              this->max_depth_left_c_star * sizeof(int*));

			for ( int i = this->cur_depth_left_c_star + 1; i < this->max_depth_left_c_star; i++ )
			{
				this->left_c_star[i] = (int*) malloc(this->nchoosek * sizeof(int));
			}
		}
		/* clear to zero for a fresh start! */
		bzero(this->left_c_star[this->cur_depth_left_c_star + 1], this->nchoosek * sizeof(int));

		if ( this->use_breadth_first_search_left )
		{
			this->recalculateLeftCStarBFS();
			this->cur_depth_left_c_star = this->cur_depth_left_c_star + 1;
			this->left_c_star_updated_to[this->cur_depth_left_c_star] = this->num_genpos;
		}
		else if ( this->use_inclusion_exclusion_left )
		{
			this->recalculateLeftCStarInclusionExclusion();
			this->cur_depth_left_c_star = this->cur_depth_left_c_star + 1;
			this->left_c_star_updated_to[this->cur_depth_left_c_star] = this->num_genpos;
		}
		else if ( this->use_incremental_method_left )
		{
			int* set = (int*) malloc(this->k * sizeof(int));

			/* this loop usually has at most one iteration */
			/* only time it has multiple is if we ran positive propagation */
			for ( int i = this->left_c_star_updated_to[this->cur_depth_left_c_star]; i < this->num_genpos; i++ )
			{
				int rank = this->genpos[i];
				lexIndexToSet(this->n, this->k, rank, set);

				this->recalculateLeftCStarIterative(set);

				this->cur_depth_left_c_star = this->cur_depth_left_c_star + 1;
				this->left_c_star_updated_to[this->cur_depth_left_c_star] = i + 1;

				if ( this->cur_depth_left_c_star + 1 >= this->max_depth_left_c_star )
				{
					this->max_depth_left_c_star = this->max_depth_left_c_star + 50;
					this->left_c_star = (int**) realloc(this->left_c_star, this->max_depth_left_c_star * sizeof(int**));
					this->left_c_star_updated_to = (int*) realloc(this->left_c_star_updated_to,
					                                              this->max_depth_left_c_star * sizeof(int*));

					for ( int i = this->cur_depth_left_c_star + 1; i < this->max_depth_left_c_star; i++ )
					{
						this->left_c_star[i] = (int*) malloc(this->nchoosek * sizeof(int));
					}
				}
				/* clear to zero for a fresh start! */
				bzero(this->left_c_star[this->cur_depth_left_c_star + 1], this->nchoosek * sizeof(int));
			}
			free(set);
		}
	}

	if ( mode >= 0 && this->use_right_c_star )
	{
		if ( this->cur_depth_right_c_star + 1 >= this->max_depth_right_c_star )
		{
			this->max_depth_right_c_star = this->max_depth_right_c_star + 50;
			this->right_c_star = (int**) realloc(this->right_c_star, this->max_depth_right_c_star * sizeof(int**));
			this->right_c_star_updated_to = (int*) realloc(this->right_c_star_updated_to,
			                                               this->max_depth_right_c_star * sizeof(int*));

			for ( int i = this->cur_depth_right_c_star + 1; i < this->max_depth_right_c_star; i++ )
			{
				this->right_c_star[i] = (int*) malloc(this->nchoosek * sizeof(int));
			}
		}
		/* clear to zero for a fresh start! */
		bzero(this->right_c_star[this->cur_depth_right_c_star + 1], this->nchoosek * sizeof(int));

		if ( this->use_breadth_first_search_right )
		{
			this->recalculateRightCStarBFS();
			this->cur_depth_right_c_star = this->cur_depth_right_c_star + 1;
			this->right_c_star_updated_to[this->cur_depth_right_c_star] = this->num_genneg;
		}
		else if ( this->use_inclusion_exclusion_right )
		{
			this->recalculateRightCStarInclusionExclusion();
			this->cur_depth_right_c_star = this->cur_depth_right_c_star + 1;
			this->right_c_star_updated_to[this->cur_depth_right_c_star] = this->num_genneg;
		}
		else if ( this->use_incremental_method_right )
		{
			int* set = (int*) malloc(this->k * sizeof(int));

			/* this loop will probably have a LOT of iterations... */
			/* NOT RECOMMENDED while using 'use_right_c_star' */
			for ( int i = this->right_c_star_updated_to[this->cur_depth_right_c_star]; i < this->num_genneg; i++ )
			{
				int rank = this->genneg[i];
				lexIndexToSet(this->n, this->k, rank, set);

				this->recalculateRightCStarIterative(set);

				this->cur_depth_right_c_star = this->cur_depth_right_c_star + 1;
				this->right_c_star_updated_to[this->cur_depth_right_c_star] = i + 1; // go one beyond

				if ( this->cur_depth_right_c_star + 1 >= this->max_depth_right_c_star )
				{
					this->max_depth_right_c_star = this->max_depth_right_c_star + 50;
					this->right_c_star = (int**) realloc(this->right_c_star,
					                                     this->max_depth_right_c_star * sizeof(int**));
					this->right_c_star_updated_to = (int*) realloc(this->right_c_star_updated_to,
					                                               this->max_depth_right_c_star * sizeof(int*));

					for ( int i = this->cur_depth_right_c_star + 1; i < this->max_depth_right_c_star; i++ )
					{
						this->right_c_star[i] = (int*) malloc(this->nchoosek * sizeof(int));
					}
				}
				/* clear to zero for a fresh start! */
				bzero(this->right_c_star[this->cur_depth_right_c_star + 1], this->nchoosek * sizeof(int));
			}
			free(set);
		}
	}

	this->recalculating_c_star_values = false;
}

void MMSBranchingManager::recalculateRightCStarBFS()
{
	for ( int rank = this->max_c_star_rank - 1; rank >= this->min_c_star_rank; rank-- )
	{
		if ( this->set_labels[rank] != LABEL_C_STAR )
		{
			continue;
		}

		this->clearMarkers();

		/* STORE VALUE */
		this->right_c_star[this->cur_depth_right_c_star + 1][rank] = this->recalculateRightCStarValueRecurse(rank);
	}
}

int MMSBranchingManager::recalculateRightCStarValueRecurse( int rank )
{
	if ( this->set_labels[rank] != LABEL_C_STAR || this->markers[rank] != 0 )
	{
		/* not in Cstar OR already marked! */
		return 0;
	}

	this->markers[rank] = 1;

	int value = 1;
	for ( int i = 0; i < this->right_degree[rank]; i++ )
	{
		int rightrank = this->right_neighbors[this->k * rank + i];

		int result = this->recalculateRightCStarValueRecurse(rightrank);
		value += result;
	}

	return value;
}

void MMSBranchingManager::recalculateRightCStarInclusionExclusion()
{
	for ( int rank = this->max_c_star_rank - 1; rank >= this->min_c_star_rank; rank-- )
	{
		if ( this->set_labels[rank] != LABEL_C_STAR )
		{
			continue;
		}

		/* RECOMPUTE! */
		int value = 1;
		int rightdeg = this->right_degree[rank];
		int* nset = (int*) malloc(this->k * sizeof(int));
		int* nset_vals = (int*) malloc(this->k * sizeof(int));

		/* INCLUSION-EXCLUSION */
		/* start with single-sized sets */
		for ( int j = 0; j < rightdeg; j++ )
		{
			int neighbor_rank = this->right_neighbors[this->k * rank + j];
			int tval = this->getRightCStarValue(neighbor_rank);
			value += tval;

		}

		for ( int i = 2; i <= this->right_degree[rank]; i++ )
		{
			int num_sets = nChooseK(this->right_degree[rank], i);

			lexIndexToSet(this->right_degree[rank], i, 0, nset);
			for ( int nrank = 0; nrank < num_sets; nrank++ )
			{
				bool all_in_c_star = true;
				for ( int j = 0; j < i; j++ )
				{
					nset_vals[j] = this->right_neighbors[this->k * rank + nset[j]];
					if ( this->set_labels[nset_vals[j]] != LABEL_C_STAR )
					{
						all_in_c_star = false;
					}
				}

				/* if any not in C*, then we definitely have 0 for this value */
				if ( all_in_c_star )
				{
					/* find the MEET of these sets */
					int meetrank = this->getMeetRank(i, nset_vals);

					int tval = this->getRightCStarValue(meetrank);

					if ( (i % 2) == 0 )
					{
						// Even sets were over counted
						value -= tval;
					}
					else
					{
						// Odd sets were under counted
						value += tval;
					}
				}

				/* get next subset of neighbors */
				getLexSuccessor(this->right_degree[rank], i, nset);
			}
		}

		free(nset);
		free(nset_vals);

		/* STORE VALUE */
		this->markers[rank] = 1;

		this->right_c_star[this->cur_depth_right_c_star + 1][rank] = value;
	}
}

void MMSBranchingManager::recalculateRightCStarIterative( int* new_set )
{
//	printf("updating Rstar for set ");
//	for ( int i = 0; i < this->k; i++ )
//	{
//		printf("%2d ", new_set[i]);
//	}
//	printf("\n");

	/* go through our sets, join with new_set, and subtract previous right_c_star value. */
	int* set = (int*) malloc(this->k * sizeof(int));
	int* meet_set = (int*) malloc(this->k * sizeof(int));

	int rank = this->max_c_star_rank - 1;
	lexIndexToSet(this->n, this->k, rank, set);

	int j = 0;
	do
	{
		if ( this->set_labels[rank] == LABEL_C_STAR )
		{
			/* compute join and rank, starting where set changed */
			for ( int i = j; i < this->k; i++ )
			{
				meet_set[i] = set[i];

				if ( new_set[i] > set[i] )
				{
					meet_set[i] = new_set[i];
				}
			}

			int meet_rank = lexIndexOfSet(this->n, this->k, meet_set);
			int meet_rstar = this->right_c_star[this->cur_depth_right_c_star][meet_rank];

			// here's where the magic happens
			this->right_c_star[this->cur_depth_right_c_star + 1][rank] =
			        this->left_c_star[this->cur_depth_right_c_star][rank] - meet_rstar;
		}

		j = getLexPredecessor(this->n, this->k, set);
		rank--;
	}
	while ( rank >= this->min_c_star_rank );

	free(set);
	free(meet_set);
}

void MMSBranchingManager::recalculateLeftCStarBFS()
{
	for ( int rank = this->min_c_star_rank; rank < this->max_c_star_rank; rank++ )
	{
		if ( this->set_labels[rank] != LABEL_C_STAR )
		{
			continue;
		}

		this->clearMarkers();

		/* STORE VALUE */
		this->left_c_star[this->cur_depth_left_c_star + 1][rank] = this->recalculateLeftCStarValueRecurse(rank);
	}
}

int MMSBranchingManager::recalculateLeftCStarValueRecurse( int rank )
{
	if ( this->set_labels[rank] != LABEL_C_STAR || this->markers[rank] != 0 )
	{
		/* not in Cstar OR already marked! */
		return 0;
	}

	this->markers[rank] = 1;

	int value = 1;
	for ( int i = 0; i < this->left_degree[rank]; i++ )
	{
		int leftrank = this->left_neighbors[this->k * rank + i];

		int result = this->recalculateLeftCStarValueRecurse(leftrank);
		value += result;
	}

	return value;
}

void MMSBranchingManager::recalculateLeftCStarInclusionExclusion()
{
	for ( int rank = this->min_c_star_rank; rank < this->max_c_star_rank; rank++ )
	{
		if ( this->set_labels[rank] != LABEL_C_STAR )
		{
			continue;
		}

		/* RECOMPUTE! */
		int value = 1;
		int* nset = (int*) malloc(this->k * sizeof(int));
		int* nset_vals = (int*) malloc(this->k * sizeof(int));

		/* INCLUSION-EXCLUSION */
		/* start with single-sized neighbor sets */
		int leftdeg = this->left_degree[rank];
		for ( int j = 0; j < leftdeg; j++ )
		{
			int neighbor_rank = this->left_neighbors[this->k * rank + j];
			int tval = this->getLeftCStarValue(neighbor_rank);
			value += tval;
		}

		/* correct for over/under counting */
		for ( int i = 2; i <= this->left_degree[rank]; i++ )
		{
			int num_sets = nChooseK(this->left_degree[rank], i);

			lexIndexToSet(this->left_degree[rank], i, 0, nset);
			for ( int nrank = 0; nrank < num_sets; nrank++ )
			{
				bool all_in_c_star = true;
				for ( int j = 0; j < i; j++ )
				{
					nset_vals[j] = this->left_neighbors[this->k * rank + nset[j]];
					if ( this->set_labels[nset_vals[j]] != LABEL_C_STAR )
					{
						all_in_c_star = false;
					}
				}

				/* if any not in C*, then we definitely have 0 for this value */
				if ( all_in_c_star )
				{
					/* find the JOIN of these sets */
					int joinrank = this->getJoinRank(i, nset_vals);

					if ( this->set_labels[joinrank] == LABEL_C_STAR )
					{

						int tval = this->getLeftCStarValue(joinrank);

						if ( (i % 2) == 0 )
						{
							// Even sets were over counted
							value -= tval;
						}
						else
						{
							// Odd sets were under counted
							value += tval;
						}
					}
				}

				getLexSuccessor(this->left_degree[rank], i, nset);
			}
		}

		free(nset);
		free(nset_vals);

		/* STORE VALUE */
		this->markers[rank] = 1;
		this->left_c_star[this->cur_depth_left_c_star + 1][rank] = value;
	}
}

void MMSBranchingManager::recalculateLeftCStarIterative( int* new_set )
{
	/* go through our sets, join with new_set, and subtract previous left_c_star value. */
//	printf("updating Lstar for set ");
//	for ( int i = 0; i < this->k; i++ )
//	{
//		printf("%2d ", new_set[i]);
//	}
//	printf("\n");
	int* set = (int*) malloc(this->k * sizeof(int));
	int* join_set = (int*) malloc(this->k * sizeof(int));

	int rank = this->min_c_star_rank;
	lexIndexToSet(this->n, this->k, rank, set);

	int j = 0;
	do
	{
		if ( this->set_labels[rank] == LABEL_C_STAR )
		{
			/* compute join and rank, starting where set changed */
			/* i = j ? */
			for ( int i = 0; i < this->k; i++ )
			{
				join_set[i] = set[i];

				if ( new_set[i] < set[i] )
				{
					join_set[i] = new_set[i];
				}
			}

			int join_rank = lexIndexOfSet(this->n, this->k, join_set);
			int join_lstar = this->left_c_star[this->cur_depth_left_c_star][join_rank];

			// here's where the magic happens
			this->left_c_star[this->cur_depth_left_c_star + 1][rank] =
			        this->left_c_star[this->cur_depth_left_c_star][rank] - join_lstar;
		}

		j = getLexSuccessor(this->n, this->k, set);
		rank++;
	}
	while ( rank < this->max_c_star_rank );

	free(set);
	free(join_set);
}

int MMSBranchingManager::getLeftCStarValue( int rank )
{
	if ( this->set_labels[rank] != LABEL_C_STAR )
	{
		return 0;
	}

	if ( !(this->recalculating_c_star_values) )
	{
		return this->left_c_star[this->cur_depth_left_c_star][rank];
	}
	else if ( this->markers[rank] != 0 )
	{
		/* if recalculating, use _next_ depth for previously calculated values */
		return this->left_c_star[this->cur_depth_left_c_star + 1][rank];
	}

	printf("-- You are asking for left_c_star of %d, but it hasn't been computed yet!\n", rank);
	return -1;
}

int MMSBranchingManager::getRightCStarValue( int rank )
{
	if ( this->set_labels[rank] != LABEL_C_STAR )
	{
		return 0;
	}

	if ( !(this->recalculating_c_star_values) )
	{
		return this->right_c_star[this->cur_depth_right_c_star][rank];
	}
	else if ( this->markers[rank] != 0 )
	{
		/* if recalculating, use _next_ depth for previously calculated values */
		return this->right_c_star[this->cur_depth_right_c_star + 1][rank];
	}

	printf("-- You have tried to request a value %d that hasn't been computed yet.!", rank);

	return -1;
}

void MMSBranchingManager::snapshotCStarValues()
{
	this->stack_left_c_star_depth.push(this->cur_depth_left_c_star);
	this->stack_right_c_star_depth.push(this->cur_depth_right_c_star);
}

void MMSBranchingManager::rollbackCStarValues()
{
	this->cur_depth_left_c_star = this->stack_left_c_star_depth.top();
	this->stack_left_c_star_depth.pop();

	this->cur_depth_right_c_star = this->stack_right_c_star_depth.top();
	this->stack_right_c_star_depth.pop();
}

void MMSBranchingManager::initGenSets()
{
	this->num_genneg = 0;
	this->max_genneg = 5000;
	this->genneg = (int*) malloc(this->max_genneg * sizeof(int*));

	this->num_genpos = 0;
	this->max_genpos = 500;
	this->genpos = (int*) malloc(this->max_genneg * sizeof(int*));

	this->stack_num_genneg.push(0);
	this->stack_num_genpos.push(0);
}

void MMSBranchingManager::freeGenSets()
{
	if ( this->genneg != 0 )
	{
		free(this->genneg);
		this->genneg = 0;
	}

	if ( this->genpos != 0 )
	{
		free(this->genpos);
		this->genpos = 0;
	}
}

void MMSBranchingManager::clearGenSets()
{
	this->num_genneg = 0;
	this->num_genpos = 0;

	while ( this->stack_num_genneg.size() > 0 )
	{
		this->stack_num_genneg.pop();
	}

	while ( this->stack_num_genpos.size() > 0 )
	{
		this->stack_num_genpos.pop();
	}
}

void MMSBranchingManager::snapshotGenSets()
{
	this->stack_num_genneg.push(this->num_genneg);
	this->stack_num_genpos.push(this->num_genpos);
}

void MMSBranchingManager::rollbackGenSets()
{
	this->num_genneg = this->stack_num_genneg.top();
	this->stack_num_genneg.pop();

	this->num_genpos = this->stack_num_genpos.top();
	this->stack_num_genpos.pop();
}

void MMSBranchingManager::addNewPositiveSets( int number )
{
	(this->num_forced_positive) += number;
}

void MMSBranchingManager::snapshotNumPos()
{
	this->stackNumPos.push(this->num_forced_positive);
}
void MMSBranchingManager::rollbackNumPos()
{
	this->num_forced_positive = this->stackNumPos.top();
	this->stackNumPos.pop();
}
void MMSBranchingManager::clearNumPos()
{
	this->num_forced_positive = 0;

	while ( this->stackNumPos.size() > 0 )
	{
		this->stackNumPos.pop();
	}
}

int MMSBranchingManager::initSetLabels()
{
	this->set_labels = (char*) malloc(this->nchoosek);
	bzero(this->set_labels, this->nchoosek);

	this->min_c_star_rank = 0;
	this->stack_min_c_star.push(0);

	this->max_c_star_rank = this->nchoosek;
	this->stack_max_c_star.push(this->nchoosek);

	this->stack_set_label_sizes.push(0);

	this->num_in_c_star = this->nchoosek;
	this->stack_num_in_c_star.push(this->nchoosek);

	return this->nchoosek;
}

void MMSBranchingManager::freeSetLabels()
{
	if ( this->set_labels != 0 )
	{
		free(this->set_labels);
		this->set_labels = 0;
	}
}

void MMSBranchingManager::branchPositiveSet( int rank )
{
	this->addMinimalPositiveSet(rank);

	/* we need to recalculate the left c star values */
// NEW FAST WAY:
	this->recalculateAllCStarValues(-1);
}

int MMSBranchingManager::addPositiveSet( int rank )
{
	if ( this->set_labels[rank] == LABEL_C_POS )
	{
		/* already marked positive */
		return 0;
	}
	else if ( this->set_labels[rank] == LABEL_C_NEG )
	{
		/* Big problem! Trying to mark a positive set to the right of a negative set*/
		/* Since we propagate negative sets, this should never happen */
		printf("-- Trying to make %10d positive but it is negative!\n", rank);
		return -1;
	}

	std::queue<int> q;

	q.push(rank);
	this->set_labels[rank] = LABEL_C_POS;

	int number_converted = 0;

	while ( q.size() > 0 )
	{
		int cur_rank = q.front();
		q.pop();

		this->stack_set_labels.push(cur_rank);
		(this->num_in_c_star) = this->num_in_c_star - 1;
		number_converted++;

		/* Walk over all left sets and also label them */
		for ( int i = 0; i < this->left_degree[cur_rank]; i++ )
		{
			int neighbor_rank = this->left_neighbors[this->k * cur_rank + i];

			if ( this->set_labels[neighbor_rank] == LABEL_C_STAR )
			{
				q.push(neighbor_rank);
				this->set_labels[neighbor_rank] = LABEL_C_POS;
			}
		}
	}

	if ( rank == this->min_c_star_rank )
	{
		bool found_min = false;
		for ( int i = rank + 1; i < this->max_c_star_rank; i++ )
		{
			if ( set_labels[i] == LABEL_C_STAR )
			{
				this->min_c_star_rank = i;
				found_min = true;
				break;
			}
		}
		if ( !found_min )
		{
			this->min_c_star_rank = this->max_c_star_rank;
		}
	}

	return number_converted;
}

void MMSBranchingManager::branchNegativeSet( int rank )
{
	this->addMaximalNegativeSet(rank);

	/* we need to recalculate the right c star values for C* ranks */
// NEW FAST WAY
	this->recalculateAllCStarValues(1);
}

int MMSBranchingManager::addNegativeSet( int rank )
{
	if ( this->set_labels[rank] == LABEL_C_NEG )
	{
		/* already marked negative */
		return 0;
	}
	else if ( this->set_labels[rank] == LABEL_C_POS )
	{
		/* Big problem! Trying to mark a negative set to the left of a positive set*/
		/* Since we propagate positive sets, this should never happen */
		printf("-- Trying to make %10d negative but it is positive!\n\t", rank);

		int* set = (int*) malloc(this->k * sizeof(int));
		lexIndexToSet(n, this->k, rank, set);
		for ( int j = 0; j < this->k; j++ )
		{
			printf("%2d ", set[j]);
		}
		printf("\n");
		free(set);
		return -1;
	}

	std::queue<int> q;

	q.push(rank);
	this->set_labels[rank] = LABEL_C_NEG;

	int number_converted = 0;

	while ( q.size() > 0 )
	{
		int cur_rank = q.front();
		q.pop();

		this->stack_set_labels.push(cur_rank);
		(this->num_in_c_star) = this->num_in_c_star - 1;
		number_converted++;

		/* Walk over all right sets and also label them */
		for ( int i = 0; i < this->right_degree[cur_rank]; i++ )
		{
			int neighbor_rank = this->right_neighbors[this->k * cur_rank + i];

			if ( this->set_labels[neighbor_rank] == LABEL_C_STAR )
			{
				q.push(neighbor_rank);
				this->set_labels[neighbor_rank] = LABEL_C_NEG;
			}
		}
	}

	if ( rank == this->max_c_star_rank - 1 )
	{
		bool found_max = false;
		for ( int i = rank - 1; i >= min_c_star_rank; i-- )
		{
			if ( set_labels[i] == LABEL_C_STAR )
			{
				this->max_c_star_rank = i + 1;
				found_max = true;
				break;
			}
		}

		if ( !found_max )
		{
			/* C^* is empty! */
			this->max_c_star_rank = min_c_star_rank;
		}
	}

	return number_converted;
}

int MMSBranchingManager::addMinimalPositiveSet( int rank )
{
	int num_added = this->addPositiveSet(rank);
	this->addNewPositiveSets(num_added);
	this->addPositiveConstraint(rank);

	/* add to genpos */
	if ( this->num_genpos >= this->max_genpos )
	{
		this->max_genpos = this->max_genpos + 100;
		this->genpos = (int*) realloc(this->genpos, this->max_genpos * sizeof(int));
	}
	this->genpos[this->num_genpos] = rank;
	this->num_genpos = this->num_genpos + 1;

	return num_added;
}

int MMSBranchingManager::addMaximalNegativeSet( int rank )
{
	int num_added = this->addNegativeSet(rank);
	this->addNegativeConstraint(rank);

	/* add to genneg */
	if ( this->num_genneg >= this->max_genneg )
	{
		this->max_genneg = this->max_genneg + 1000;
		this->genneg = (int*) realloc(this->genneg, this->max_genneg * sizeof(int));
	}
	this->genneg[this->num_genneg] = rank;
	this->num_genneg = this->num_genneg + 1;

	return num_added;
}

void MMSBranchingManager::snapshotSetLabels()
{
	this->stack_set_label_sizes.push(this->stack_set_labels.size());
	this->stack_min_c_star.push(this->min_c_star_rank);
	this->stack_max_c_star.push(this->max_c_star_rank);
	this->stack_num_in_c_star.push(this->num_in_c_star);
}

void MMSBranchingManager::rollbackSetLabels()
{
	int to_size = this->stack_set_label_sizes.top();

	int num_reset = 0;
	while ( this->stack_set_labels.size() > to_size )
	{
		num_reset++;
		int rank = this->stack_set_labels.top();

		this->set_labels[rank] = LABEL_C_STAR;

		this->stack_set_labels.pop();
	}

//	printf("-- reset %d sets to C*\n", num_reset);
	this->stack_set_label_sizes.pop();

	this->min_c_star_rank = this->stack_min_c_star.top();
	this->stack_min_c_star.pop();

	this->max_c_star_rank = this->stack_max_c_star.top();
	this->stack_max_c_star.pop();

	this->num_in_c_star = this->stack_num_in_c_star.top();
	this->stack_num_in_c_star.pop();
}

void MMSBranchingManager::clearSetLabels()
{
	bzero(this->set_labels, this->nchoosek);
	this->min_c_star_rank = 0;
	this->max_c_star_rank = this->nchoosek;
	this->num_in_c_star = this->nchoosek;

	/*  clear stacks */
	while ( this->stack_set_labels.size() > 0 )
	{
		this->stack_set_labels.pop();
	}

	while ( this->stack_set_label_sizes.size() > 0 )
	{
		this->stack_set_label_sizes.pop();
	}

	this->stack_set_label_sizes.push(0);

	while ( this->stack_min_c_star.size() > 0 )
	{
		this->stack_min_c_star.pop();
	}
	this->stack_min_c_star.push(0);

	while ( this->stack_max_c_star.size() > 0 )
	{
		this->stack_max_c_star.pop();
	}
	this->stack_max_c_star.push(this->nchoosek);

	while ( this->stack_num_in_c_star.size() > 0 )
	{
		this->stack_num_in_c_star.pop();
	}
	this->stack_num_in_c_star.push(this->nchoosek);
}

void MMSBranchingManager::initAll()
{
	this->found_solution = false;
	this->found_solution_at = -1;
	this->broke_optimal = 0;

	int data_allocated = 0;
	data_allocated += this->initMarkers();
	data_allocated += this->initSetLabels();
	data_allocated += this->initNeighbors();
	data_allocated += this->initLeftShift();
	data_allocated += this->initRightShift();
	data_allocated += this->initCStarValues();
	this->num_forced_positive = 0;
	this->initGenSets();

	this->initLP();

	if ( this->root != 0 )
	{
		delete this->root;
	}
	this->root = new MMSNode(0);
}

void MMSBranchingManager::freeAll()
{
	this->freeMarkers();
	this->freeNeighbors();
	this->freeLeftShift();
	this->freeRightShift();
	this->freeCStarValues();
	this->freeGenSets();
	this->num_forced_positive = 0;
	this->freeSetLabels();

	if ( this->lp_disabled == false )
	{
		this->freeLP();
	}
}

void MMSBranchingManager::snapshot()
{
	this->snapshotCStarValues();
	this->snapshotNumPos();
	this->snapshotSetLabels();
	this->snapshotLP();
	this->snapshotGenSets();
}

void MMSBranchingManager::rollback()
{
	this->rollbackCStarValues();
	this->rollbackNumPos();
	this->rollbackSetLabels();
	this->rollbackLP();
	this->rollbackGenSets();
	this->found_solution = false;
}

int MMSBranchingManager::getBranchSet()
{
	switch ( this->branching_rule )
	{
		case BRANCH_MIDDLE:
			return this->getBranchSetMiddle();

		case BRANCH_BALANCED:
			return this->getBranchSetBalanced();

		case BRANCH_OPTHALF:
			return this->getBranchSetOptHalf();

		case BRANCH_LEFTHALF:
			return this->getBranchSetLeftHalf();

		case BRANCH_LMAXPOS:
			return this->getBranchSetMaxLPositive();

		case BRANCH_LMINNEG:
			return this->getBranchSetMinLNegative();

		default:
			return this->min_c_star_rank;
	}
}

int MMSBranchingManager::getBranchSetMiddle()
{
	int num_sets = 0;
	for ( int i = this->min_c_star_rank; i < this->max_c_star_rank; i++ )
	{
		if ( this->set_labels[i] == LABEL_C_STAR )
		{
			num_sets++;

			if ( num_sets >= (int) (this->num_in_c_star / 2) )
			{
				return i;
			}
		}
	}

	/* we should never get here */
	return -1;
}

int MMSBranchingManager::getBranchSetBalanced()
{
	int cur_index = -1;
	int max_min = 0;
	for ( int i = this->min_c_star_rank; i < this->max_c_star_rank; i++ )
	{
		if ( this->set_labels[i] == LABEL_C_STAR )
		{
			// parameterized balance!
			int min = this->left_c_star[this->cur_depth_left_c_star][i] * this->left_divisor;
			if ( this->right_c_star[this->cur_depth_right_c_star][i] < min )
			{
				min = this->right_c_star[this->cur_depth_right_c_star][i];
			}

			if ( min > max_min )
			{
				max_min = min;
				cur_index = i;
			}
		}
	}

	return cur_index;
}

int MMSBranchingManager::getBranchSetOptHalf()
{
	int cur_index = this->min_c_star_rank;
	double closest_to_half = -1.0;
	double sum = 0.0;
	int* set = (int*) malloc(this->k * sizeof(int));

	for ( int i = this->min_c_star_rank; i < this->max_c_star_rank; i++ )
	{
		if ( this->set_labels[i] == LABEL_C_STAR )
		{
			lexIndexToSet(n, this->k, i, set);

			sum = 0;
			for ( int j = 0; j < this->k; j++ )
			{
				sum += this->solution[set[j]];
			}

			/* distance to -0.5 */
			double dist = abs(sum + 0.5);

			if ( closest_to_half < 0 || dist < closest_to_half )
			{
				closest_to_half = dist;
				cur_index = i;
			}
		}
	}

	if ( closest_to_half < 0.4999 )
	{
		this->broke_optimal = this->broke_optimal + 1;
		//printf("-- getBranchOptHalf makes this inequality infeasible both ways!\n");
	}

	return cur_index;
}

int MMSBranchingManager::getBranchSetMaxLPositive()
{
	int cur_index = -1;
	int cur_lstar = 0;
	double sum = 0.0;
	int* set = (int*) malloc(this->k * sizeof(int));
	lexIndexToSet(this->n, this->k, this->min_c_star_rank, set);

	for ( int i = this->min_c_star_rank; i < this->max_c_star_rank; i++ )
	{
		if ( this->set_labels[i] == LABEL_C_STAR )
		{
			sum = 0;
			for ( int j = 0; j < this->k; j++ )
			{
				sum += this->solution[set[j]];
			}

			if ( sum >= 0 && this->left_c_star[this->cur_depth_left_c_star][i] > cur_lstar )
			{
				cur_lstar = this->left_c_star[this->cur_depth_left_c_star][i];
				cur_index = i;
			}
		}
		getLexSuccessor(this->n, this->k, set);
	}

	free(set);

	if ( cur_index < 0 )
	{
		return this->getBranchSetMinLNegative();
	}

	return cur_index;
}

int MMSBranchingManager::getBranchSetMinLNegative()
{
	int cur_index = -1;
	int cur_lstar = this->nchoosek;
	double sum = 0.0;
	int* set = (int*) malloc(this->k * sizeof(int));
	lexIndexToSet(this->n, this->k, this->min_c_star_rank, set);

	for ( int i = this->min_c_star_rank; i < this->max_c_star_rank; i++ )
	{
		if ( this->set_labels[i] == LABEL_C_STAR )
		{
			sum = 0;
			for ( int j = 0; j < this->k; j++ )
			{
				sum += this->solution[set[j]];
			}

			if ( sum < 0 && this->left_c_star[this->cur_depth_left_c_star][i] < cur_lstar )
			{
				cur_lstar = this->left_c_star[this->cur_depth_left_c_star][i];
				cur_index = i;
			}
		}
		getLexSuccessor(this->n, this->k, set);
	}

	free(set);

	if ( cur_index < 0 )
	{
		return this->getBranchSetMaxLPositive();
	}
	return cur_index;
}

int MMSBranchingManager::getBranchSetLeftHalf()
{
	int max_left = 0;
	for ( int i = this->min_c_star_rank; i < this->max_c_star_rank; i++ )
	{
		if ( this->set_labels[i] == LABEL_C_STAR )
		{
			int val = this->left_c_star[this->cur_depth_left_c_star][i];
			if ( val > max_left )
			{
				max_left = val;
			}
		}
	}

	/* default: max / 2 */
	int goal = (int) (max_left / this->left_divisor);

	int cur_index = -1;
	int cur_val = 0;

	for ( int i = this->min_c_star_rank; i < this->max_c_star_rank; i++ )
	{
		if ( this->set_labels[i] == LABEL_C_STAR )
		{
			int val = this->left_c_star[this->cur_depth_left_c_star][i];
			if ( cur_val < val && val <= goal )
			{
				cur_val = val;
				cur_index = i;
			}
		}
	}

	return cur_index;
}

bool MMSBranchingManager::propagate()
{
	bool result = true;

	bool updated = true;

	switch ( this->propagation_mode )
	{
		case PROPAGATE_SGAC:
			result = this->propagateSGAC();
			break;

		case PROPAGATE_BRANCH_SGAC:
			result = this->propagateBranchSGAC();
			break;

		case PROPAGATE_GAC:
			result = this->propagateGAC();
			break;
	}

	if ( !result || this->num_forced_positive >= this->mms_target )
	{
		/* we have too many positive! */
		return false;
	}

	if ( this->lp_disabled == false )
	{
		result = this->isLPFeasible();
	}

	bool pos_plus_negative = true;
	while ( result && this->propagate_positive && pos_plus_negative )
	{
		pos_plus_negative = false;
		int old_num_positive = this->num_forced_positive;

		result = this->propagatePositive();

		if ( this->num_forced_positive > old_num_positive )
		{
			updated = true;
		}

		if ( !result )
		{
			return false;
		}

		if ( result && updated )
		{
			int cstar_size_before = this->num_in_c_star;
			/* one more stab at the negatives */
			switch ( this->propagation_mode )
			{
				case PROPAGATE_SGAC:
					result = this->propagateSGAC();
					break;

				case PROPAGATE_BRANCH_SGAC:
					result = this->propagateBranchSGAC();
					break;

				case PROPAGATE_GAC:
					result = this->propagateGAC();
					break;
			}

			if ( cstar_size_before > this->num_in_c_star )
			{
				pos_plus_negative = true;
			}
		}

		if ( this->lp_disabled == false )
		{
			result = this->isLPFeasible();
		}
	}

	if ( !result )
	{
		if ( this->write_proof )
		{
			/** TODO: write the full LP? Write the dual? */
			printf("\tThe linear program is infeasible.\t\t%% PROOF\n");
		}

		return false;
	}

	return result;
}

bool MMSBranchingManager::propagateSGAC()
{
	for ( int rank = this->min_c_star_rank; rank < this->max_c_star_rank; rank++ )
	{
		if ( this->set_labels[rank] == LABEL_C_STAR )
		{
			this->snapshot();

			/* try making rank be positive */
			this->branchPositiveSet(rank);

			this->propagateGAC();

			if ( this->isLPFeasible() == false )
			{
				this->rollback();

				this->branchNegativeSet(rank);
			}
			else
			{
				this->rollback();
			}
		}
	}

	return true;
}

bool MMSBranchingManager::propagateBranchSGAC()
{
	bool result = this->propagateGAC();

	if ( !result )
	{
		return false;
	}

	int rank = this->getBranchSet();

	if ( rank < 0 )
	{
		return false;
	}

	this->snapshot();

	bool outproof = this->write_proof;
	this->write_proof = false;

	/* try making rank be positive */
	this->branchPositiveSet(rank);

	result = this->propagateGAC();

	if ( !result || this->isLPFeasible() == false )
	{
		this->rollback();

		this->write_proof = outproof;

		if ( this->write_proof )
		{
			printf("If the sum $ ");
			/* */
			printf("$ is positive, then we find in one GAC-step infeasibility and therefore it must be non-negative.\n");
		}

		this->branchNegativeSet(rank);

		/* recurse! */
		result = this->propagateBranchSGAC();
	}
	else
	{
		this->rollback();
	}

	this->write_proof = outproof;
	return result;
}

bool MMSBranchingManager::propagateGAC()
{
	/* Add negative sets if their positive values become too much */
	bool added_negative = false;

	int num_added = 0;

	added_negative = false;

	for ( int rank = this->min_c_star_rank; rank < this->max_c_star_rank; rank++ )
	{
		if ( this->set_labels[rank] == LABEL_C_STAR )
		{
			if ( this->num_forced_positive == 0 && this->contains0[rank] )
			{
				/* can do Baranyai trick when no positive sets yet! */
				if ( this->left_c_star[this->cur_depth_left_c_star][rank]
				        >= this->mms_target - this->gvals[this->n - this->k] )
				{
					if ( num_added == 0 && this->write_proof )
					{
						printf("\n\nThe following sums must be strictly negative or we have our target number of non-negative sets:\n\\begin{align*}\n\t");
					}

					num_added++;
					/* we can make this set be negative since it contains 0 and the Baranyai technique will get us over the target */
					/* this only works when we haven't forced anything to be positive */
					added_negative = true;

					this->addMaximalNegativeSet(rank);

					if ( this->write_proof )
					{
						if ( num_added > 1 && ((num_added - 1) % 3) == 0 )
						{
							printf("\\\\\n\t");
						}

						int* set = (int*) malloc(this->k * sizeof(int));
						lexIndexToSet(n, this->k, rank, set);
						for ( int i = 0; i < this->k; i++ )
						{
							if ( i == this->k / 2 )
							{
								printf(" & ");
							}
							printf("x_{%2d}", 1 + set[i]);

							if ( i < this->k - 1 )
							{
								printf(" + ");
							}

						}

						free(set);
						if ( (num_added % 3) != 0 )
						{
							printf("\t&\t");
						}
					}
				}
			}
			else if ( this->left_c_star[this->cur_depth_left_c_star][rank] + this->num_forced_positive
			        >= this->mms_target )
			{
				if ( num_added == 0 && this->write_proof )
				{
					printf("\n\nThe following sums must be strictly negative or we have our target number of non-negative sets:\n\\begin{align*}\n\t");
				}

				num_added++;
				added_negative = true;

				this->addMaximalNegativeSet(rank);
				if ( this->write_proof )
				{
					if ( num_added > 1 && ((num_added - 1) % 3) == 0 )
					{
						printf("\\\\\n\t");
					}
					int* set = (int*) malloc(this->k * sizeof(int));
					lexIndexToSet(n, this->k, rank, set);
					for ( int i = 0; i < this->k; i++ )
					{
						if ( i == this->k / 2 )
						{
							printf(" & ");
						}
						printf("x_{%2d}", 1 + set[i]);

						if ( i < this->k - 1 )
						{
							printf(" + ");
						}

					}
					free(set);

					if ( (num_added % 3) != 0 )
					{
						printf("\t&\t");
					}
				}
			}
		}
	}

	/* only recalculate here. */
	if ( added_negative )
	{
		this->recalculateAllCStarValues(1);
	}

	if ( this->write_proof && num_added > 0 )
	{
		printf("\\end{align*}\n\n");
	}

	return true;
}

bool MMSBranchingManager::propagatePositive()
{
	/* Add negative sets if their positive values become too much */
	bool added_positive = false;
	bool ever_added_positive = false;

	int num_added = 0;
	added_positive = false;

	for ( int rank = this->max_c_star_rank - 1; rank >= this->min_c_star_rank; rank-- )
	{
		if ( this->set_labels[rank] == LABEL_C_STAR )
		{
			this->snapshotLP();
			this->addNegativeConstraint(rank);
			bool result = this->isLPFeasible();
			this->rollbackLP();

			if ( !result )
			{
				if ( num_added == 0 && this->write_proof )
				{
					printf("\n\nThe following sums must be non-negative or else the LP becomes infeasible:\n\\begin{align*}\n\t");
				}

				ever_added_positive = true;
				added_positive = true;
				num_added++;

				if ( this->write_proof )
				{
					if ( num_added > 1 && ((num_added - 1) % 3) == 0 )
					{
						printf("\\\\\n\t");
						fflush(stdout);
					}

					int* set = (int*) malloc(this->k * sizeof(int));
					lexIndexToSet(n, this->k, rank, set);
					for ( int i = 0; i < this->k; i++ )
					{
						if ( i == this->k / 2 )
						{
							printf("&");
						}
						printf("x_{%2d}", 1 + set[i]);

						if ( i < this->k - 1 )
						{
							printf(" + ");
						}

					}

					free(set);
					if ( (num_added % 3) != 0 )
					{
						printf("\t&\t");
					}

					fflush(stdout);
				}

				this->addMinimalPositiveSet(rank);

				if ( this->num_forced_positive >= this->mms_target )
				{
					if ( this->write_proof )
					{
						printf("\n\\end{align*}\n");
						printf("These positive sets now force at least %d non-negative $%d$-sums, and our target was $%d$ non-negative $%d$-sums.\n\n",
						       this->num_forced_positive, this->k, this->mms_target, this->k);
					}

					return false;
				}
			}
		}
	}

	if ( this->write_proof && num_added > 0 )
	{
		printf("\n\\end{align*}\nThese positive sets now force at least %d non-negative $%d$-sums.\n\n",
		       this->num_forced_positive, this->k);
	}

	if ( ever_added_positive )
	{
		/* negative means update left c-star values */
		this->recalculateAllCStarValues(-1);
	}

	return true;
}

MMSBranchingManager::MMSBranchingManager( int n, int k, bool write_proof ) :
		SearchManager()
{
	this->n = n;
	this->k = k;
	this->r = n % k;
	this->c = (n - this->r) / k;
	this->nchoosek = nChooseK(n, k);
	this->write_proof = write_proof;
	this->left_divisor = 2.0;

	this->lp_disabled = false;
	this->use_branchless = false;
	this->use_stochastic = false;
	this->done_branchless_or_stochastic = false;

	this->use_breadth_first_search_left = false;
	this->use_inclusion_exclusion_left = false;
	this->use_incremental_method_left = true;
	this->use_breadth_first_search_right = false;
	this->use_inclusion_exclusion_right = true;
	this->use_incremental_method_right = false;

	this->use_right_c_star = true;

	this->gvals = (int*) malloc(n * sizeof(int));

	bzero(this->gvals, n * sizeof(int));
	for ( int i = k; i < n; i++ )
	{
		int m = i - (i % k);
		this->gvals[i] = nChooseK(m - 1, k - 1);
	}
	this->mms_target = nChooseK(n - 1, k - 1);
	this->Baranyai_value = nChooseK((this->c - 1) * k - 1, k - 1);
	this->num_forced_positive = 0;
	this->propagation_mode = PROPAGATE_GAC;
	this->propagate_positive = false;

	this->initAll();
}

/**
 * Default constructor
 */
/**
 * Destructor
 */
MMSBranchingManager::~MMSBranchingManager()
{
	this->freeAll();
	free(this->gvals);
}

/**
 * pushNext -- deepen the search to the next child
 * 	of the current node.
 *
 * @return the label for the new node. -1 if none.
 */
LONG_T MMSBranchingManager::pushNext()
{

	/* get 'top' node */
	SearchNode* parent = 0;
	if ( this->stack.size() > 0 )
	{
		parent = (SearchNode*) this->stack.back();
	}
	else
	{
		parent = (SearchNode*) this->root;
	}

	if ( parent->curChild >= 1 )
	{
		if ( this->write_proof )
		{
			printf("\n\\end{mycases}\n");
		}
		return -1;
	}

	LONG_T child = parent->curChild + 1;

	return this->pushTo(child);
}

/**
 * pushTo -- deepen the search to the specified child
 * 	of the current node.
 *
 * @param child the specified label for the new node
 * @return the label for the new node. -1 if none, or failed.
 */
LONG_T MMSBranchingManager::pushTo( LONG_T child )
{
	/* get 'top' node */
	MMSNode* parent = 0;
	if ( this->stack.size() > 0 )
	{
		parent = (MMSNode*) this->stack.back();
	}
	else
	{
		parent = (MMSNode*) this->root;
	}
	parent->curChild = child;

	this->snapshot();

	// Do Branchless or Stochastic?
	if ( this->searchDepth >= this->jobDepth )
	{
		// time for branchless/stochastic
		if ( this->use_branchless )
		{
			if ( this->done_branchless_or_stochastic )
			{
				return -1;
			}

			bool result = this->doBranchlessSearch();
			this->done_branchless_or_stochastic = true;

			return -1;
		}
		else if ( this->use_stochastic )
		{
			if ( this->done_branchless_or_stochastic )
			{
				return -1;
			}

			bool result = this->doStochasticSearch();
			this->done_branchless_or_stochastic = true;

			return -1;
		}
	}

	if ( parent->branch_rank < 0 )
	{
		parent->branch_rank = this->getBranchSet();

		if ( parent->branch_rank < 0 )
		{
			if ( this->write_proof )
			{
				printf("There are no more sets to select for branching ($\\cutstar = \\varnothing$).\n\n");
				printf(" min_c_star = %d\tmax_c_star = %d\tnum_in_c_star=%d\n\n", this->min_c_star_rank,
				       this->max_c_star_rank, this->num_in_c_star);
			}
			return -1;
		}

		if ( this->write_proof )
		{
			printf("We now consider whether or not the sum $");
			int* set = (int*) malloc(this->k * sizeof(int));
			lexIndexToSet(n, this->k, parent->branch_rank, set);
			for ( int i = 0; i < this->k; i++ )
			{
				printf("x_{%3d}", 1 + set[i]);

				if ( i < this->k - 1 )
				{
					printf(" + ");
				}
			}
			free(set);
			printf("$ is strictly negative or non-negative.\n");
			printf("\\begin{mycases}\n");
		}
	}

	if ( child == 0 )
	{
		if ( this->write_proof )
		{
			printf("\t\\mycase{The sum $");
			int* set = (int*) malloc(this->k * sizeof(int));
			lexIndexToSet(n, this->k, parent->branch_rank, set);
			for ( int i = 0; i < this->k; i++ )
			{
				printf("x_{%3d}", 1 + set[i]);

				if ( i < this->k - 1 )
				{
					printf(" + ");
				}
			}
			free(set);
			printf("$ is strictly negative.}\n");
		}
		this->branchNegativeSet(parent->branch_rank);
	}
	else if ( child == 1 )
	{

		if ( this->write_proof )
		{
			printf("\t\\mycase{The sum $");
			int* set = (int*) malloc(this->k * sizeof(int));
			lexIndexToSet(n, this->k, parent->branch_rank, set);
			for ( int i = 0; i < this->k; i++ )
			{
				printf("x_{%3d}", 1 + set[i]);

				if ( i < this->k - 1 )
				{
					printf(" + ");
				}
			}
			free(set);
			printf("$ is non-negative.}\n");
		}
//		printf("\n-- pushTo: Adding positive set %10d.\n", branch_rank);
		this->branchPositiveSet(parent->branch_rank);

	}
	else
	{
		/* incorrect! */
		this->rollback();
		return -1;
	}

	MMSNode* cnode = new MMSNode(child);

	this->stack.push_back(cnode);

	return child;
}

/**
 * pop -- remove the current node and move up the tree.
 *
 * @return the label of the node after the pop.
 * 		This return value is used for validation purposes
 * 		to check proper implementation of push*() and pop().
 */
LONG_T MMSBranchingManager::pop()
{
	if ( this->stack.size() > 0 )
	{
		SearchNode* node = this->stack.back();
		this->stack.pop_back();
		LONG_T child = node->label;
		delete node;

		if ( this->found_solution_at >= this->searchDepth )
		{
			//printf("-- removing solution at depth %d\n", this->searchDepth);
			this->found_solution_at = -1;
		}

		if ( this->broke_optimal > 0 )
		{
			this->broke_optimal = this->broke_optimal - 1;
		}

		this->rollback();
		return child;
	}
	else
	{
		return -1;
	}
}

/**
 * prune -- Perform a check to see if this node should be pruned.
 *
 * @return 0 if no prune should occur, 1 if prune should occur.
 */
int MMSBranchingManager::prune()
{
	if ( (this->found_solution || this->found_solution_at >= 0) )
	{
		/* we found one already, and aren't using OPTHALF to find many examples */
		printf("--pruning for found_solution\n");
		return 1;
	}

	if ( this->num_forced_positive >= this->mms_target )
	{
		/* we have too many positive! */
		if ( this->write_proof )
		{
			printf("\tWe have %d non-negative sets in $\\genpos_{%d}$, but our target is %d.\t\t\n",
			       this->num_forced_positive, this->searchDepth, this->mms_target);
		}
		printf("--pruning for num_forced_positive\n");
		return 1;
	}

	bool result = this->propagate(); // This also tests the LP

	if ( !result )
	{
		printf("--pruning for propagate.\n");
		return 1;
	}

	return 0;
}

/**
 * isSolution -- Perform a check to see if a solution exists
 * 		at this point.
 *
 *  By "Solution" we of course mean a "counterexample" to the conjecture.
 *  Test the LP feasible point for a counterexample.
 *
 * @return 0 if no solution is found, 1 if a solution is found.
 */
int MMSBranchingManager::isSolution()
{
	bool result = this->isLPOptimalACounterexample();

	if ( result )
	{
		if ( this->found_solution_at < 0 )
		{
			this->found_solution = true;
			this->found_solution_at = this->searchDepth;
			printf("-- found solution at depth %d\n", this->searchDepth);
		}

		return 1;
	}

	return 0;
}

/**
 * writeSolution -- create a buffer that contains a
 * 	description of the solution.
 *
 * Write the counterexample, as well as the number of non-negative k-sets and how that violates expectation
 *
 * @return a string of data allocated with malloc().
 * 	It will be deleted using free().
 */
char* MMSBranchingManager::writeSolution()
{
	this->printLPOptimal();
	printf("T MIN NUMBER_OF_SETS %d VERSUS %d\n", this->sol_pos_sets, nChooseK(this->n - 1, this->k - 1));

	return 0;
}

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
char* MMSBranchingManager::writeStatistics()
{
	/* TODO: collect and write statistics */
	if ( this->found_solution == false )
	{
		printf("T MIN NUMBER_OF_SETS %d VERSUS %d\n", this->mms_target, nChooseK(this->n - 1, this->k - 1));
	}
	else
	{
		printf("T MIN NUMBER_OF_SETS %d VERSUS %d\n", this->sol_pos_sets, nChooseK(this->n - 1, this->k - 1));
	}

	return super::writeStatistics();
}

/**
 * ClearAll will restart the search from scratch, but won't reallocate.
 * It just clears the data to zero.
 */
void MMSBranchingManager::clearAll()
{
	this->clearMarkers();
	this->clearCStarValues();
	this->clearGenSets();
	this->clearNumPos();

	this->found_solution = false;
	this->found_solution_at = -1;
	this->num_forced_positive = 0;

	this->done_branchless_or_stochastic = false;

	this->clearSetLabels();
	this->clearLP();

	if ( this->strong_mode )
	{
		/* we need to add the set 0, n-k+1, ..., n-1 to cut_negative */
		int* set = (int*) malloc(this->k * sizeof(int));

		set[0] = 0;
		for ( int i = 1; i < this->k; i++ )
		{
			set[i] = this->n - this->k + i;
		}
		int index = lexIndexOfSet(this->n, this->k, set);
		this->branchNegativeSet(index);
		this->snapshot();

		free(set);
		set = 0;
	}
}

/**
 * setBranchingRule
 */
void MMSBranchingManager::setBranchingRule( char rule )
{
	this->branching_rule = rule;
}

void MMSBranchingManager::setLeftDivisor( double divisor )
{
	this->left_divisor = divisor;
}

/**
 * setPropagationMode
 */
void MMSBranchingManager::setPropagationMode( char mode )
{
	this->propagation_mode = mode;
}

/**
 * setPropagationPositive
 */
void MMSBranchingManager::setPropagationPositive( bool proppos )
{
	this->propagate_positive = proppos;
}

/**
 * setProofMode
 */
void MMSBranchingManager::setProofMode( bool mode )
{
	this->write_proof = mode;
}

/**
 * setTarget
 *
 * Change the target number of non-negative k-sums to avoid.
 *
 * Useful for computing g(n,k) when it is strictly less than (n-1 choose k-1).
 */
void MMSBranchingManager::setTarget( int target )
{
	this->mms_target = target;
}

/**
 * setGValue
 *
 * Store a previously-computed value of g(m,k) for some other m < n.
 */
void MMSBranchingManager::setGValue( int m, int value )
{
	if ( m < this->k || m >= this->n )
	{
		return;
	}

	this->gvals[m] = value;
}

/**
 * setStrongMode
 *
 * Change whether or not to test for strong pair (n,k).
 */
void MMSBranchingManager::setStrongMode( bool mode )
{
	this->strong_mode = mode;
}

/**
 * useRightCStar()
 *
 * Should we use right c-star values or not?
 */
void MMSBranchingManager::useRightCStar( bool useRCS )
{
	this->use_right_c_star = useRCS;
}

/**
 * useIncrementalMethod()
 *
 * Set us up for the Incremental Method (and appropriately changes the other settings).
 */
void MMSBranchingManager::useIncrementalMethodRight()
{
	this->use_incremental_method_right = true;
	this->use_inclusion_exclusion_right = false;
	this->use_breadth_first_search_right = false;
}/**
 * useIncrementalMethod()
 *
 * Set us up for the Incremental Method (and appropriately changes the other settings).
 */
void MMSBranchingManager::useIncrementalMethodLeft()
{
	this->use_incremental_method_left = true;
	this->use_inclusion_exclusion_left = false;
	this->use_breadth_first_search_left = false;
//this->use_right_c_star = false; // ??? automatic?
}

/**
 * useInclusionExclusion()
 *
 * Set us up for the Inclusion/Exclusion Method (and appropriately changes the other settings).
 */
void MMSBranchingManager::useInclusionExclusionRight()
{
	this->use_incremental_method_right = false;
	this->use_inclusion_exclusion_right = true;
	this->use_breadth_first_search_right = false;
}
void MMSBranchingManager::useInclusionExclusionLeft()
{
	this->use_incremental_method_left = false;
	this->use_inclusion_exclusion_left = true;
	this->use_breadth_first_search_left = false;
}
/**
 * useBFS()
 *
 * Set us up for the BFS Method (and appropriately changes the other settings).
 */
void MMSBranchingManager::useBFSRight()
{
	this->use_incremental_method_right = false;
	this->use_inclusion_exclusion_right = false;
	this->use_breadth_first_search_right = true;
}
void MMSBranchingManager::useBFSLeft()
{
	this->use_incremental_method_left = false;
	this->use_inclusion_exclusion_left = false;
	this->use_breadth_first_search_left = true;
}

bool MMSBranchingManager::testAgainst( MMSBranchingManager* obj1, MMSBranchingManager* obj2 )
{
	bool result = true;

	if ( this->min_c_star_rank != obj1->min_c_star_rank || this->min_c_star_rank != obj2->min_c_star_rank )
	{
		printf("PROBLEM AT MIN_C_STAR: %10d %10d %10d\n", this->min_c_star_rank, obj1->min_c_star_rank,
		       obj2->min_c_star_rank);
	}
	if ( this->max_c_star_rank != obj1->max_c_star_rank || this->max_c_star_rank != obj2->max_c_star_rank )
	{
		printf("PROBLEM AT MIN_C_STAR: %10d %10d %10d\n", this->max_c_star_rank, obj1->max_c_star_rank,
		       obj2->max_c_star_rank);
	}

	int* set = (int*) malloc(this->k * sizeof(int));
	lexIndexToSet(this->n, this->k, 0, set);
	for ( int i = 0; i < this->nchoosek; i++ )
	{
		if ( this->left_c_star[this->cur_depth_left_c_star][i] != obj1->left_c_star[obj1->cur_depth_left_c_star][i]
		        || this->left_c_star[this->cur_depth_left_c_star][i]
		                != obj2->left_c_star[obj2->cur_depth_left_c_star][i] )
		{
			result = false;
			printf("LEFT_C_STAR BAD AT INDEX %5d (SET: ", i);
			for ( int j = 0; j < this->k; j++ )
			{
				printf("%3d ", set[j]);
			}

			printf(") : VALUES %10d %10d %10d (previously %10d; depth = %d)\n",
			       this->left_c_star[this->cur_depth_left_c_star][i], obj1->left_c_star[obj1->cur_depth_left_c_star][i],
			       obj2->left_c_star[obj2->cur_depth_left_c_star][i],
			       obj2->left_c_star[obj2->cur_depth_left_c_star - 1][i], obj2->cur_depth_left_c_star);
		}

		getLexSuccessor(this->n, this->k, set);
	}

	lexIndexToSet(this->n, this->k, 0, set);
	for ( int i = 0; i < this->nchoosek; i++ )
	{
		if ( this->right_c_star[this->cur_depth_right_c_star][i] != obj1->right_c_star[obj1->cur_depth_right_c_star][i]
		        || this->right_c_star[this->cur_depth_right_c_star][i]
		                != obj2->right_c_star[obj2->cur_depth_right_c_star][i] )
		{
			result = false;
			printf("RIGHT_C_STAR BAD AT INDEX %5d (SET: ", i);
			for ( int j = 0; j < this->k; j++ )
			{
				printf("%3d ", set[j]);
			}

			printf(") : VALUES %10d %10d %10d (previously %10d; cur_depth = %d)\n",
			       this->right_c_star[this->cur_depth_right_c_star][i],
			       obj1->right_c_star[obj1->cur_depth_right_c_star][i],
			       obj2->right_c_star[obj2->cur_depth_right_c_star][i],
			       obj2->right_c_star[obj2->cur_depth_right_c_star - 1][i], obj2->cur_depth_right_c_star);
		}

		getLexSuccessor(this->n, this->k, set);
	}
	free(set);
	return result;
}

/**
 * disableLP ONLY FOR TESTING PURPOSES!
 */
void MMSBranchingManager::disableLP()
{
	this->lp_disabled = true;
}

/**
 * doBranchlessSearch
 *
 *
 */
bool MMSBranchingManager::doBranchlessSearch()
{
	time_t start_time = time(NULL);
	clock_t start_clock = clock();

	bool result = true;
	long long int num_samples = 0;
	while ( result && time(NULL) - start_time < this->killtime )
	{
		result = this->propagateGAC();

		if ( result )
		{
			result = this->isLPFeasible();

			if ( !result && this->write_proof )
			{
				printf("Now the LP is infeasible!\n\n");
			}

			if ( this->isLPOptimalACounterexample() )
			{
				this->printLPOptimal();
				printf("T SUM NUM_SAMPLES %lld\n", num_samples);
				printf("T SUM NUM_POSITIVE_SETS %d\n\n", this->num_genpos);
				time_t end_time = time(NULL);
				if ( end_time - start_time > 120 )
				{
					printf("T SUM SEARCH_TIME %ld\n\n", (time(NULL) - start_time));
				}
				else
				{
					printf("T SUM SEARCH_TIME %2.4lf\n\n",
					       ((double) (clock() - start_clock)) / (double) CLOCKS_PER_SEC);
				}
				this->writeCompleteJob(stdout);
				return true;
			}
		}

		if ( !result )
		{
			/* AHA! */
//			printf("We have discovered that (%d,%d) is a GOOD pair!\n\n", this->n, this->k);
			printf("T SUM NUM_SAMPLES %lld\n", num_samples);
			printf("T SUM NUM_POSITIVE_SETS %d\n\n", this->num_genpos);
			time_t end_time = time(NULL);
			if ( end_time - start_time > 120 )
			{
				printf("T SUM SEARCH_TIME %ld\n\n", (time(NULL) - start_time));
			}
			else
			{
				printf("T SUM SEARCH_TIME %2.4lf\n\n", ((double) (clock() - start_clock)) / (double) CLOCKS_PER_SEC);
			}
			this->writeCompleteJob(stdout);
			return false;
		}

		result = this->propagatePositive();
	}
	printf("T SUM NUM_SAMPLES %lld\n", num_samples);
	printf("T SUM NUM_POSITIVE_SETS %d\n\n", this->num_genpos);
	time_t end_time = time(NULL);
	if ( end_time - start_time > 120 )
	{
		printf("T SUM SEARCH_TIME %ld\n\n", (time(NULL) - start_time));
	}
	else
	{
		printf("T SUM SEARCH_TIME %2.4lf\n\n", ((double) (clock() - start_clock)) / (double) CLOCKS_PER_SEC);
	}
	this->writePartialJob(stdout);

	return true;
}

/**
 * doStochasticSearch
 *
 * Do ALL negative propagations, then randomly select sets to check if they make the LP infeasible
 * when negative. Make these sets positive, propagate negatives, and repeat!
 *
 */
bool MMSBranchingManager::doStochasticSearch()
{
	time_t start_time = time(NULL);
	clock_t start_clock = clock();

	this->snapshot();

	srand((unsigned) start_time);

	bool result = true;
	long long int num_samples = 0;
	while ( result )
	{
		result = this->propagateGAC();

		if ( this->num_in_c_star == 0 )
		{
			printf("T SUM NUM_SAMPLES %lld\n", num_samples);
			printf("T SUM NUM_POSITIVE_SETS %d\n\n", this->num_genpos);
			time_t end_time = time(NULL);
			if ( end_time - start_time > 120 )
			{
				printf("T SUM STOCHASTIC_SEARCH_TIME %ld\n\n", (time(NULL) - start_time));
			}
			else
			{
				printf("T SUM STOCHASTIC_SEARCH_TIME %2.4lf\n\n",
				       ((double) (clock() - start_clock)) / (double) CLOCKS_PER_SEC);
			}
			this->writeCompleteJob(stdout);
			this->rollback();
			return false;
		}

		if ( result )
		{
			result = this->isLPFeasible();

			if ( !result && this->write_proof )
			{
				printf("Now the LP is infeasible!\n\n");
			}

			if ( this->isLPOptimalACounterexample() )
			{
				this->printLPOptimal();
				printf("T SUM NUM_SAMPLES %lld\n", num_samples);
				printf("T SUM NUM_POSITIVE_SETS %d\n\n", this->num_genpos);
				time_t end_time = time(NULL);
				if ( end_time - start_time > 120 )
				{
					printf("T SUM STOCHASTIC_SEARCH_TIME %ld\n\n", (time(NULL) - start_time));
				}
				else
				{
					printf("T SUM STOCHASTIC_SEARCH_TIME %2.4lf\n\n",
					       ((double) (clock() - start_clock)) / (double) CLOCKS_PER_SEC);
				}
				this->writeCompleteJob(stdout);
				this->rollback();
				return true;
			}
		}

		if ( !result )
		{
			/* AHA! */
			printf("T SUM NUM_SAMPLES %lld\n", num_samples);
			printf("T SUM NUM_POSITIVE_SETS %d\n\n", this->num_genpos);
			time_t end_time = time(NULL);
			if ( end_time - start_time > 120 )
			{
				printf("T SUM STOCHASTIC_SEARCH_TIME %ld\n\n", (time(NULL) - start_time));
			}
			else
			{
				printf("T SUM STOCHASTIC_SEARCH_TIME %2.4lf\n\n",
				       ((double) (clock() - start_clock)) / (double) CLOCKS_PER_SEC);
			}
			this->writeCompleteJob(stdout);
			this->rollback();
			return false;
		}

		bool found_pos = false;
		while ( !found_pos && (time(NULL) - start_time) < this->killtime )
		{
			int rank = -1;
			num_samples++;

			// we find a random rank...
			int in_c_star_index = rand() % this->num_in_c_star;
			int count = 0;
			for ( int i = this->min_c_star_rank; i < this->max_c_star_rank; i++ )
			{
				if ( this->set_labels[i] == LABEL_C_STAR )
				{
					count++;
					if ( count >= in_c_star_index )
					{
						rank = i;
						break; // out of for loop
					}
				}
			}

			this->snapshotLP();
			this->addNegativeConstraint(rank);
			bool result = this->isLPFeasible();
			this->rollbackLP();

			if ( !result )
			{
				if ( this->write_proof )
				{
					printf("\n\nThe following sum must be non-negative or else the LP becomes infeasible:\n\\begin{align*}\n\t");
					int* set = (int*) malloc(this->k * sizeof(int));
					lexIndexToSet(n, this->k, rank, set);
					for ( int i = 0; i < this->k; i++ )
					{
						if ( i == this->k / 2 )
						{
							printf("&");
						}
						printf("x_{%2d}", 1 + set[i]);

						if ( i < this->k - 1 )
						{
							printf(" + ");
						}

					}
					printf("\n\\end{align*}\n");
					free(set);
				}

				this->addMinimalPositiveSet(rank);
				found_pos = true;

				if ( this->num_forced_positive >= this->mms_target )
				{
					if ( this->write_proof )
					{
						printf("These positive sets now force at least %d non-negative $%d$-sums, and our target was $%d$ non-negative $%d$-sums.\n\n",
						       this->num_forced_positive, this->k, this->mms_target, this->k);
					}

					printf("T SUM NUM_SAMPLES %lld\n", num_samples);
					printf("T SUM NUM_POSITIVE_SETS %d\n\n", this->num_genpos);
					time_t end_time = time(NULL);
					if ( end_time - start_time > 120 )
					{
						printf("T SUM STOCHASTIC_SEARCH_TIME %ld\n\n", (time(NULL) - start_time));
					}
					else
					{
						printf("T SUM STOCHASTIC_SEARCH_TIME %2.4lf\n\n",
						       ((double) (clock() - start_clock)) / (double) CLOCKS_PER_SEC);
					}
					this->writeCompleteJob(stdout);
					this->rollback();
					return false;
				}
				else if ( this->write_proof )
				{
					printf("These positive sets now force at least %d non-negative $%d$-sums.\n\n",
					       this->num_forced_positive, this->k);
				}
			}
		}
	}

	// write it as PARTIAL for use in supercomputer!
	this->writePartialJob(stdout);

	printf("T SUM NUM_SAMPLES %lld\n", num_samples);
	printf("T SUM NUM_POSITIVE_SETS %d\n\n", this->num_genpos);
	time_t end_time = time(NULL);
	if ( end_time - start_time > 120 )
	{
		printf("T SUM STOCHASTIC_SEARCH_TIME %ld\n\n", (time(NULL) - start_time));
	}
	else
	{
		printf("T SUM STOCHASTIC_SEARCH_TIME %2.4lf\n\n", ((double) (clock() - start_clock)) / (double) CLOCKS_PER_SEC);
	}
	this->rollback();
	return true;
}

void MMSBranchingManager::useBranchlessSearch( bool branchless )
{
	this->use_branchless = branchless;
}

void MMSBranchingManager::useStochasticSearch( bool stochastic )
{
	this->use_stochastic = stochastic;
}
