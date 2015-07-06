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
 * shiftcalculations.cpp
 *
 * Compute certain values associated with the shift poset.
 *
 *  Created on: Nov 27, 2012
 *      Author: stolee
 */

#include <stdlib.h>
#include <string.h>
#include "translation.hpp"
#include "shiftcalculations.hpp"

/**
 * fill_leftshift_values_* computes the values L_k(S) for all
 * k-subsets of [n] and places them in the data array by lexicographic
 * rank. The suffix of the function definitions describes the
 * method used:
 * 		- form : use Chowdhury's formula
 * 		- lex  : use the lexicographical streaming method
 */
void fill_leftshift_values_form( int n, int k, int* data )
{
	int** table = (int**) malloc(sizeof(int*) * (k + 1));
	for ( int i = 0; i < k; i++ )
	{
		table[i] = (int*) malloc(sizeof(int) * nChooseK(n, i));
	}
	table[k] = data;

	table[0][0] = 0;

	int* set = (int*) malloc(sizeof(int) * k);

	for ( int i = 1; i <= k; i++ )
	{
		for ( int j = 0; j < i; j++ )
		{
			set[j] = j;
		}

		int index = 0;
		table[i][index] = 1;

		/** IN LEX ORDER! **/
		while ( getLexSuccessor(n, i, set) >= 0 )
		{
			index++;

			int value = nChooseK(set[i - 1] + 1, i) - nChooseK(set[i - 1] - set[0], i);

			for ( int l = 1; l <= i - 2; l++ )
			{
				int lindex = lexIndexOfSet(n, l, set);
				value -= table[l][lindex] * nChooseK(set[i - 1] - set[l], i - l);
			}

			table[i][index] = value;
		}
	}

	for ( int i = 0; i < k; i++ )
	{
		free(table[i]);
	}
	free(table);
	free(set);
}

/**
 * fill_leftshift_values_* computes the values L_k(S) for all
 * k-subsets of [n] and places them in the data array by lexicographic
 * rank. The suffix of the function definitions describes the
 * method used:
 * 		- form : use Chowdhury's formula
 * 		- lex  : use the lexicographical streaming method
 */
void fill_leftshift_values_lex( int n, int k, int* data )
{
	int* set = (int*) malloc(k * sizeof(int));
	int index = 0;
	data[0] = 1;

	int** a = (int**) malloc((k + 1) * sizeof(int*));
	for ( int j = 0; j <= k; j++ )
	{
		a[j] = (int*) malloc((k + 1) * sizeof(int));
		for ( int l = 0; l <= k; l++ )
		{
			a[j][l] = 1;
		}
	}

	int j = 0;
	while ( (j = getLexSuccessor(n, k, set)) >= 0 )
	{
		/* we increased our position */
		index++;

		/* updates! */
		if ( j > 0 )
		{
			for ( int l = j; l <= k; l++ )
			{
				a[l][l] = a[j - 1][l] + a[j - 1][j - 1];
			}
		}
		else
		{
			for ( int l = j; l <= k; l++ )
			{
				a[l][l] = a[0][l] + 1;
			}
		}

		data[index] = a[k][k];

		/* further updates */
		if ( j != 0 && set[j - 1] + 2 == set[j] )
		{
			for ( int l = 0; l <= k; l++ )
			{
				a[j - 1][l] = a[l][l];
			}
		}
	}

	for ( int j = 0; j <= k; j++ )
	{
		free(a[j]);
	}
	free(a);
	free(set);
}

/**
 * update_lstar_values_lex
 *
 * @param n the number of variables.
 * @param k the size of the subsets
 * @param lstar_data The current lstar values, will contain new values.
 * @param cutstar_markers marks when a set is in cutstar or not, will be updated if S join T = S.
 * @param new_set is the new set being added to the positive sets
 */
void update_lstar_values_lex( int n, int k, int* lstar_data, int* cutstar_markers, int* new_set )
{
	int* cur_set = (int*) malloc(sizeof(int) * k);
	int* join_set = (int*) malloc(sizeof(int) * k);
	for ( int i = 0; i < k; i++ )
	{
		cur_set[i] = n - k + i;
		join_set[i] = new_set[i];
	}
	int index = nChooseK(n, k) - 1;
	int j = 0;
	int min_updated_j = k;

	while ( (j = getLexPredecessor(n, k, cur_set)) >= 0 )
	{
		if ( j < min_updated_j )
		{
			min_updated_j = j;
		}
		if ( cutstar_markers[index] == 0 )
		{
			continue;
		}

		/* update join set to match */
		for ( int l = min_updated_j; l < k; l++ )
		{
			if ( new_set[l] <= cur_set[l] )
			{
				join_set[l] = new_set[l];
			}
			else
			{
				join_set[l] = cur_set[l];
			}
		}
		int join_rank = lexIndexOfSet(n, k, join_set);

		if ( join_rank == index )
		{
			cutstar_markers[index] = 1;
		}
		else if ( cutstar_markers[join_rank] == 0 )
		{
			/* since join_rank isn't in cutstar (YET by revlex) we can use lstar_data[join_rank] */
			/* update L* */
			lstar_data[index] = lstar_data[index] - lstar_data[join_rank];
		}
		/* if not above cases, no change to lstar */
	}

	free(cur_set);
	free(join_set);
}

/**
 * fill_rightshift_values_* computes the values R_k(S) for all
 * k-subsets of [n] and places them in the data array by lexicographic
 * rank. The suffix of the function definitions describes the
 * method used:
 * 		- form : use Chowdhury's formula (on the reversed set)
 * (Observe that the lex-order method does not work for R_k(S).)
 */
void fill_rightshift_values_form( int n, int k, int* data )
{
	exit(0);
}

