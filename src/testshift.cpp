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
 * testshift.cpp
 *
 * Test the shift calculations in
 * shiftcalculations.cpp
 *
 *  Created on: Nov 27, 2012
 *      Author: stolee
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "translation.hpp"
#include "shiftcalculations.hpp"

int main( int argc, char** argv )
{
	bool do_bfs = false;
	bool do_form = true;
	bool do_lex = true;

	int k = 3;
	int n = 20;
	int N = 50;

	clock_t start_clock;
	clock_t end_clock;

	time_t start_time;
	time_t end_time;

	for ( int i = 0; i < argc; i++ )
	{
		if ( i < argc - 1 )
		{
			if ( strcmp(argv[i], "-k") == 0 )
			{
				k = atoi(argv[i + 1]);
			}
			if ( strcmp(argv[i], "-n") == 0 )
			{
				n = atoi(argv[i + 1]);
			}
			if ( strcmp(argv[i], "-N") == 0 )
			{
				N = atoi(argv[i + 1]);
			}
		}
	}

	initBinomialTable(N, k);
	int NchooseK = nChooseK(N, k);
	int* data_form = (int*) malloc(sizeof(int) * NchooseK);
	int* data_lex = (int*) malloc(sizeof(int) * NchooseK);

	for ( int i = n; i <= N; i++ )
	{
		printf("Starting test for n = %d, k = %d.\n", i, k);
		if ( do_form )
		{
			printf("Testing FormulaMethod...\n");
			start_clock = clock();
			start_time = time(NULL);
			fill_leftshift_values_form(i, k, data_form);
			end_clock = clock();
			end_time = time(NULL);

			double float_time = ((double) (end_clock - start_clock)) / (double) CLOCKS_PER_SEC;
			printf("...FormulaMethod(n=%06d,k=%d) done in %4.3lf seconds.\n", i, k, float_time);
		}
		if ( do_lex )
		{
			printf("Testing Lex Streaming...\n");
			start_clock = clock();
			start_time = time(NULL);
			fill_leftshift_values_lex(i, k, data_lex);
			end_clock = clock();
			end_time = time(NULL);

			double float_time = ((double) (end_clock - start_clock)) / (double) CLOCKS_PER_SEC;
			printf("...LexStreaming(n=%06d,k=%d) done in %4.3lf seconds.\n", i, k, float_time);
		}

		continue;
		/* TIME TO VERIFY CORRECTNESS */
		int ichoosek = nChooseK(i, k);
		int* set = (int*) malloc(sizeof(int) * k);
		for ( int j = 0; j < ichoosek; j++ )
		{
			if ( data_form[j] != data_lex[j] )
			{
				printf("!!! wrong value at index %05d: form = %05d. lex = %05d. Set = ", j, data_form[j], data_lex[j]);
				lexIndexToSet(n, k, j, set);
				for ( int a = 0; a < k; a++ )
				{
					printf("%04d ", set[a]);
				}
				printf("\n");
			}
		}
		free(set);
	}

	free(data_form);
	free(data_lex);

	cleanBinomialTable();
}

