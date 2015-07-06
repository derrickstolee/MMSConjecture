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
 * CStarTest.cpp
 *
 *  Created on: Jul 26, 2012
 *      Author: stolee
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "translation.hpp"
#include "MMSBranchingManager.hpp"

int main( int argc, char** argv )
{

	int n = 0;
	int k = 0;
	char prop_mode = PROPAGATE_GAC;
	char branch_rule = BRANCH_BALANCED;

	for ( int i = 1; i < argc - 1; i++ )
	{
		if ( strcmp(argv[i], "-N") == 0 )
		{
			n = atoi(argv[i + 1]);
		}
		if ( strcmp(argv[i], "-K") == 0 )
		{
			k = atoi(argv[i + 1]);
		}
		if ( strcmp(argv[i], "--prop") == 0 )
		{
			if ( strcmp(argv[i + 1], "none") == 0 )
			{
				prop_mode = PROPAGATE_NONE;
			}
			else if ( strcmp(argv[i + 1], "gac") == 0 )
			{
				prop_mode = PROPAGATE_GAC;
			}
			else if ( strcmp(argv[i + 1], "bsgac") == 0 )
			{
				prop_mode = PROPAGATE_BRANCH_SGAC;
			}
			else if ( strcmp(argv[i + 1], "sgac") == 0 )
			{
				prop_mode = PROPAGATE_SGAC;
			}
		}
		if ( strcmp(argv[i], "--rule") == 0 )
		{
			if ( strcmp(argv[i + 1], "middle") == 0 )
			{
				branch_rule = BRANCH_MIDDLE;
			}
			else if ( strcmp(argv[i + 1], "balanced") == 0 )
			{
				branch_rule = BRANCH_BALANCED;
			}
			else if ( strcmp(argv[i + 1], "opthalf") == 0 )
			{
				branch_rule = BRANCH_OPTHALF;
			}
		}
	}

	if ( n <= 0 || k <= 0 )
	{
		printf("Usage: CStarTest.exe -N # -K # \n");
		return 0;
	}

	MMSBranchingManager* manager = new MMSBranchingManager(n, k);

	manager->importArguments(argc, argv);
	manager->setPropagationMode(prop_mode);
	manager->setBranchingRule(branch_rule);

	int nchoosek = nChooseK(n, k);

	manager->addPositiveSet(nchoosek / 4);
	manager->addNegativeSet(3 * nchoosek / 4);

	clock_t start_new_clock = clock();
	manager->recalculateAllCStarValues(0);

	int* lvalues = (int*) malloc(nchoosek * sizeof(int));
	int* rvalues = (int*) malloc(nchoosek * sizeof(int));

	for ( int i = 0; i < nchoosek; i++ )
	{
		lvalues[i] = manager->getLeftCStarValue(i);
		rvalues[i] = manager->getRightCStarValue(i);
	}
	clock_t end_new_clock = clock();


	printf("NEW method took %10ld clock cycles.\n",end_new_clock - start_new_clock);
	clock_t start_old_clock = clock();
	for ( int i = 0; i < nchoosek; i++ )
	{
		manager->recalculateLeftCStarValue(i);
		int cur_l = manager->getLeftCStarValue(i);
		if ( cur_l != lvalues[i] )
		{
			printf("-- Left Values do not agree at rank %5d: (new) %10d != %10d (old)\n", i, lvalues[i], cur_l);
		}

		manager->recalculateRightCStarValue(i);
		int cur_r = manager->getRightCStarValue(i);
		if ( cur_r != rvalues[i] )
		{
			printf("--Right Values do not agree at rank %5d: (new) %10d != %10d (old)\n", i, rvalues[i], cur_r);
		}
	}
	clock_t end_old_clock = clock();

	printf("OLD method took %10ld clock cycles.\n",
	       end_old_clock - start_old_clock);

	delete manager;

	return 0;
}

