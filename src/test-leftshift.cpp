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
 * test-leftshift.cpp
 *
 *  Created on: Dec 7, 2012
 *      Author: stolee
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "translation.hpp"
#include "MMSBranchingManager.hpp"

//#define TEST_AROUND(OPERATION,MSG) start_time = time(NULL); OPERATION; end_time = time(NULL); printf("%s took %4ld seconds\n", MSG, end_time-start_time);
#define TEST_AROUND(OPERATION,MSG) start_clock = clock(); OPERATION; end_clock = clock(); printf("%s took %4.6lf seconds\n", MSG, ((double)(end_clock - start_clock)) / (double) CLOCKS_PER_SEC);

int main( int argc, char** argv )
{

	int n = 0;
	int k = 0;
	char prop_mode = PROPAGATE_GAC;
	bool strong_mode = false;
	bool prop_pos = false;
	bool write_proof = false;
	char branch_rule = BRANCH_LEFTHALF;
	int target = 0;
	bool use_right_shift = false;

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
		if ( strcmp(argv[i], "--strong") == 0 )
		{
			strong_mode = true;
			if ( target == 0 )
			{
				target = nChooseK(n - 1, k - 1) + 1;
			}
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
			else if ( strcmp(argv[i + 1], "pos") == 0 )
			{
				prop_pos = true;
			}
		}
		if ( strcmp(argv[i], "--proof") == 0 )
		{
			write_proof = true;
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
			else if ( strcmp(argv[i + 1], "lefthalf") == 0 )
			{
				branch_rule = BRANCH_LEFTHALF;
			}
		}
		if ( strcmp(argv[i], "--target") == 0 )
		{
			target = atoi(argv[i + 1]);
		}
	}

	if ( strcmp(argv[argc - 1], "--proof") == 0 )
	{
		write_proof = true;
	}
	if ( strcmp(argv[argc - 1], "--strong") == 0 )
	{
		strong_mode = true;
		if ( target == 0 )
		{
			target = nChooseK(n - 1, k - 1) + 1;
		}
	}

	if ( n <= 0 || k <= 0 )
	{
		printf("Usage: mmsconj.exe -N # -K # --prop [none|gac|bsgac|sgac] [--strong] [TreeSearch arguments]\n");
		return 0;
	}

	initBinomialTable(n, k);

	bool use_bfs = false;

	MMSBranchingManager* manager_bfs = new MMSBranchingManager(n, k, write_proof);
	MMSBranchingManager* manager_inex = new MMSBranchingManager(n, k, write_proof);
	MMSBranchingManager* manager_lex = new MMSBranchingManager(n, k, write_proof);

	MMSBranchingManager* manager = manager_bfs;

	manager->disableLP();
	manager->useRightCStar(use_right_shift);
	manager->importArguments(argc, argv);
	manager->setPropagationMode(prop_mode);
	manager->setPropagationPositive(prop_pos);
	manager->setBranchingRule(branch_rule);

	if ( target > 0 )
	{
		manager->setTarget(target);
	}
	if ( strong_mode )
	{
		manager->setStrongMode(strong_mode);
	}

	manager = manager_inex;

	manager->disableLP();
	manager->useRightCStar(use_right_shift);
	manager->importArguments(argc, argv);
	manager->setPropagationMode(prop_mode);
	manager->setPropagationPositive(prop_pos);
	manager->setBranchingRule(branch_rule);

	if ( target > 0 )
	{
		manager->setTarget(target);
	}
	if ( strong_mode )
	{
		manager->setStrongMode(strong_mode);
	}

	manager = manager_lex;

	manager->disableLP();
	manager->useRightCStar(use_right_shift);
	manager->importArguments(argc, argv);
	manager->setPropagationMode(prop_mode);
	manager->setPropagationPositive(prop_pos);
	manager->setBranchingRule(branch_rule);

	if ( target > 0 )
	{
		manager->setTarget(target);
	}
	if ( strong_mode )
	{
		manager->setStrongMode(strong_mode);
	}

	time_t start_time = 0;
	time_t end_time = 0;
	clock_t start_clock = 0;
	clock_t end_clock = 0;
	int count = 0;

	int* set = (int*) malloc(k * sizeof(int));

	manager_bfs->useBFSLeft();
	manager_inex->useInclusionExclusionLeft();
	manager_lex->useIncrementalMethodLeft();

	/* NOW We make choices and test lstar values! */

	/* TEST */
	if ( (use_bfs && manager_bfs->testAgainst(manager_inex, manager_lex) == false)
	        || (manager_inex->testAgainst(manager_lex, manager_lex) == false) )
	{
		count++;
		printf("FAIL AT %d.\n", count);
		exit(1);
	}

	if ( use_bfs )
	{
		TEST_AROUND(manager_bfs->propagate(), "BFS propagate");
	}
	TEST_AROUND(manager_inex->propagate(), "INEX propagate");
	TEST_AROUND(manager_lex->propagate(), "LEX propagate");

	/* TEST */
	if ( (use_bfs && manager_bfs->testAgainst(manager_inex, manager_lex) == false)
	        || (manager_inex->testAgainst(manager_lex, manager_lex) == false) )
	{
		count++;
		printf("FAIL AT %d.\n", count);
		exit(1);
	}

	int branch;

	for ( int i = 0; i < 10; i++ )
	{
		branch = manager_inex->getBranchSet();
//		printf("Branching at %d : ", branch);
//
//		lexIndexToSet(n, k, branch, set);
//		for ( int j = 0; j < k; j++ )
//		{
//			printf("%3d ", set[j]);
//		}
//		printf("\n");

		if ( use_bfs )
		{
			TEST_AROUND(manager_bfs->branchPositiveSet(branch), "BFS branch");
		}
		TEST_AROUND(manager_inex->branchPositiveSet(branch), "INEX branch");
		TEST_AROUND(manager_lex->branchPositiveSet(branch), " LEX branch");
		printf("\n");

		/* TEST */
		if ( (use_bfs && manager_bfs->testAgainst(manager_inex, manager_lex) == false)
		        || (manager_inex->testAgainst(manager_lex, manager_lex) == false) )
		{
			count++;
			printf("FAIL AT %d.\n", count);
			exit(1);
		}

		if ( use_bfs )
		{
			TEST_AROUND(manager_bfs->propagate(), "BFS propagate");
		}
		TEST_AROUND(manager_inex->propagate(), "INEX propagate");
		TEST_AROUND(manager_lex->propagate(), "LEX propagate");

		/* TEST */
		if ( (use_bfs && manager_bfs->testAgainst(manager_inex, manager_lex) == false)
		        || (manager_inex->testAgainst(manager_lex, manager_lex) == false) )
		{
			count++;
			printf("FAIL AT %d.\n", count);
			exit(1);
		}

	}

	delete manager_bfs;
	delete manager_inex;
	delete manager_lex;

	cleanBinomialTable();

	return 0;
}

