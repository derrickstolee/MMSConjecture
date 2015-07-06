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
 * mmsconj.cpp
 *
 *  Created on: Jul 24, 2012
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
	bool strong_mode = false;
	bool prop_pos = false;
	bool write_proof = false;
	bool use_right_c_star = true;
	char branch_rule = BRANCH_BALANCED;
	int target = 0;
	double divisor = 2;
	bool do_stochastic = false;
	bool do_branchless = false;

	for ( int i = 1; i < argc; i++ )
	{
		if ( strcmp(argv[i], "--stochastic") == 0 )
		{
			srand((unsigned) time(NULL));
			do_stochastic = true;
//			use_right_c_star = false;
		}
		if ( strcmp(argv[i], "--branchless") == 0 )
		{
			do_branchless = true;
//			use_right_c_star = false;
		}
	}

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
			else if ( strcmp(argv[i + 1], "lmaxpos") == 0 )
			{
				branch_rule = BRANCH_LMAXPOS;
			}
			else if ( strcmp(argv[i + 1], "lminneg") == 0 )
			{
				branch_rule = BRANCH_LMINNEG;
			}

			/* numerical parameter? */
			if ( i < argc - 2 && '0' <= argv[i + 2][0] && argv[i + 2][0] <= '9' )
			{
				divisor = atof(argv[i + 2]);
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

	initBinomialTable(2 * n, k + 1);

	MMSBranchingManager* manager = new MMSBranchingManager(n, k, write_proof);

	manager->importArguments(argc, argv);
	manager->setPropagationMode(prop_mode);
	manager->setPropagationPositive(prop_pos);
	manager->setBranchingRule(branch_rule);
	manager->setLeftDivisor(divisor);
	manager->useBranchlessSearch(do_branchless);
	manager->useStochasticSearch(do_stochastic);

	for ( int i = 0; i < argc - 2; i++ )
	{
		if ( strcmp(argv[i], "--method") == 0 )
		{
			if ( strcmp(argv[i + 1], "bfs") == 0 )
			{
				manager->useBFSLeft();
			}
			else if ( strcmp(argv[i + 1], "inex") == 0 )
			{
				manager->useInclusionExclusionLeft();
			}
			else if ( strcmp(argv[i + 1], "lex") == 0 )
			{
				manager->useIncrementalMethodLeft();
			}

			if ( strcmp(argv[i + 2], "bfs") == 0 )
			{
				manager->useBFSRight();
			}
			else if ( strcmp(argv[i + 2], "inex") == 0 )
			{
				manager->useInclusionExclusionRight();
			}
			else if ( strcmp(argv[i + 2], "lex") == 0 )
			{
				manager->useIncrementalMethodRight();
			}
			else
			{
				/* put NONE here, for example */
				manager->useRightCStar(false);
			}
		}
	}

	if ( target > 0 )
	{
		manager->setTarget(target);
	}
	if ( strong_mode )
	{
		manager->setStrongMode(strong_mode);
	}

	for ( int i = 0; i < argc; i++ )
	{
		if ( argv[i][0] == '-' && argv[i][1] == 'g' && argv[i][2] == 'm' )
		{
			/** we have something -gm#v#*/
			int vpos = 3;
			while ( argv[i][vpos] != 0 && argv[i][vpos] != 'v' )
			{
				vpos++;
			}

			if ( argv[i][vpos] == 'v' )
			{
				// success
				argv[i][vpos] = 0;
				int mval = atoi(argv[i] + 3);
				int val = atoi(argv[i] + vpos + 1);

				manager->setGValue(mval, val);
			}
		}
	}

//	if ( do_branchless )
//	{
//		manager->doBranchlessSearch();
//	}
//	else if ( do_stochastic )
//	{
//		manager->doStochasticSearch();
//	}
//	else
//	{
	while ( manager->readJob(stdin) >= 0 )
	{
		manager->clearAll();

		manager->propagate();

		/** doSearch will load a job AND run branchless/stochastic, if necessary */
		manager->doSearch();
	}
//}

	delete manager;

	cleanBinomialTable();

	return 0;
}

