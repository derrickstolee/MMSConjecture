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

#include "MMSBranchingManager.hpp"

int main( int argc, char** argv )
{

	int n = 0;
	int k = 0;
	char prop_mode = PROPAGATE_GAC;
	char branch_rule = BRANCH_MIDDLE;
	int target = 0;

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
		if ( strcmp(argv[i], "--target") == 0 )
		{
			target = atoi(argv[i + 1]);
		}
	}

	if ( n <= 0 || k <= 0 )
	{
		printf("Usage: mmsconj.exe -N # -K # --prop [none|gac|bsgac|sgac] [TreeSearch arguments]\n");
		return 0;
	}

	MMSBranchingManager* manager = new MMSBranchingManager(n, k);

	manager->importArguments(argc, argv);
	manager->setPropagationMode(prop_mode);
	manager->setBranchingRule(branch_rule);
	manager->setTarget(target);

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

	while ( manager->readJob(stdin) >= 0 )
	{
		manager->propagate();

		manager->doSearch();

		manager->clearAll();
	}

	delete manager;

	return 0;
}

