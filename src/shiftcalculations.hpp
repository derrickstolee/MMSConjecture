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
 * shiftcalculations.hpp
 *
 * Compute certain values associated with the shift poset.
 *
 *  Created on: Nov 27, 2012
 *      Author: stolee
 */

#ifndef SHIFTCALCULATIONS_HPP_
#define SHIFTCALCULATIONS_HPP_


/**
 * fill_leftshift_values_* computes the values L_k(S) for all
 * k-subsets of [n] and places them in the data array by lexicographic
 * rank. The suffix of the function definitions describes the
 * method used:
 * 		- form : use Chowdhury's formula
 * 		- lex  : use the lexicographical streaming method
 */
void fill_leftshift_values_form(int n, int k, int* data);
void fill_leftshift_values_lex(int n, int k, int* data);

/**
 * update_lstar_values_lex
 *
 * @param n the number of variables.
 * @param k the size of the subsets
 * @param lstar_data The current lstar values, will contain new values.
 * @param cutstar_markers marks when a set is in cutstar or not, will be updated if S join T = S.
 * @param new_set is the new set being added to the positive sets
 */
void update_lstar_values_lex(int n, int k, int* lstar_data, int* cutstar_markers, int* new_set);

/**
 * fill_rightshift_values_* computes the values R_k(S) for all
 * k-subsets of [n] and places them in the data array by lexicographic
 * rank. The suffix of the function definitions describes the
 * method used:
 * 		- form : use Chowdhury's formula (on the reversed set)
 * (Observe that the lex-order method does not work for R_k(S).)
 */
void fill_rightshift_values_form(int n, int k, int* data);

#endif /* SHIFTCALCULATIONS_HPP_ */
