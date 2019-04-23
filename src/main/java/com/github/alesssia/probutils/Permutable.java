/*
 * 	 Permutable.java
 *
 *   This is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   It is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this code.  If not, see <http://www.gnu.org/licenses/>.
 *
 *   For any bugs or problems found, please contact us at
 *   alessia.visconti@kcl.ac.uk
 */

package com.github.alesssia.probutils;

/**
	Defines an objects that can be subjected to permutation. 
    
	Objects using Adaptive must extend this class.

	@author      Alessia Visconti <alessia.visconti@kcl.ac.uk>
	@version     1.0     
*/

public class Permutable
{
	/**
		Performs a permutation test and returns the p-value
		obtained by the permuted model.
		
		It fails always since this method is not implemented and
		must be overridden by the actual class to use.
		
		@return p-value of the permuted model
		
	*/
     public double permute()
     {
        assert true : "Permutable class do not implement permute. This may be a bug";
		
		//1.1 is just a random value greater than every p-value
        return 1.1;
     }
}