/*
 * 	 FamilySet.java
 *
 *   PopPAnTe is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   PopPAnTe is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with PopPAnTe.  If not, see <http://www.gnu.org/licenses/>.
 *
 *   For any bugs or problems found, please contact us at
 *   alessia.visconti@kcl.ac.uk
 */

package com.github.alesssia.poppante;

import java.util.*;

/**
	Describe the population in the pedigree file as
	number of people.

	When the analysis mode is "pedcheck", it also 
	describe errors associated with the pedigree file.
	It is an helper for the pedigree file reading.

	@author      Alessia Visconti <alessia.visconti@kcl.ac.uk>
	@version     1.0               
*/

class FamilySet
{
	/**  The number of members in the population */
	int people; 		
	/**  Errors in the pedigree file */				 
	String message;						

	/** 
		Constructor. 
		
		Sets the population size.
		
		@param n the number of people
	*/	
	public FamilySet(int n)
	{
		people = n;
		message = null;
	}
	
	/** 
		Constructor. 
		
		Sets the population size and a string listing the errors associated
		with the pedigree file.
		
		@param n the number of people
		@param msg the error message
	*/	
	public FamilySet(int n, String msg)
	{
		people = n;
		message = msg;
	}
	
	/**
		Returns the number of people in the population.
		
		@return the number of people
	*/
	public int people()
	{
		return people;
	}
	
	/**
		Returns the errors associated with the pedigree file.
		
		@return the errors
	*/
	public String message()
	{
		return message;
	}

}