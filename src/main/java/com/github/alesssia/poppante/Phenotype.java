/*
 * 	 Phenotype.java
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

import java.util.Arrays;
import java.util.Vector;
import java.util.Collections;
import org.apache.commons.lang.ArrayUtils;

/**
	Represent a phenotype trait.

	It can store any quantitative or dichotomic trait.

	@author      Alessia Visconti <alessia.visconti@kcl.ac.uk>
	@version     1.0     
*/

class Phenotype
{
	/** Phenotype name */
	private final String name;  		
	/** List of people (position) having a 
		missing value for this marker */
	private Vector<Long> missing;	
	
	/**
		Constructor. 
		
		@param n the phenotype name
	*/
	public Phenotype(String n)
	{
		name = n;
		missing = new Vector();
	}
	
	/**
		Returns the phenotype name.
		
		@return the phenotype name
	*/
	public String name()
	{
		return name;
	}
	
	
	/**
		Adds a person (position) to the vector of people
		having a missing value for this phenotype
		
		@param position the person position
		@see com.github.alesssia.poppante.MyFileReader
	*/
	
	public void addMissing(long position)
	{
		missing.add(position);
	}
	
	/**
		Returns the vector of people (positions)
		having a missing value for this phenotype
		
		@return positions of people having a missing value
		@see com.github.alesssia.poppante.DataManager
	*/
	public long[] missing()
	{
		//FIXME: bottlneck
		Collections.sort(missing);
		Object [] tmp = missing.toArray();
		return ArrayUtils.toPrimitive(Arrays.copyOf(tmp, tmp.length, Long[].class));
	}
	
	/** 
		Deletes the missing list
	*/
	public void resetMissing()
	{
		missing = null;
	}
}
