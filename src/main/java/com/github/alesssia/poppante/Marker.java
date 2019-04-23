/*
 * 	 Marker.java
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
 *   mario.falchi@kcl.ac.uk  
 *   alessia.visconti@kcl.ac.uk
 */

package com.github.alesssia.poppante;

import java.util.*;
import org.apache.commons.lang.ArrayUtils;
import java.util.Collections;

/**
	Represent a marker along with its coordinate.

	A marker is a methylation site, but it can
	also be a SNP.

	@author      Alessia Visconti <alessia.visconti@kcl.ac.uk>
	@author      Mario Falchi     <mario.falchi@kcl.ac.uk>
	@author      Hariklia Eleftherohorinou <hariklia.eleftherohorinou06@imperial.ac.uk>
	@version     1.0     
*/

class Marker
{
	/** Marker name */
	private final String name; 
	/** Chromosome */ 		
	private final String chromosome; 	
	/** Position (bp) */
	private final long position;	
	/** List of people (position) having a 
		missing value for this marker */
	private Vector<Long> missing;	


	/**
		Constructor. 
	
		@param name the marker name
		@param chromosome the chromosome 
		@param position the position where
	*/
	public Marker (String name, String chromosome, long position)
	{
		this.chromosome = chromosome;
		this.name = name;
		this.position = position;
		
		missing = new Vector();
	}
	
	
	/**
		Constructor. 
	
		@param name the marker name
	*/
	public Marker (String name)
	{
		this.chromosome = null;
		this.name = name;
		this.position = -1;
		
		missing = new Vector();
	}

	/**
		Returns the marker name.
	
		@return the marker name
	*/
	public String name()
	{
		return name;
	}

	/**
		Returns the marker chromosome.
	
		@return the chromosome
	*/
	public String chromosome()
	{
		return chromosome;
	}

	/**
		Returns the marker position
	
		@return the position
	*/
	public long position()
	{
		return position;
	}


	/**
		Extracts from a list of markers those that are within
		a given window from the current position. 
		
		The methods do not return a list of Markers but the
		a list of Integer representing the positions of the
		selected Markers in the supplied list.
		It is used in the region-based test evaluations.
	
		@precondition the window for the region-based test should have been defined
	
		@param list the list of Markers
		@return the list of indices of Markers that are within the window
		@throws Exception if a precondition is not satisfied
		@see com.github.alesssia.poppante.DataManager
		@see com.github.alesssia.poppante.PCA
	*/
	public Vector<Integer> getMarkersInWindow(Vector<Marker> list) throws Exception
	{
		assert Constants.region != Utilities.INVALID_I : "Internal error: region-base testing called when single predictor is required.";

		Vector<Integer> indices = new Vector<Integer>();
	
		for(int i=0; i<list.size(); i++)
		{
			Marker that = list.get(i);
		
			//check the chromosome (if they are on different chromosomes continues)
			if (!this.chromosome.equals(that.chromosome()))
				continue;
		
			boolean inRegion = that.position() >= this.position && that.position() <= (this.position + Constants.region);
			if (inRegion) 
				indices.add(i);
		}

		return indices;
	}
	
	
	/**
		Adds a person (position) to the vector of people
		having a missing value for this marker
		
		@param position the person position
		@see com.github.alesssia.poppante.MyFileReader
	*/
	
	public void addMissing(long position)
	{
		missing.add(position);
	}
	
	/**
		Returns the vector of people (positions)
		having a missing value for this marker
		
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
