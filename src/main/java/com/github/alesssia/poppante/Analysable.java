/*
 *   Analysable.java
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

/**
	Represents the list of analisable families and
	family members (given a pair <phenotype, methylation site>).

	It is represented by two objects:
	- a list of analysable families, with their property
	- a vector of booleans representing which people are
	  analysable. The order of this vector follows the order
	  of the positionTable is the DataManager object, that
	  codifies the position of every other data structure.


	@author      Alessia Visconti <alessia.visconti@kcl.ac.uk>
	@version     1.0               
	@see		 com.github.alesssia.poppante.Test
	@see		 com.github.alesssia.poppante.VC
	@see		 com.github.alesssia.poppante.DataManager
*/

public class Analysable
{
	/** List of families that can be analysed, 
		along with the number of analysable 
		family members, and the start end position 
		in the datamanager's data structures */
	Vector<AnalysableFamily> analysableFamiles;
	/** Whether an individual is analysable 
		It follows the order of any other object 
		stored in datamanager (that is the order
		of the positionTable) */
	private	boolean[] isAnalysable;
	/** Number of analysable individuals */
	private int numAnalysable;
	
	/** 
		Constructor. 
		
		Initialises the data structures
	*/
	public Analysable()
	{
		analysableFamiles = new Vector<>();
		isAnalysable = null;
		numAnalysable = 0;
	}
		
	/**
		Adds a new analysable family
		
		@param family the family to add
	*/
	public void addFamily(AnalysableFamily family)
	{
		analysableFamiles.add(family);
		numAnalysable += family.numAnalysable();
	}
	
	/**
		Sets the which individuals are analysable
		
		@param a the array of position that can be analysed
	*/
	public void setIsAnalysable(boolean[] a)
	{
		isAnalysable = a.clone();
	}

	/**
		Returns the number of analysable families
		
		@return the number of analysable families
	*/
	public int numFamilies()
	{
		return analysableFamiles.size();
	}
	
	/**
		Returns the number of analysable individuals
		
		@return the number of analysable individuals
	*/
	public int numAnalysable()
	{
		return numAnalysable;
	}
	
	/**
		Returns whether the person in the i-th position is analysable
		
		@param i the position
		@return whether the person is analysable
	*/
	public boolean isAnalysable(int i)
	{
		return isAnalysable[i];
	}
	
	/**
		Returns the information of the i-th analysable family
		
		@param i the position
		@return the i-th analysable family
	*/
	public AnalysableFamily analysableFamiles(int i)
	{
		return analysableFamiles.get(i);
	}
	
}


/**
	Helper class for Analysable.

	Represents an analysable family as its ID and
	describes some of its properties, that is:

	- the number of analysable members
	- the position where this family starts in the
	  position Table, information that is used by
	  the VC to indicise the members' data in the
	  DataManager's data structures.

	@author      Alessia Visconti <alessia.visconti@kcl.ac.uk>
	@version     1.0               
	@see		 com.github.alesssia.poppante.DataManager
*/
	
class AnalysableFamily
{

	/** Family ID */
	private final String famID;
	/** Number of analysable member */
	private final int numAnalysable;
	/** Position in the datamanager's data 
		structure where this family starts */
	private final int startPos;
	
	/** 
		Constructor. 
		
		Defines a new analysable family
		
		@param id the family's ID
		@param n the number of analysable member
		@param start position in the datamanager's data structure where this family starts
	*/
	public AnalysableFamily(String id, int n, int start)
	{
		famID = id;
		numAnalysable = n;
		startPos = start;
	}
	
	/**
		Returns the family's ID
		
		@return the family's ID
	*/
	public String famID()
	{
		return famID;
	}
	
	/**
		Returns the number of analysable family's member
		
		@return the  the number of analysable family's member
	*/
	public int numAnalysable()
	{
		return numAnalysable;
	}
	
	/**
		Returns the position in the datamanager's data 
		structure where this family starts
		
		@return the start position
	*/
	public int startPos()
	{
		return startPos;
	}
	
}
