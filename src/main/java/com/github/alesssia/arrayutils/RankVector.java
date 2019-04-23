/*
 * 	 RankVector.java
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

package com.github.alesssia.arrayutils;

import java.util.*;
import com.github.alesssia.probutils.*;

/**
	Represents the rank of a vector of Double.

	@author      Alessia Visconti <alessia.visconti@kcl.ac.uk>
	@version     1.0     

	@see com.github.alesssia.poppante.ArrayComparator
*/

public class RankVector
{
	/** Vector */
    private final Double[] vector;	
	/** Rank Vector */
    private Integer[] rank;	
	/** Keeps the order in the Ranks */
	private List rankList;	

	/**
		Constructor. 
		
		Creates a rank vector from a given vector of Double.
		
		@param v the vector of Double
	*/
    public RankVector(Double[] v)
    {
		vector = v;
		rank = null;
    }

	/**
		Returns the Double value at the i-th position.
		
		@param i the position
		@return the Double value
	*/
    public Double get(int i)
    {
		return vector[i];
    }

	/**
		Returns the rank of the Double value at the i-th position
		
		@param i the position
		@return the rank of Double value
	*/
    public Integer getRank(int i)
    {
		return rank[i];
    }
	
	/**
		Returns the Double value at the i-th position
		
		@param i the position
		@return the Double value
	*/
    public Double getAtRank(int i)
    {
		int index =  rankList.indexOf(i);
		return vector[index];
    }
	
	/**
		Returns the position of the i-th ranked object
		
		@param i the position
		@return the Double value
	*/
    public int whereIsRank(int i)
    {
		return rankList.indexOf(i);
    }


	/**
		Evaluates the ranking of the Double vector.
	*/
    public void rank()
    {
		ArrayComparator<Double> comparator = new ArrayComparator<>(vector);
		rank = comparator.createIndexArray();
		
		//the ranklist is used to know who is at which rank (see BH)
		rankList = Arrays.asList(rank);
    }
	
	/**
		Evaluates the quantile normalisation (normal inverse transformation)
		for the rank vector. 
		
		@throws RuntimeException if the p-value of one of the transformed value is outside range 
		@return the transformed values
	*/
	public double[] transform() throws RuntimeException
	{
		int n = vector.length;
	    double scale = 1.0 / n;
	    double[] zs = new double[n];
    		
	    int j = -1;
	    for(int i=0; i < n; i++)
	    {
			for(j = i; j < n-1; j++)
			    if (!Objects.equals(vector[rank[i]], vector[rank[j]]))
					break;
					
				if (!Objects.equals(vector[rank[i]], vector[rank[j]]))
					j--;
					
			//ninv may throw an Exception if here the pvalue is outside range in Normal Inverse Transformation	
			double z = MyNormalDistribution.ninv(((i + j) * 0.5 + 0.5) * scale);
			for (int k = i; k <= j; k++)
			    zs[k] = z;

			i=j;
	    }
		
		//sots the transformed values according to the input
		double[] zsSort = new double[zs.length];
		for(int i=0; i<n; i++)
			zsSort[i] = zs[rank[i]];
			
		return zsSort;
	}

}

/**
	Sorts arrays and creates Index Arrays.

	@author      Alessia Visconti <alessia.visconti@kcl.ac.uk>
	@version     1.0     
        
    @param <T extends Comparable> the type of the value being boxed

	@see com.github.alesssia.arrayutils.RankVector
*/

class ArrayComparator<T extends Comparable<? super T>>
{
	/** Array to compare */	
    private final T[] array;		
	/** Sorted array */
    private SortedMap<T, Vector<Integer>> sortedArray; 

	/**
		Constructor. 
		
		Initialises the data structures.
		
		@param array the array
	*/
    public ArrayComparator(T[] array)
    {
        this.array = array;
		sort();
    }
	
	/**
		Sorts the array.
	*/
	private void sort()
	{
        sortedArray = new TreeMap<>();
        for(int i = 0 ; i < array.length ; i ++)
		{
		    Vector<Integer> tmp = sortedArray.get(array[i]);
		    if (tmp == null)
				tmp = new Vector<>();

		    tmp.add(i);
		    sortedArray.put(array[i], tmp);
		}
	}
	
	/**
		Created an index array from a sorted array.
		
		@precondition the index array must be sorted
		
		@return the Index Array
	*/	
    public Integer[] createIndexArray()
    {
		assert sortedArray != null : "Internal error: index array needs a sorted array.";
		
        Integer[] indexArray = new Integer[array.length];
        int i = 0;

        for(T key : sortedArray.keySet())
		{
		    Vector<Integer> tmp = sortedArray.get(key);
            for (Integer t : tmp) 
			{
                int index = t;
                indexArray[index] = i;
                i++;
            }
        }
	
        return indexArray;
    }

}




