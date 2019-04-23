/*
 * 	 BH.java
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

package com.github.alesssia.statutils;

import com.github.alesssia.arrayutils.*;
import java.util.Arrays;

/** 
	Implements the Benjamini–Hochberg procedure (BH step-up procedure) to
	control the false discovery rate. 

	References for the method can be found at:
	Benjamini, Yoav, and Yosef Hochberg. "Controlling the false discovery rate:
    a practical and powerful approach to multiple testing." Journal of the Royal 
    Statistical Society. Series B (Methodological) (1995): 289-300.

	@author      Alessia Visconti <alessia.visconti@kcl.ac.uk>
	@version     1.0     
*/

public class BH extends FDR
{
	/**  The sorted p-values */
    private RankVector pvaluesRank;
	
	
	/**
		Constructor.
		
		It creates an empty object
	*/	
	public BH()
	{
		super();
		pvaluesRank = null;
	}
	
	/**
		Constructor.
		
		It initialises the data structures.
		
		@param pv the pvalues vector
	*/	
	public BH(Double[] pv)
	{
		super(pv);
        pvaluesRank = null;
	}
	
	/**
		Implementes the BH procedures.
		
		p-values are ranked and then multiplied by the number of tests and divided by its assigned
		rank to give the adjusted p-values.
	*/	
    @Override
    public void adjust()
    {
		pvaluesRank = new RankVector(pvalues);
		pvaluesRank.rank();
		
		final int n = pvalues.length;
		double min = Double.MAX_VALUE; 
        for(int i=n-1; i>=0; i--)
		{
			//rank vector counts from zero
			double adj = (double)n/(i+1)*pvaluesRank.getAtRank(i);
			//I use the min to implemen the cumulative minimum
			if (adj < min)
				min = adj;
			else
				adj = min;
			
			adjPvalues[pvaluesRank.whereIsRank(i)] = adj;
		}	
    }
}	
	



