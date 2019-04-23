/*
 * 	 FDR.java
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

/**
	Implements False Discovery Rate (FDR) control

	@author      Alessia Visconti <alessia.visconti@kcl.ac.uk>
	@version     1.0               
 */
public abstract class FDR
{
	/**  The unsorted p-values */
    protected Double[] pvalues;
	/** The adjusted p-values */
	protected double[] adjPvalues;
	
	/**
		Constructor.
		
		It creates an empty object
	*/	
	public FDR()
	{
		pvalues = null;
		adjPvalues = null;
	}
	             	
	/**
		Constructor.
		
		Initialises the FDR engine.
		
		@param pv the p-values vector
	*/	
    public FDR(Double[] pv)  
	{  
		pvalues = pv;
		adjPvalues = new double[pvalues.length];
	}
	
	
	/**
		Evaluates the adjusted p-values.
	*/	
	abstract public void adjust();
	
	/**
		Returns the adjusted p-value in the i-th position
	    
		@param i the position
	    @return the adjusted p-value
	*/	
	public double adjPvalue(int i)
	{
		return adjPvalues[i];
	}
}
	



