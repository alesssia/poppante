/*
 * 	 Decomposer.java
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

package com.github.alesssia.algebrautils;

/**
	Performs the matrix decomposition.

	@author      Alessia Visconti <alessia.visconti@kcl.ac.uk>
	@version     1.0               
 */
public abstract class Decomposer
{
	/** Number close the smallest representable number */	
	public static final double FPMIN   = 1.0e-320;   
 
	/**  The size of the original (square) matrix. */
    protected int size;             
	
	/**
		Constructor.
		
		Initialises an empty object.
	*/	
    public Decomposer()  
	{  
		size = -1;
	}
	
	
	/**
		Decomposes the matrix.
		
		@param a the matrix to decompose
	*/	
	abstract protected void decompose(double[][] a);
	
	/**
		Solves a system of linear equations of the form 
		x=Ab.
		
		It verifies that the dimensions of systems agree. 
		
		@precondition the dimensions of the systems must agree
		
		@param a the matrix A
		@param b the vector b
		@return the solution of the system
	*/
	abstract public double[] solve(double[][] a, double[] b);
	
	/**
		Evaluates the natural logarithm of the square root of 
		determinant of the decomposition.
		
		It corresponds to ln(det(A))/2 or ln(sqrt(det(A)))
		
		@return the logarithm of the square root of the determinant
	*/	
    abstract public double lnDet();
	
}
	



