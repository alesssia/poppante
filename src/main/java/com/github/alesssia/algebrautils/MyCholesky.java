/*
 * 	 MyCholesky.java
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

import Jama.*;

/** 
	Performs the Cholesky decomposition.

	@author      Alessia Visconti <alessia.visconti@kcl.ac.uk>
	@version     1.0
*/

public class MyCholesky extends Decomposer
{
	/**  The result of the decomposition */
	CholeskyDecomposition choleskyDecomposition;  
	
	/**
		Constructor. 
		
		Initialises an empty object.	
	*/	
    public MyCholesky()
    {
		super();
		choleskyDecomposition = null;
    }
    
	/**
		{@inheritDoc}
		
		It can be applied to a square matrix only.
		
		@precondition the matrix must be square
	*/
    @Override
    protected void decompose(double[][] a) 
    {
		//Cholesky requires a square matrix
        assert a.length == a[0].length : "MyCholesky: the Cholesky can't be applied: the matrix is not square";
		
		choleskyDecomposition = new CholeskyDecomposition(new Matrix(a));
    }
	
	/**
		Gets the Cholesky factor L of the decomposed matrix A.
	
		@precondition a matrix A should have been decomposed
		
		@return L, the Cholesky factor
	*/
	public double[][] getL()
	{
		assert choleskyDecomposition != null : "MyCholesky: Cholesky Decomposition not initialised.";
		
		return choleskyDecomposition.getL().getArray();
	}
	
	/**
		Converts the solution.
		
		It recode the solution of the system from Jama.Matrix
		to an array of double.
		
		@param X the matrix
		@return an array
	*/
	private double[] convert(Matrix X)
	{
		double[][] Xv = X.getArray();
		double[] x = new double[Xv.length];
		for (int i=0; i<Xv.length; i++)
			x[i] = Xv[i][0];
	
		return x;	
	}

	/**
		{@inheritDoc}
	*/
    @Override
    public double[] solve(double[][] a, double[] b) 
    {
		assert a.length == b.length : "MyCholesky: system dimensions must agree";
		 
		decompose(a);
		Matrix X = choleskyDecomposition.solve(new Matrix(b, b.length));
		return convert(X);
    }

	/**
		{@inheritDoc}
	*/	
    @Override
    public double lnDet() 
    {
		assert choleskyDecomposition != null : "MyCholesky: Cholesky Decomposition not initialised.";
		
		//I do not need to perform the square root explicitly since 
		//it is already considered in the decomposition.
		return Math.log( (choleskyDecomposition.getL()).det()   );
    }
	
}
