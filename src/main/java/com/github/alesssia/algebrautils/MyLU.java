/*
 * 	 MyLU.java
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

import no.uib.cipr.matrix.*;
import java.util.Arrays;

/** 
	Performs the LU decomposition.

	@author      Alessia Visconti <alessia.visconti@kcl.ac.uk>
	@version     1.0     
*/

public class MyLU extends Decomposer
{
	
	/**  The result of the decomposition */
	private DenseLU luDecomposition;  
	
	/**
		{@inheritDoc}
	*/	
    public MyLU()
    {
		super();
		luDecomposition = null;
    }
    
	/**
		{@inheritDoc}
		
		It can be applied to a square matrix only.
		
		@precondition the matrix must be square
	*/
    @Override
    protected void decompose(double[][] a) 
    {
		//LU requires a square matrix
        assert a.length == a[0].length : "MyLU: the LU decomposition can't be applied: the matrix is not square";
		
		size = a.length;
		luDecomposition = DenseLU.factorize( new DenseMatrix(a) );
    }

	/**
		{@inheritDoc}
	*/
    @Override
    public double[] solve(double[][] a, double[] b) 
    {
		assert a.length == b.length : "MyLU: system dimensions must agree";
		 
		decompose(a);
		return luDecomposition.solve(new DenseMatrix(new DenseVector(b))).getData();
    }
	

	/**
		{@inheritDoc}
	*/	
    @Override
    public double lnDet() 
    {
		assert luDecomposition != null : "MyLU: LU decomposition not initialised.";
		assert size != -1 : "MyLU: LU dimension not set.";
		
		//getData returns an array describing:
		// - in the lower triangular matrix the L part of the decomposition 
		//   (since the diagonal of L is set to 1 this information is not 
		//   represented). 
		// - in the upper triangular matrix the U part of the decomposition 
		double[] decomposition = luDecomposition.getL().getData();
		
		//To solve this, I will use the properties of logarithmic functions. In fact,
		//what I need to evaluate is the natural logarithm of the square root of 
		//determinant of the decomposition, that is ln(sqrt(det(A))), or ln(det(A))/2,
		//where det(A) = (-1)^S\prod(L_{ii})\prod(U_{ii}). 
		//First I can observe that L_{ii} == 1 for each i, and thus det(A) can be 
		//evaluated as: det(A) = (-1)^S\prod(U_{ii}). Then I can transform
		//ln(det(A))/2=log(\prod(U_{ii}))/2=\Sigma(log(U_{ii}))/2.
		// FIXME in this implementation (-1)^S has been ignored.
		double result = 0.0;
		//Of the vector representation of the decomposition, I will need only the diagonal elements
		//that are described by the following recursive procedure
		// - n_0 = 0
		// - n_n = n_(n-1) + dim + 1
		for (int i=0; i<decomposition.length; i+=size+1 )
			result += Math.log(decomposition[i]);
		
		return result/2;
	}
}
