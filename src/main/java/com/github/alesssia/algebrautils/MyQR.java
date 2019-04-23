/*
 * 	 MyQR.java
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

/** 
	Performs the QR decomposition.

	@author      Alessia Visconti <alessia.visconti@kcl.ac.uk>
	@version     1.0     
*/

public class MyQR extends Decomposer
{
	/**  The result of the decomposition */
	private QR qrDecomposition; 
	
	/**
		{@inheritDoc}
	*/	
	public MyQR()
	{
		super();
		qrDecomposition = null;
	}
	
	
	/**
		{@inheritDoc}
	*/	
    @Override
    protected void decompose(double[][] a)
    {
        size = a.length;
		qrDecomposition = QR.factorize( new DenseMatrix(a) );
    }
	
	/**
		{@inheritDoc}
	*/	
    @Override
	public double[] solve(double[][] a, double[] b)
	{
		assert a.length == b.length : "MyQR: system dimensions must agree";
		
		//the matrix decomposition is evaluated only 
		//for the lnDet calculation. In fact MJT does not
		//have a specific solver codified for the QR, but
		//it uses the general Matrix solver. This general
		//solver decides whether to use LU or QR at runtime
		//(and since our matrix is squared LU is used)
		//However the usage of QR may change the determinant
		//the minimization step and thus the convergence of
		//the system.
        decompose(a); 
		
		//Solving the system using the general matrix solver.
		DenseMatrix A = new DenseMatrix(a);
        Matrix B = new DenseMatrix(new DenseVector(b));
        Matrix X = new DenseMatrix(new DenseVector(new double[b.length]));
		X = A.solve(B, X);
		
		// return new double[dim];
		return ((DenseMatrix)X).getData();
	}
	
	/**
		{@inheritDoc}
	*/			
    @Override
    public double lnDet() 
    {
		assert qrDecomposition != null : "MyQR: QR Decomposition not initialised.";
		assert size != -1 : "MyQR: QR dimension not set.";
		
		//To solve this, I will use the properties of logarithmic functions. In fact,
		//what I need to evaluate is the natural logarithm of the square root of 
		//determinant of the decomposition, that is ln(sqrt(det(A))), or ln(det(A))/2,
		//where det(A) = det(R)*det(Q), Since Q is unitary, |\det(Q)|=1. Thus,	
		//|\det(A)|=|\det(R)|=|\prod(R_{ii})|. Then I can transform
		//ln(det(A))/2=log(\prod(R_{ii}))/2=\Sigma(log(R_{ii}))/2.
		UpperTriangDenseMatrix R = qrDecomposition.getR();
		double result = 0.0;
		for (int i=0; i<size; i++)
			//FIXME: I am forcing an ABS here
			result += Math.log(Math.abs(R.get(i,i)));

        //qrDecomposition.getQ() is a bottleneck and it is only needed to
        //compute the sign. I am thus commeting it FIXME
		// det *= Math.sqrt(Math.abs(qrDecomposition.getQ().det())); //I still need to consider the unity part
			
		return result/2;
    }
}
	



