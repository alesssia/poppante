/*
* 	 MyPCA.java
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

import java.util.Arrays;
import Jama.*;
 
/** 
	Performs the Principal Component Analysis using the 
	SVD decomposition and evaluate the proportion of overall 
	variability each principal components accounts for.

	The input data matrix n-by-m, having n variables (in rows) and
	m observations (in columns) is expressed as a data matrix m-by-n, 
	having m observations (in row) and n principal components (in cols).

	The number of observations m must be greater or equals to the 
	number of variables n.

	The proportion of overall variability is expressed in [0,1], as
	any proportion.

	It uses	the "Jama" library for the SVD evaluation.

	@author      Alessia Visconti <alessia.visconti@kcl.ac.uk>
	@version     1.0     
*/

public class MyPCA 
{
	/** Input data matrix  */
	private double[][] input; 
	/** The value used to represent missing  */
	private double missingValue; 
	/** Whether the input matrix includes missing values  */
	private boolean hasMissing; 
	/** Transformed marix used for PCA  */
	private double[][] data; 
	/** Eigenvector matrix whose column are sorted by eigenvalues */
	private double[][] pca; 
	/** Proportion of variance (in [0,1]) */
	private double[] proportionOfVariance;
	
	/** Small number used to compare double */	
	public static final double DELTA = 1.0e-5;  
 
	/**
		Constructor.
		
		Initialises the data structure.
		
		The row of the data matrix describes the variable,
		while the columns describes the observations.
		It allows matrix with missing values that will
		be transformed during the PCA
		
		@precondition the number of observations m must be greater or equals to the 
		number of variables n.
		
		@param a the input matrix
		@param missing the value used to represent missing values
	*/
	public MyPCA(double[][] a, double missing) 
	{
		assert a[0].length >= a.length : "MyPCA : the number of observations is smaller than the number of variables to estimate";
		
		input = a;
		
		//The input matrix is tranpsosed (it ensures
		//also that input matrix is not modified)
		Matrix m = new Matrix(input);
		data = m.transpose().getArray();
		
		missingValue = missing;
		hasMissing = true;
		
        pca= null;
		proportionOfVariance = null;
	}
	
	/**
		Constructor.
		
		Initialises the data structure.
		
		The row of the data matrix describes the variable,
		while the columns describes the observations.
		
		@precondition the number of observations m must be greater or equals to the 
		number of variables n.
		
		@param a the input matrix
		@param missing the value used to represent missing values
	*/
	public MyPCA(double[][] a) 
	{
		assert a[0].length >= a.length : "MyPCA : the number of observations is smaller than the number of variables to estimate";
		
		input = a;
		
		//The input matrix is tranpsosed (it ensures
		//also that input matrix is not modified)
		Matrix m = new Matrix(input);
		data = m.transpose().getArray();
		
		missingValue = Double.MAX_VALUE;
		hasMissing = false;
		
        pca= null;
		proportionOfVariance = null;
	}
	
	/**
		Returns the modified data used internally
		by the PCA.
		
		If called before evaluatePCA() data is equal to the transposed
		input matrix, otherwise it correspond to the transformed
		data matrix, having the missing values set to zero, and 
		being transpose, centered and scaled.
		
		@return the data matrix used for the PCA
	*/
	public double[][] data()
	{
		return data;
	}
	
	/**
		Returns the proportion of overall variability explained
		by each principal components. 
		
		The position in the represents the principal component.
		
		@return the input data
	*/
	public double[] proportionOfVariance()
	{
		if (proportionOfVariance == null)
			return null;
		
		return proportionOfVariance.clone();
	}
	
	
	/**
		Performs the Principal Component Analysis.
		
		Before the SVD decomposition the input data matrix is 
		transposed, the missing data transformed to zero, and
		the new matrix is scaled and centred 
	*/
	public void evaluatePCA() 
	{
		//The input matrix is scaled and centered
		normalise();
		
		//reset of missing values
		if (hasMissing)
			resetMissing();
		
		//The input matrix is decomposed an the left singular vectors
		//are assigne to the PCA
		SingularValueDecomposition svd = new SingularValueDecomposition(new Matrix(data));
		pca = svd.getU().getArray();
	}
	
	/**
		Returns all the principal components
		
		@return the PCA results
	*/
	public double[][] allPCs()
	{
		return pca;
	}
	
	/**
		Returns the i-th principal component.
		Principal components starts from 1.
		
		@precondition the i-th principal component should exist
		
		@param i which principal component
		@return the principal component
	*/
	public double[] getPC(int i)
	{
		assert i > 0 && i <= pca[0].length : "MyPCA: the number of principal components does not exist";
		
		//PCs matrix colunts from 0, so the
		//ith PC is the i-1-th col
		i--;
		
		//PCs are stored in columns
		double[] pc = new double[pca.length];
		for (int j=0; j<pca.length; j++)
       		pc[j] = pca[j][i];
		
		return pc;
	}
	
	/**
		Returns the first n-th principal components,
		that is the principal components from 1 to n.
		
		Rows represents observations, columns principal
		 components
		
		@param n the last principal component to extract
		@return the selected principal components
	*/
	public double[][] getPCs(int n)
	{
		assert n > 0 && n <= pca[0].length : "MyPCA: the number of principal components does not exist";
		
		//all PCs requested
		if (n == pca[0].length)
			return pca;
				
		double[][] pcs = new double[pca.length][n];
		for (int i=0; i<n; i++)
			for (int j=0; j<pca.length; j++)
				pcs[j][i] = pca[j][i];
		
		return pcs;
	}
	
	/** 
		Set the missing values to zero.
		
		Since missing values are not allowed, this function 
		transforms all the value codified as missing (by using 
		given value) to zero.
	*/
	private void resetMissing()
	{
		for (int i=0; i<data.length; i++)
			for (int j=0; j<data[0].length; j++)
				if (Math.abs(data[i][j] - missingValue) < DELTA)
					 data[i][j] = 0;
	}
	
 
	
 
 	/**
		Centers and scales the data.
	*/
	private void normalise() 
	{
		double[] colmean = colmean();
        double[] colsd = colsd(colmean);         	
		
		for (int i = 0; i < data.length; i++) 
			for (int j = 0; j < data[0].length; j++) 
			{
				//Does not modify missingValue
				if (Math.abs(data[i][j] - missingValue) < DELTA) 
                                {}
                                //Avoids division by zero
                                 else if (colsd[j] == 0)
					data[i][j] = 0.0;
				//scales and centers
				else
					data[i][j] = (data[i][j]-colmean[j])/colsd[j];
			}
	}
 
   
	/**
		Evaluates the column means.
		
		Values set to missing (if any) do not contribute 
		to the mean evaluation.
	
		@return the column mean
	*/
   	private double[] colmean()
	{
		//It should not consider missing values
		double [] size = new double[data[0].length];
		double [] mean = new double[data[0].length];
		
		for (int i = 0; i < data.length; i++) 
			for (int j = 0; j < data[0].length; j++) 
				if(!hasMissing || Math.abs(data[i][j] - missingValue) > DELTA)
				{
					mean[j] += data[i][j];
					size[j]++;
				}
		
		for (int j = 0; j < mean.length; j++) 
			mean[j] /= size[j];
		
		return mean;
	}
	
	/**
		Evaluates the column standard deviation
		
		Values set to missing (if any) do not contribute 
		to the mean evaluation.
	
        @param colmean the column mean
		@return the column standard deviation
	*/
	private double[] colsd(double[] colmean)
	{
		//It should not consider missing values
		double [] size = new double[colmean.length];
		double [] sd = new double[colmean.length];
		
		for (int i = 0; i < data.length; i++) 
			for (int j = 0; j < data[0].length; j++) 
				if(!hasMissing || Math.abs(data[i][j] - missingValue) > DELTA)
				{
					sd[j] += (data[i][j]-colmean[j])*(data[i][j]-colmean[j]);
					size[j]++;
				}
					
		for (int j = 0; j < sd.length; j++) 
			sd[j] = Math.sqrt(sd[j]/(size[j]-1));

		return sd;			
	}
	
	/**
		Evaluates the (sample) variance of the vector v.
		
		@param v a vector
		@return the sample variance of the vector
	*/
   	private static double variances(double[] v)
	{
		double mean = 0.0;
		for (int i=0; i<v.length; i++)
			mean += v[i];
		mean /= v.length;
		
	    
	    double s = 0.0;
		for (int i=0; i<v.length; i++)
        	s += (mean-v[i])*(mean-v[i]);
        
		return  s/(v.length-1);
	}
	
	
	/**
		Evaluates the proproportion of overall variability 
		each principal component accounts for.
		
		@precondition the PCA has been evaluated
	*/
	public void evaluateProportionOfVariance()
	{
		assert pca != null : "MyPCA : PCA not evaluated yet";

		//Transpose the PCA to simplify further calculations
		Matrix m = new Matrix(pca);
		double[][] tmp = m.transpose().getArray();

		//gets the variances of the PCs
		proportionOfVariance  = new double[tmp.length];
		double sum = 0.0;
		for (int i=0; i<tmp.length; i++)
		{
			proportionOfVariance [i] = variances(tmp[i]);
			sum += proportionOfVariance[i];
		}
		
		//evaluates the proportions 
		for (int i=0; i<proportionOfVariance.length; i++)
			proportionOfVariance[i] /= sum;
	}
	
	/**
		Returns the number of PCs necessary to account for at
		least the given proportion of variance.
		
		@param minProportionOfVariance the minim variance one would account for
		@return the number of PCs necessary to account that proportion of variance
	*/
	public int howMany(double minProportionOfVariance)
	{
		int numPC = 0;
		double cumulativeVariance = 0.0;
		while (cumulativeVariance < minProportionOfVariance)
		{
			cumulativeVariance += proportionOfVariance[numPC];
			numPC++;
		}
		
		return numPC;
	}

}
 
