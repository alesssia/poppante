/*
 * 	 MyLinearRegression.java
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
 *   along with this code. If not, see <http://www.gnu.org/licenses/>.
 *
 *   For any bugs or problems found, please contact us at
 *   alessia.visconti@kcl.ac.uk
 */

package com.github.alesssia.algebrautils;

import org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression;


/**
	Performs a multivariate linear regression using Ordinary 
		Least Square.

	Specifically it estimates the beta and the residuals of multiple 
	linear regression:
		Y = beta*X + R
	where Y is a vector of length n and X is a matrix nxm, and R
    is a vector of residuals.

	@author      Alessia Visconti <alessia.visconti@kcl.ac.uk>
	@version     1.0
*/
public class MyLinearRegression
{
	/** The object used to solve the system */
	private OLSMultipleLinearRegression regressor; 
	
	/** The beta of multiple linear regression */
	private double[] beta;
	
	/**
		Constructor.
		
		Initialises the data structures as in OLSMultipleLinearRegression.newSampleData,
		that is assumes that rows are concatenated with y values first in each row. 
		
		For example, an input data array containing the sequence of values 
		(1, 2, 3, 4, 5, 6, 7, 8, 9) with nobs = 3 and nvars = 2 creates a regression 
		dataset with two independent variables, as below:

		   y   x[0]  x[1]
		   --------------
		   1     2     3
		   4     5     6
		   7     8     9
		
		
		@precondition the number of observations should be greater than  the number 
		of predictors (x's columns) 
		
		@param data data from a flat input array
		@param nobs the number of observations
		@param nvars the number of predictors
		@see org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression.newSampleData
	*/
	public MyLinearRegression(double[] data, int nobs, int nvars)
	{
		assert data.length == nobs * (nvars+1) : "MyLinearRegression: not well formed input matrix";
		assert nobs > nvars : "MyLinearRegression: not enough observations for this many predictors";
		
		regressor = new OLSMultipleLinearRegression();
		regressor.newSampleData(data, nobs, nvars);
		
		beta = null;
	}
	
	/**
		Return the beta of multiple linear regression.
		
		@precondition the beta should have been evaluated
		
		@return the beta of multiple linear regression
	*/
	public double[] beta()
	{
        if (beta == null)
			return null;
        
        return beta.clone();
	}
	
	/**
		Calculates the regression coefficients using Ordinary 
		Least Square.
	*/
	public void regress()
	{
		beta = regressor.estimateRegressionParameters();
	}
	
	
	
	/** 
		Estimates the residuals R, that is it evaluates:
			R = Y - betaX
		
		If not pre-velauated, it also calculates the regression 
		coefficients.
		
		@return the residuals of multiple linear regression
	*/
	public double[] residuals()
	{
		if (beta == null)
			regress();
			 
        return regressor.estimateResiduals();
	}
	
	
	
	
}
