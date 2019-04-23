/*
 *   NullModel.java
 *
 *   PopPAnTE is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   PopPAnTE is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with PopPAnTE.  If not, see <http://www.gnu.org/licenses/>.
 *
 *   For any bugs or problems found, please contact us at
 *   mario.falchi@kcl.ac.uk  
 *   alessia.visconti@kcl.ac.uk
 */

package com.github.alesssia.poppante;

/**
	Represent a Null model as its log likelihood and its degree 
	of freedom.

	Null models are used to speed up the computation when multiple
	methylation sites and phenotypes show the same patter of missing
	values, that is: the same people have missing values for the 
	same methylation sites and phenotype values.

	@author      Alessia Visconti <alessia.visconti@kcl.ac.uk>
	@version     1.0               
	@see		 com.github.alesssia.poppante.VC
	@see		 com.github.alesssia.poppante.DataManager
*/

public class NullModel
{
	/** This is the unique code identifying the 
		missingness pattern */
	private String missingnessPattern;
	/**  Log logLikelihood of the model */
	private double logLikelihood;   
	/**  Degrees of freedom of the model */
	private int df; 
	/** Environmental and/or generic variance of the
		null model */
	private double[] variances;
	/** -likelihood of every equations
		belonging to the system */
	private double[] details;
		
	
	/**
		Constructor. 
		
		Creates an empty object
	*/
	public NullModel()
	{
		missingnessPattern = null;
		logLikelihood = Utilities.INVALID_D;
		df = Utilities.INVALID_I;
		variances = null;
		details = null;
	}
	
	/**
		Constructor. 
		
		Initialises the unique code identifying the
		missingness pattern.
		
		@param pattern the unique code 
	*/
	public NullModel(String pattern)
	{
		missingnessPattern = pattern;
		logLikelihood = Utilities.INVALID_D;
		df = Utilities.INVALID_I;
		variances = null;
		details = null;
	}
	
	/**
		Returns the unique code identifying the
		missingness pattern.
		
		@return the unique code
	*/
	public String missingnessPattern()
	{
		return missingnessPattern;
	}
	
	/**
		Returns the log Likelihood of the model
		
		@return the log Likelihood
	*/
	public double logLikelihood()
	{
		return logLikelihood;
	}
	
	/**
		Returns the degrees of freedom of the model
		
		@return the degrees of freedom
	*/
	public int df()
	{
		return df;
	}
	
	/**
		Returns the variance of the model.
		
		When association test is performed both generic and environmental
		variances are available. When the heritability test is performed
		only the latter is returned,
		
		@return the degrees of freedom
	*/
	public double[] variances()
	{
		return variances.clone();
	}
	
	/**
		Returns -likelihood of every equations belonging to the system
		
		@return the -logLikelihoods
	*/
	public double[] details()
	{
		return details.clone();
	}
	
	/**
		Returns 
		
		@return true if the model has already been evaluated, false otherwise 
	*/
	public boolean isSet()
	{
		return variances != null;
	}
	
	/**
		Sets the unique code identifying the
		missingness pattern.
		
		@param pattern the unique code
	*/
	public void setMissingnessPattern(String pattern)
	{
		missingnessPattern = pattern;
	}
	
	/**
		Sets the log Likelihood of the model
		
		@param value the log Likelihood
	*/
	public void setlogLikelihood(double value)
	{
		logLikelihood = value;
	}
	
	/**
		Sets the degrees of freedom of the model
		
		@param value the degrees of freedom
	*/
	public void setDf(int value)
	{
		df = value;
	}
	
	/**
		Sets the variance of the model.
		
		@param values the variance(s)
	*/
	public void setVariances(double[] values)
	{
		variances = values;
	}
	
	/**
		Sets the -likelihood of every equations belonging to the system
		
		@param logLikelihoods the -logLikelihoods
	*/
	public void setDetails(double[] logLikelihoods)
	{
		details = logLikelihoods.clone();
	}
	
}