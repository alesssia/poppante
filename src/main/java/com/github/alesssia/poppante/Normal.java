/*
 * 	 Normal.java
 *
 *   PopPAnTe is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   PopPAnTe is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with PopPAnTe.  If not, see <http://www.gnu.org/licenses/>.
 *
 *   For any bugs or problems found, please contact us at
 *   mario.falchi@kcl.ac.uk  
 *	 hariklia.eleftherohorinou06@imperial.ac.uk
 *   alessia.visconti@kcl.ac.uk
 */

package com.github.alesssia.poppante;

import com.github.alesssia.algebrautils.*;

/**
	Represents a equation. 
	
	It can be used to describe both full and null model.

	References for the implementation can be found at:
		Press, William H., et al. Numerical recipes in C. 
		Vol. 2. Cambridge: Cambridge university press, 1996.

	@author      Alessia Visconti <alessia.visconti@kcl.ac.uk>
	@author      Mario Falchi     <mario.falchi@kcl.ac.uk>
	@author      Hariklia Eleftherohorinou <hariklia.eleftherohorinou06@imperial.ac.uk>
	@version     1.0  
	@see com.github.alesssia.poppante.NormalSet
	@see com.github.alesssia.poppante.VCs
*/

class Normal 
{
	/** Beta coefficients */
	private double[] means;
	/** Variances */
	private double[] variances;
	/** Variant components (random effect) arrays.
		
		It contains two components: the intercepts
		and the kinship information */
	private double[][][] varComponents;
	/** Fixed effects.
		
		They include the tested value (full model only)
		and the covariates (always) */	
	private double[][] linearModel; //tested values are the methylation sites for association only
	/** Outcomes to be predicted.*/
	private double[] scores;	//Methylation sites for heritability, phenotypes for association
	/** The object used to decompose and solve the system */
	private final Decomposer decomposition;
	/** The likelihood of the system */
	private double likelihood;
	/** A constant to add to the log likelihood.
		
		It is equal to Math.log(2 * Math.PI) * -0.5 * dim
		where dim is the number of observations. */
	private double constant;
	/** Number of observations */
	private int dim;
	/** Number of fixed effects in the system.*/
    private int linearMD;
	/** Number of random effects in the system.*/
    private final int vcCount;
	
	/** Whether a constant has been added to the likelihood */
	private boolean includeLikelihoodConstant;
	/** Whether the model has been evaluated 
		
		It used to decide if add the constant when the enableConstant 
		and disableConstant are called.*/
	private boolean init;
	
	/**
		Constructor.
		
		It initialises the data structure.
		
		@param lmd number of fixed effects
		@param vc number of random effects
	*/
	public Normal(int vc, int lmd)
	{
		//decomposition is chosen ad hoc
		//Pedigree-based kinship is always positive definite, so I can use
		//Cholesky decomposition -- that can be used also if the external
		//kinship has been mofified with a bending' procedure to modify
		//the eigenvalues of non-positive definite matrices (that is the
		//default behaviour)
		if (Constants.kinship == null || Constants.decomposition == null)
			decomposition = new MyCholesky();
		else if (Constants.decomposition.equals("QR"))
			decomposition = new MyQR();
		else
			decomposition = new MyLU();
				
		
		
		vcCount = vc;
		linearMD = lmd;
        dim = -1;
	}
	
	/**
		Updates the number of fixed effects.
		
		@param n the new number of fixed effects
	*/
	public void resizeLM(int n)
	{
		linearMD = n;
	}

	/**
		Returns the system likelihood.
		
		@return the likelihood
	*/	
	public double likelihood()
	{
		return likelihood;
	}

	/**
		Returns the vector of variances.
		
		@return the vector of variances
	*/	
	public double[] variances()
	{
		return variances;
	}
	
	/**
		Returns the vector of beta coefficients.
		
		@return the vector of beta coefficients
	*/	
	public double[] means()
	{
		return means;
	}

	/**
		Returns the outcomes.
		
		@return the outcomes
	*/	
	public double[] scores() 
	{
		return scores;
	}
	
	/**
		Returns the i-th outcome.
		
		@param i the position
		@return the outcomes
	*/	
	public double scores(int i)
	{
		return scores[i];
	}
	
	/**
		Sets the i-th outcome to value
		
		@param i the position
		@param value the value
	*/
	public void setScores(int i, double value)
	{
		scores[i] = value;
	}
	
	/**
		Returns the i-th predictor.
		
		@param i the position
		@return the predictor
	*/	
	public double[] linearModel(int i)
	{
		return linearModel[i];
	}
	
	/**
		Sets the j-th value of the i-th predictor to value
		
		@param i the position of the predictor
		@param j the position of the value
		@param value the value
	*/
	public void setLinearModel(int i, int j, double value)
	{	
		linearModel[i][j] = value;
	}
	
	/**
		Returns the i-th random effect.
		
		@param i the position
		@return the random effect
	*/	
	public double[][] varComponents(int i)
	{
		return varComponents[i];
	}
	
	/**
		Returns the [j-th, k-th] value of the i-th random effect to value
		
		@param i the position of the random effect
		@param j the position in the random effect
		@param k the position in the random effect
		@return the value
	*/
	public double varComponents(int i, int j, int k)
	{
		return this.varComponents[i][j][k];
	}
	
	/**
		Sets the [j-th, k-th] value of the i-th effect to value
		
		@param i the position of the random effect
		@param j the position in the random effect
		@param k the position in the random effect
		@param value the value
	*/
	public void setVarComponents(int i, int j, int k, double value)
	{
		//It makes a full matrix out of a triangular one
		varComponents[i][j][k] = value;
		varComponents[i][k][j] = value;
	}
	
	/**
		Returns the number of observations
		
		@return the number of observations 
	*/
	public int dim()
	{
		return dim;
	}
	
	/**
		Returns the j-th predictor.
		
		It is used by the permutation tests.
		
		@param j the position
		@return the predictor
		@see com.github.alesssia.poppante.VC
	*/
	public double[] predictor(int j)
	{
		double[] predictor = new double[dim];
		for(int i=0; i<dim; i++)
			predictor[i] = linearModel[i][j];
	    
		return predictor;
	}
	
	/**
		Sets the j-th predictor.
		
		It is used by the permutation tests.
		
		@param predictor the predictor
		@param j the position
		@see com.github.alesssia.poppante.VC
	*/
	public void setPredictor(double[] predictor, int j)
	{
		for(int i=0; i<dim; i++)
			linearModel[i][j] = predictor[i];
	}
	
	
	/**
		Initialises the data structure
	
		@param dim the number of observations
		@see com.github.alesssia.poppante.VC
	*/
	public void prepare(int dim) 
	{		
		includeLikelihoodConstant = false;
		init = false; 
		
		this.dim = dim;
		
		// initialises dara structure
		constant = Math.log(2 * Math.PI) * -0.5 * dim;
		
		scores = new double[dim];				 		
		means = new double[linearMD]; 		
		variances = new double[vcCount];	
		linearModel = Utilities.set(dim, linearMD, 1.0);
		
		varComponents = new double[vcCount][dim][dim];
		varComponents[0] = new double[dim][dim];	
		for(int i=0; i<dim; i++)
			varComponents[0][i][i] = 1.0;
	}
	
	/**
		Sets the beta coefficients and the variances
		
		@param mu the beta coefficients
		@param sigma the variances
	*/
	public void setParameters(double[] mu, double[] sigma)
	{
		means = mu;
		variances = sigma;
	} 

	/**
		Adds a constants to the likelihood.
	*/
	public void enableConstant()
	{
		if (!includeLikelihoodConstant)
		{
			includeLikelihoodConstant = true;
			if (init) 
				likelihood += constant;
		}
	}

	/**
		Subtracts a constants to the likelihood.
	*/
	public void disableConstant()
	{
		if (includeLikelihoodConstant)
		{
			includeLikelihoodConstant = false;
			if (init) 
				likelihood -= constant;
		}
	}
	
	/**
		Calculate the residuals of the linear model
		
		@return the residuals
	*/
	private double[] calculateResiduals() 
	{
		double[] residuals = new double[dim];
			
		for (int i = 0; i < dim; i++) 
		{
			for (int j = 0; j < linearMD; j++)
				residuals[i] += linearModel[i][j] * means[j];
			residuals[i] -= scores[i];
		}
		
		return residuals;
	}

	/**
		Calculate the covariance matrix
		
		@return the covariance matrix
	*/
	private double[][] calculateCovariances() 
	{
		double[][] varMatrix = new double[dim][dim];
		
		for (int j = 0; j < vcCount; j++) 
			for (int r = 0; r < dim; r++) 
				for (int c = r; c < dim; c++)
				{		
					varMatrix[r][c] += variances[j] * varComponents[j][r][c];
					varMatrix[c][r] = varMatrix[r][c]; //Decomposer wants a symmetric matrix
				}
		
		return varMatrix;
	}

	/**
		Evaluates the equation.
		
		@return the likelihood of the system
		@throws InfiniteLikelihoodException if the likelihood (that actually is the determinant of the
				decomposition) is Infinity
	*/
	public double evaluate() throws InfiniteLikelihoodException
	{
		likelihood = includeLikelihoodConstant ? constant : 0.0;
		
		//solves the equation x = varMatrix * residuals
		double[]   residuals = calculateResiduals();
		double[][] varMatrix = calculateCovariances();
		
		likelihood -= 0.5 * Utilities.innerProduct(residuals, decomposition.solve(varMatrix, residuals));
		likelihood -= decomposition.lnDet();
		
		//There may be problems when computing the determinant of a large matrix in floating point arithmetic,
		//due to accuracy issues, and this generates -/+Infinity values
		//See here for details: https://scicomp.stackexchange.com/questions/23267/is-there-any-rapid-way-to-calculate-the-determinant-of-nxn-covariance-matrix/23268#23268
		if (Double.isInfinite(Math.abs(likelihood))) 
			throw new InfiniteLikelihoodException("Warning : matrix decomposition failed");
				
		init = true;
		return likelihood;
	}

	/**
		Removes the c-th predictor from the system
		
		@param c the position
	*/
	public void deleteColumn(int c) 
	{
		int nr = linearModel.length;
		int nc = linearModel[0].length - 1;
		double[][] tmpLinearModel = new double[nr][nc];
		for (int i = 0; i < nr; i++) 
		{
			int count = 0;
			for (int j = 0; j < nc; j++) 
				if (j != c)
					tmpLinearModel[i][count++] = linearModel[i][j];
		}
		
		linearModel = tmpLinearModel;
		means = new double[nc];
	}
	
}


/**
	Runtime exception. 

	Raised when the determinat of the decomposition matrix cannot be evaluated
	
	@author      Alessia Visconti <alessia.visconti@kcl.ac.uk>
	@version     1.0             
 */
class InfiniteLikelihoodException extends RuntimeException
{
    private static final long serialVersionUID = 1L;
	/**
		Constructor. 
		
		Initialises the exception message.
		
		@param msg the exception message
	*/
	public InfiniteLikelihoodException(String msg) 
	{
		super(msg);
	}
}

