/*
 * 	 NormalSet.java
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

/**
	Represents a set of linear mixed models within a variant 
	component framework.
	
	Each family correspond to a different linear mixed model,
	that is represented by a single equation.

	References for the implementation can be found at:
		Press, William H., et al. Numerical recipes in C. 
		Vol. 2. Cambridge: Cambridge university press, 1996.

	@author      Alessia Visconti <alessia.visconti@kcl.ac.uk>
	@author      Mario Falchi     <mario.falchi@kcl.ac.uk>
	@author      Hariklia Eleftherohorinou <hariklia.eleftherohorinou06@imperial.ac.uk>
	@version     1.0     
	@see com.github.alesssia.poppante.Normal	
*/

class NormalSet 
{
	
	 /** The fractional convergence tolerance to be achieved in the Amoeba minimiser. */	 
	public static final double PRECISION = 1.0E-8; 	

	/** Number of linear equations in the set */
	private final int size;			
	/** Number of fixed effects in the model.
		
		The number of fixed effects is the same in each
		linear equation that belong to the set. */
	private int linearMD;		
	/** Number of random effects in the model.
		
		The number of random effects is the same in each
		linear equation that belong to the set. */
	private final int vcCount;		
	
	/** Linear mixed models.
		
		Each of them describes one of the family in
		the analysis. */
	private final Normal[] sets;
	/** Variances. */
	private double[] variances; 
	/** Beta coefficients */
	private double[] means;		
	/** Likelihood of the model */
	private double likelihood;	


	/**
		Constructor. 
		
		It initialises the data structure.
		
		@param setCount number of equations.
		@param lmd number of fixed effects
		@param vc number of random effects
	*/
	public NormalSet(int setCount, int lmd, int vc) 
	{
		size = setCount;
		linearMD = lmd;
		vcCount = vc;
		
		sets = new Normal[size];
		for (int i = 0; i < size; i++)
			sets[i] = new Normal(vcCount, linearMD);
		
		means = new double[linearMD]; 
		variances = null;
		
		likelihood = 0.0;
	}
	
	/**
		Returns the equation at the given index.
		
		@param index the index
		@return the equation
	*/	
	public Normal sets(int index) 
	{
		return sets[index];
	}
	
	/**
		Returns the number of equations.
		
		@return the number of equations
	*/
	public int size()
	{
		return size;
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
		Returns the variance at the given position.
		
		@param i the position
		@return the variance
	*/	
	public double variances(int i)
	{
	     return variances[i];
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
		Returns the beta coefficient at the given position.
		
		@param i the position
		@return the beta coefficient
	*/	
	public double means(int i)
	{
	     return means[i];
	}
		
	/**
		Solves the set of equations.
		
		It uses the Nelder–Mead method (aka Amoeba method) to
		minimise the equations in the multi-dimensional space. 
		Its convergence is not guarantee.
		
		@throws NotConvergencingException if the Amoeba minimiser or the SVD decomposition can't converge
		@see com.github.alesssia.poppante.Amoeba
		@see com.github.alesssia.poppante.SVD
		@see com.github.alesssia.poppante.Normal
	*/
	public void solve() throws NotConvergencingException
	{
		editLinearDegenerates();
			
		//Amoeba is initialised to have a size of 
		//#variables in the liner system + #components in the VC
		Amoeba solver = new Amoeba(this, countParameters());
		
		//directions are set to have diagonal 1
		solver.reset(1);
		solver.setPoint(getStartingPoint());  
		
		solver.minimize(PRECISION); 
		
		double scale = 2.0;
		int counter = 0;
		double currentmin;
		double lastmin;
		
		do 
		{		
			double varMax = solver.point(linearMD);

			for (int i = linearMD + vcCount - 1; i > linearMD; i--) 
				if (solver.point(i) > varMax)
					varMax = solver.point(i);
			
			for (int i = linearMD + vcCount - 1; i >= linearMD; i--) 
				if (solver.point(i) < (varMax - 5.0))
					solver.setPoint(i, (varMax - 5.0));
			
			lastmin = solver.fmin();
			
			scale *= 0.5;
			
			solver.reset(scale);
			
			solver.minimize(PRECISION);
			counter++;
			currentmin = solver.fmin();
			
		} while (currentmin > PRECISION && ((lastmin-currentmin)/currentmin) > PRECISION);

		selectPoint(solver.point());
	}

	/**
		Edits the linear degenerates.
		
		@throws NotConvergencingException if the SVD decomposition can't converge
		@see com.github.alesssia.poppante.SVD
	*/
	private void editLinearDegenerates() throws NotConvergencingException
	{
		int rows = countObservations();
		int cols = linearMD;
		double[][] m = new double[rows][cols]; 	
		double[] b = new double[rows]; 			
		
		int total = 0;
		for (int i = 0; i < size; i++) 
		{ 
			for (int j = 0; j < sets[i].dim(); j++) 
			{ 
				b[total] = sets[i].scores(j);
				m[total] = sets[i].linearModel(j);
				total++;
			}
		}
		
		SVD engine = new SVD();
		
		int c = 1;
		while (c <= cols) 
		{	
			engine.decompose(m, rows, c);
						
			engine.edit();
			boolean redundant = false;
			
			for (int j = 0; j < c; j++) 
				if (engine.w(j) == 0.0)
					redundant = true;
			
			if (redundant) 
			{
				for (int i = 0; i < size; i++)
					sets[i].deleteColumn(c-1);

				double[][] tmpM = new double[rows][cols - 1];
				
				for (int i = 0; i < rows; i++) 
				{
					int count = 0;
					for (int j = 0; j < cols; j++) 
					{
						if (j != c - 1)
							tmpM[i][count++] = m[i][j];
					}
				}
				
				m = tmpM;
				cols--; 
				
				if (c > cols)
					engine.decompose(m, rows, cols);
					
			} 
			else
				c++;
		}
					
		variances = new double[vcCount];
		engine.backSubst(b);
		
		int tmp = means.length; //engine may change the lenght of the means vector
		means = engine.x(); 
		if (means.length != tmp)
		{
			tmp = means.length;
			linearMD = tmp; 
			for (int i = 0; i < size; i++)
				sets[i].resizeLM(tmp);
		}

 		double residualVar = (Math.abs(engine.RSS(m, b)) + Utilities.FPMIN) / rows;
		double[] var_weights = new double[vcCount]; 
		int[] var_counts = new int[vcCount]; 
		double[] scaled_variances = new double[vcCount];

		for (int i = 0; i < size; i++) 
		{
			for (int j = 0; j < sets[i].dim(); j++) 
			{
				double scaled_variances_sum = 0.0;
				for (int v = 0; v < vcCount; v++) 
				{
					if (sets[i].varComponents(v).length > 0)
						scaled_variances[v] = sets[i].varComponents(v,j,j); 
					else
						scaled_variances[v] = 0.0;
					
					scaled_variances_sum += scaled_variances[v];
				}
				
				double value = 1.0 / scaled_variances_sum;
				for (int v = 0; v < vcCount; v++)
					scaled_variances[v] *= value;
				
				for (int v = 0; v < vcCount; v++) {
					if (scaled_variances[v] > 0.0) 
					{
						var_weights[v] += scaled_variances[v];
						var_counts[v]++;
					}
				}
			}	
		}
		
		for (int v = 0; v < vcCount; v++) 
			variances[v] = residualVar * var_weights[v] / var_counts[v];

	}

	/**
		Selects the starting points for the Amoeba minimiser.
		
		The starting points correspond to the current beta 
		coefficients and to the logarithm of the current
		variances.
		
		@return the starting point
		@see com.github.alesssia.poppante.Amoeba
	*/
	private double[] getStartingPoint() 
	{
		
		double[] startPoint = new double[linearMD + vcCount];
        System.arraycopy(means, 0, startPoint, 0, linearMD);
		
		for(int i=0; i<vcCount; i++)
			startPoint[i+linearMD] = Math.log(variances[i]);

		return startPoint;	
	}
	
	/**
		Sets the beta coefficients and the variances to the 
		values resulting from the Amoeba minimiser. 
		
		These set values are propagated to each equation
		belonging to the set.
		
		Variances greater than 16 are set to 1.0E7 while 
		variances smaller than -16 are set to 1.0E-7.
		
		@param v the starting point
		@see com.github.alesssia.poppante.Amoeba
	*/
	private void selectPoint(double[] v) 
	{	
		int j = linearMD;
        System.arraycopy(v, 0, means, 0, j);
		
		for(int i=0; i<vcCount; i++)
		{
			if(v[j] > 16 || v[j] < -16)
				variances[i] = v[j] > 0 ? 1.0E7 : 1.0E-7;
			else
				variances[i] = Math.exp(v[j]);
			j++;
		}
		
		for(int i=0; i<size; i++)
			sets[i].setParameters(means, variances);
	}


	/**
		Returns the number of parameters in the systems.
		
		It correspond to fixed and random effects.
		
		@return the number of parameters
	*/
	public int countParameters() 
	{
		return linearMD + vcCount;
	}

	/**
		Returns the number of observations in the systems.
		
		It counts the observations that belong to each
		equation.
		
		@return the number of observations
	*/
	private int countObservations() 
	{
		int rows = 0;
		for (int i = 0; i < size; i++)
			rows += sets[i].dim();
		
		return rows;
	}

	/**
		Calls the method evaluate after having set the beta 
		coefficients and the variances 
		
		@param point the new beta coefficients and variances
		@return the likelihood of the system
		@see com.github.alesssia.poppante.Amoeba
	*/
	public double callEvaluate(double[] point)
	{
		selectPoint(point);
		return evaluate();
	}
	
	/**
		Adds a constants to the likelihood of the equations
		belonging to the system.

		@see com.github.alesssia.poppante.Normal
	*/
	public void enableConstant() 
	{
		for (int i = 0; i < size; i++)
			sets[i].enableConstant();
	}
	
	/**
		Substracts a constants to the likelihood of the equations
		belonging to the system.

		@see com.github.alesssia.poppante.Normal
	*/
	public void disableConstant() 
	{
		for (int i = 0; i < size; i++)
			sets[i].disableConstant();
	}
	
	/**
		Evaluates the equations belonging to the
		system and returns the likelihood of the
		solutions.		
		
		@return the likelihood of the system
		@see com.github.alesssia.poppante.Normal
	*/
	public double evaluate()
	{
		likelihood = 0.0;
		for (int i = 0; i < size; i++)
			likelihood += sets[i].evaluate();   //operator is NORMAL_MUL_LK (see QTDT implementation)
												// sum of the logLikelihood of the single model
		
		return (-likelihood);
	}
	
	/**
		Returns the -likelihood of every equations
		belonging to the system
		
		@return the -likelihood of every equations in the system
		@see com.github.alesssia.poppante.Normal
		@see com.github.alesssia.poppante.VC
	*/
	public double[] evaluateDetails() 
	{
		double[] LK = new double[size];

		for (int i = 0; i < size; i++)
			LK[i] = (-sets[i].likelihood());
		
		return LK;
	}

}

