/*
 * 	 Amoeba.java
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
 *   hariklia.eleftherohorinou06@imperial.ac.uk
 *   alessia.visconti@kcl.ac.uk
 */

package com.github.alesssia.poppante;

/** 
	Performs the Nelder–Mead method (aka Amoeba method) for minimising an 
	objective function in a multi-dimensional space. 

	The algorithm is provided with a starting guess, that is, an N-vector 
	of independent variables as the first point to try. The algorithm is 
	then supposed to make its own way	downhill through the unimaginable 
	complexity of an N-dimensional topography, until it encounters a (local, 
	at least) minimum.

	References for the implementation can be found at:
		Press, William H., et al. Numerical recipes in C. 
		Vol. 2. Cambridge: Cambridge university press, 1996.
	
	@author      Alessia Visconti <alessia.visconti@kcl.ac.uk>
	@author      Mario Falchi     <mario.falchi@kcl.ac.uk>
	@author      Hariklia Eleftherohorinou <hariklia.eleftherohorinou06@imperial.ac.uk>
	@version     1.0 
*/

class Amoeba 
{
	/** Maximum number of iterations for the minimiser */
	public static final long   CYCLEMAX  = 50000L;	 
	/**  Precision around zero */ 	
	public static final double ZEPS      = 3.0e-10;   
	
	/** The linear system */
	private double[][] directions;
	
	/** Geometrical representation of an N-dimension problem.
		
		If the problem has N- dimension the simplex has
		N+1 points */
	private final double[][] simplex;
	
	/** The variables to optimise */
	private double[] point;
	
	/** Points in the simplex */
	private final double[] y;
	
	/** Stores the simplex values */
	private double[] psum;
	
	/** Trial points */
	private final double[] ptry;
	
	/** The minimum likelihood reached */
	private double fmin;

	/**  Linear system to optimise */
	private final NormalSet caller;

	/**
		Constructor. 
		
		It initialises the data structures.
		
		@param caller the linear mixed model to optimise
		@param n the dimension of the system	
	*/
	public Amoeba(NormalSet caller, int n) 
	{
		//caller calls NormalSet.evaluate() to compute the 
		//overall likelihood 
		this.caller = caller;
		
		point = new double[n];
		directions = new double[n][n];
		
		simplex = new double[n+1][n];
		y = new double[n+1];
		psum = new double[n];
		ptry = new double[n];
		
		fmin = Double.MAX_VALUE;
	}
	
	/**
		Returns the points.
		
		@return the optimised parameters
	*/
	public double[] point()
	{
		return point;
	}
	
	/**
		Returns the i-th point.
		
		@param i the position
		@return the point
	*/
	public double point(int i)
	{
		return point[i];
	}
	
	/**
		Returns the minimum likelihood reached.
		
		@return the minimum likelihood
	*/
	public double fmin()
	{
		return fmin;
	}
	
	/**
		Sets the i-th point to value
		
		@param i the position
		@param value the value
	*/
	public void setPoint(int i, double value)
	{
		point[i] = value;
	}
	
	/**
		Setes the points to values
		
		@param values the values
	*/
	public void setPoint(double[] values)
	{
		point = values;
	}

	/**
		Resets the directions.
		
		It scales them and set the minimum to the maximum 
		possible value.
		
		@param scale the scale value for diagonal evaluation
	*/
	public void reset(double scale) 
	{
		int n = directions.length;
		
		//all values but the diagonal are set to 0
		directions = new double[n][n];
		for (int i = 0; i < n; i++) 
			directions[i][i] = 1.0 * scale;
			
		fmin = Double.MAX_VALUE;
	}

	/** 
		Minimises in N -dimensions by downhill simplex method.
		
		@param ftol the fractional convergence tolerance to be achieved in the function value
		@return the minimum value
		@throws NotConvergencingException if it can't converge
	*/
	public double minimize(double ftol) throws NotConvergencingException
	{
		if (point.length == 0)
        {
            fmin = caller.callEvaluate(point);
            return fmin;
        }
		int nvertex = point.length + 1;
		
		for (int i = 0; i < point.length; i++) 
		{
			for (int j = 0; j < simplex[i].length; j++)
				simplex[i][j] = point[j] + directions[i][j];
			
			//y[i] is the -likelihood of the NormalSet
			//where the the means are set to simplex[i]
			//and the variance to e^simplex[i]
			y[i] = caller.callEvaluate(simplex[i]);
			
			if (y[i] < fmin)
				fmin = y[i];
		}
	
		
		System.arraycopy(point, 0, simplex[nvertex-1], 0, point.length);
			
		y[nvertex-1] = caller.callEvaluate(simplex[nvertex-1]);		
		if (y[nvertex - 1] < fmin)
			fmin = y[nvertex - 1];
		
		long cycleCount = nvertex;
		
		psum = new double[psum.length];
		for (int m = 0; m < nvertex; m++) 
			for (int j = 0; j < psum.length; j++)
				psum[j] += simplex[m][j];
			
		while (true) 
		{
			
			int inhi, ilo, ihi;
			
			if (y[0] > y[1])
			{
				ilo = inhi = 1;
				ihi = 0;
			}
			else
			{
				ilo = inhi = 0;
				ihi = 1;
			}

			for (int i = 2; i < nvertex; i++) 
			{
				if (y[i] <= y[ilo])
					ilo = i;
				else if (y[i] > y[ihi]) 
				{
					inhi = ihi;
					ihi = i;
				} 
				else if (y[i] > y[inhi])
					inhi = i;
			}
					
			double rtol = 2.0 * Math.abs(y[ihi] - y[ilo]) / (Math.abs(y[ihi]) + Math.abs(y[ilo]) + ZEPS);
			
			if (rtol < ftol) 
			{
				System.arraycopy(simplex[ilo], 0, point, 0, point.length);
				fmin = y[ilo];
				return fmin;
			}
			
			if (cycleCount > CYCLEMAX)
				throw new NotConvergencingException("Warning : Amoeba couldn't converge in " + CYCLEMAX + " cycles");
			
			cycleCount += 2L;
			double ytry = amotry(ihi, -1.0);
				
			if (ytry <= y[ilo])
				amotry(ihi, 2.0);
			else if (ytry >= y[inhi]) 
			{
				double ysave = y[ihi];
				ytry = amotry(ihi, 0.5);
				if (ytry >= ysave) 
				{
					for (int i = 0; i < nvertex; i++) 
					{
						if (i != ilo) 
						{
							for (int s = 0; s < simplex[i].length; s++)
								simplex[i][s] += simplex[ilo][s];
							for (int s = 0; s < simplex[i].length; s++)
								simplex[i][s] *= 0.5;
							
							y[i] = caller.callEvaluate(simplex[i]);
						}
					}
					cycleCount += point.length;
					
					psum = new double[psum.length];
					for (int m = 0; m < nvertex; m++)
						for (int s = 0; s < simplex[m].length; s++)
							psum[s] += simplex[m][s];		
					
				}
			}
			else
				cycleCount--;
		
		}
	}

	/**
		Evaluate a trial point.
		
		It extrapolates by a factor through the face of the simplex across 
		from the high point, tries it, and replaces the high 
		point if the new point is better.
		
		@param ihi the highest (worst) point
		@param factor the factor
		@return the trial point
	*/
	private double amotry(int ihi, double factor)
	{
		double fac = (1.0 - factor) / point.length;
		
		for (int s = 0; s < psum.length; s++)
			ptry[s] = fac * psum[s];
		
		double value = factor-fac;
		for (int s = 0; s < simplex[ihi].length; s++)
			ptry[s] += value * simplex[ihi][s];
		
		double ytry = caller.callEvaluate(ptry);
		
		if (ytry < y[ihi]) 
		{
			y[ihi] = ytry;
			for (int s = 0; s < psum.length; s++) 
				psum[s] -= simplex[ihi][s];
			for (int s=0; s<psum.length; s++)
			{
				simplex[ihi][s] = ptry[s];
				psum[s] +=  simplex[ihi][s];
			}
		}
				
		return ytry;
	
	}
}


/**
	Runtime exception. 

	Raised when the a computation does not converge.

	@author      Alessia Visconti <alessia.visconti@kcl.ac.uk>
	@version     1.0             
 */
class NotConvergencingException extends RuntimeException
{
    private static final long serialVersionUID = 1L;
	/**
		Constructor. 
		
		Initialises the exception message.
		
		@param msg the exception message
	*/
	public NotConvergencingException(String msg) 
	{
		super(msg);
	}
}
