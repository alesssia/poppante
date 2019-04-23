/*
 * 	 Adaptive.java
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

package com.github.alesssia.probutils;

/**
	Implements an adaptive permutation approach. 

	References for the implementation can be found at:
		Che, R., Jack, J. R., Motsinger-Reif, A. A., and
	    Brown, C. C. (2014). An adaptive permutation approach 
	    for genome-wide association study: evaluation and 
	    recommendations for use. BioData Min, 7(9).
	and at: 
		Besag, J., and Clifford, P. (1991). Sequential monte 
		carlo p-values. Biometrika, 78(2), 301-304.

	
	@author      Alessia Visconti <alessia.visconti@kcl.ac.uk>
	@version     1.0     
    @param <T extends Permutable> the type of the value being boxed
*/

public class Adaptive<T extends Permutable>
{
	/** p-value that controls the experiment-wise error rate (EWER). */
	private final double alpha;
	/** Significance of the threshold that equals the standard error in estimation,
		when the desired p-value is equal to alpha.  */
	private final double c;
	/** Object to which the permutation test is performed.
		
		It must implement the permute function. This function
		performs the permutation and returns the test statistic.
	*/
	private final T obj;	

	/**
		Constructor.
		
		@param alpha p-value one would like to estimate
		@param c desired precision level of p-value estimation
		@param obj the object to assess
	 */
	public Adaptive(double alpha, double c, T obj)
	{
		this.alpha = alpha;
		this.c = c;
		this.obj = obj;
	}
	
	/**
		Performs the adaptive permutation.
		
		@param pstar the p-value of the test statistic
		@return the p-value for the adaptive permutation
        @param <T extends Permutable> the type of the value being boxed
	 */
	public <T> double adapt(double pstar)
	{
		//selects the permutation parameters
		int r = r();
		int b = b();
		
		
		int R = 0; //counts the number of successes
		for(int B=1; B<=b; B++)
		{
			if (obj.permute() <= pstar)
				R++;
		
			if (R == r)
				return ((double)R)/B;
		}
		
		return ((double)(R+1))/(b+1);
	}
	

	/**
		Chooses the maximal number of permutations.
	
		It works for either adaptive or standard permutation.
		The number of permutations is evaluated as a function
		of the p-value one would like to estimate (alpha) and 
		of desired precision level of p-value estimation (c)
		where the standard error of the estimated p-value is 
		c*alpha.

		@return the maximal number of permutations
	*/	
	private int b() 
	{		
		return (int)Math.ceil( (alpha*(1-alpha)) / (Math.pow(alpha*c, 2)) );
	}


	/**
		Determines the number of tests that should be sampled in
		adaptive permutation.
	
		The number of tests is chosen such that the desired level 
		of significance (alpha) will be sampled with a 68% CI 
		(about 1 standard error) contained within a specified 
		level of precision (c).

		@return the maximal number of permutations
	*/	
	private int r()
	{
		double error = alpha * c;
		int R = 0;
		boolean done = false;
		while (!done)
		{
			R++;
			MyNegativeBinomialDistribution nb = new MyNegativeBinomialDistribution(R, alpha);
			double a = R / (nb.qnbinom(0.1586553) + R);
			double b = R / (nb.qnbinom(0.8413447) + R);
			if (Math.max(Math.abs(a-alpha), Math.abs(b-alpha)) < error)
				 done = true;
		}
	
		return R;
	}
	
}
