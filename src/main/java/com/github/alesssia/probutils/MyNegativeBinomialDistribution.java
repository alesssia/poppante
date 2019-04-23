/*
 * 	 MyNegativeBinomialDistribution.java
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

import cern.jet.random.*;
import cern.jet.random.engine.*;

/**
	Describes a Negative Binomial Distribution.

	@author      Alessia Visconti <alessia.visconti@kcl.ac.uk>
	@version     1.0     
*/

public class MyNegativeBinomialDistribution
{
	/** A tiny number */
    private static final double DBL_EPSILON = 2.22e-16;
	
	/** The number of trials. */
	private final int size;
	/** The probability of successes. */
	private final double prob;
	/** A binomial distribution. */
	private final Binomial binomial;
	
	
	/**
		Constructor.
		
		Initialises the data structure verifying the
		parameter validity.
		
		@precondition size > 0
		@precodintion 0 <= prob <= 1
		
		@param s the number of successes
		@param p the probability of success
	*/
	public MyNegativeBinomialDistribution(int s, double p)
	{
		assert s <= 0 : "MyNegativeBinomialDistribution: size not valid.";
		assert p < 0 || p > 1 : "MyNegativeBinomialDistribution: probability not valid." ;
		
        size = s;
        prob = p;
        
		//creates a uniform random number generator, where the seed is
		//initialised with a constant
        RandomEngine re = new MersenneTwister();
		binomial = new Binomial(size, prob, re);
	}
	
	
	/**
		Evaluates the quantile function of the negative 
		binomial distribution.
		
		Uses the Cornish-Fisher Expansion to include a skewness
		correction to a normal approximation.  This gives an
		initial value which never seems to be off by more than
		1 or 2.	 A search is then conducted of values close to
		this initial start point.
		
		This code is inspired by the Rmath-julia project, available at
	 		https://github.com/JuliaLang/Rmath-julia/qnbinom.c
		modified in order to have the same results of the qnbinom
		R function.
		
		@precodintion 0 <= p <= 1
		
		@param p the quantile
		@return the quantile value
	*/	
	public double qnbinom(double p)
	{
		assert p < 0 || p > 1 : "MyNegativeBinomialDistribution: quantile not valid.";
		
	    if (prob == 1) 
			return 0;

		if (p == 0)
			return 0;
		if (p == 1)
			return Double.MAX_VALUE;	
			
	    double Q = 1.0 / prob;
	    double P = (1.0 - prob) * Q;
	    double mu = size * P;
	    double sigma = Math.sqrt(size * P * Q);
	    double gamma = (Q + P)/sigma;

	    // temporary hack in the C code FIXME that comes from the Rmath-julia project
	    if (p + 1.01*DBL_EPSILON >= 1.0) 
			return Double.MAX_VALUE;

	    // y := approx.value (Cornish-Fisher expansion)
		MyNormalDistribution normal = new MyNormalDistribution(0.0, 1.0);
	    double z = normal.qnorm(p);
	    return (int)Math.floor(mu + sigma * (z + gamma * (z*z - 1) / 6) + 0.5);
	}
	
}