/*
 * 	 MyNormalDistribution.java
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
	Describes a Normal Distribution.

	@author      Alessia Visconti <alessia.visconti@kcl.ac.uk>
	@version     1.0     
*/

public class MyNormalDistribution
{
	/** A tiny number */
    public static final double DBL_EPSILON = 2.22e-16;
    
	/** Mean the Normal Distribution */
	double mu;
	/** Standard deviation of the Normal Distribution */
	double sigma;

	
	/**
		Constructor.
		
		Initialises the data structure.
		
		@precondition s >= 0
		@param m the mean
		@param s the standard deviation
	*/	
	public MyNormalDistribution(double m, double s)
	{
		assert (s < 0.0) : "MyNormalDistribution: sigma not valid.";
		
		sigma = s;
		mu = m;
	}
	
	
	/**
		Evaluates the quantile function of the negative 
		binomial distribution.
		
		Uses the Cornish-Fisher Expansion to include a skewness
		correction to a normal approximation.  This gives an
		initial value which never seems to be off by more than
		1 or 2.	 A search is then conducted of values close to
		this initial start point.
		
		This code is inspired by the Rmath-julia project, available at: 		
			https://github.com/JuliaLang/Rmath-julia/qnorm.c
		
		@precondition 0 <= p <= 1
		
		@param p the quantile
		@return the quantile value
		@throws RuntimeException if the p-value is outside range in Normal Inverse Transformation
	*/	
	public double qnorm(double p) throws RuntimeException
	{		
		assert (p < 0.0 || p > 1.0) : "MyNormalDistribution: quantile not valid.";
		
		if (p == 0)
			return Double.MIN_VALUE;
		if (p == 1)
			return Double.MAX_VALUE;
		
		if(sigma == 0)	
			return mu;
		
		return ninv(p);
	}
	
	/**
		Evaluates the quantile normalisation (inverse normal transformation). 
		
		This code is the java Equivalent of: 
			Wichura, M. J. (1988) Algorithm AS 241: The percentage points of
			the normal distribution.  _Applied Statistics_, *37*, 477-484.
		It has been translated from the C Equivalent of an unknown developer.
		
		@param p the value to be transformed
		@return the transformed value
		@throws RuntimeException if the p-value is outside range in Normal Inverse Transformation
	*/	
	public static double ninv(double p) throws RuntimeException
	{
		
	   double SPLIT1 = 0.425;
	   double SPLIT2 = 5.0;
	   double CONST1 = 0.180625;
	   double CONST2 = 1.6;
	   double HUGE_VAL = Double.MAX_VALUE;
	
	   double a[] = new double[] { 3.3871328727963666080E0, 1.3314166789178437745E2, 1.9715909503065514427E3,1.3731693765509461125E4, 4.5921953931549871457E4, 6.7265770927008700853E4, 3.3430575583588128105E4, 2.5090809287301226727E3 };
	
	   double b[] = new double[]{4.2313330701600911252E1, 6.8718700749205790830E2,  5.3941960214247511077E3,  2.1213794301586595867E4, 3.9307895800092710610E4, 2.8729085735721942674E4, 5.2264952788528545610E3 };
	
	   double c[] = new double[] { 1.42343711074968357734E0, 4.63033784615654529590E0, 5.76949722146069140550E0, 3.64784832476320460504E0, 1.27045825245236838258E0, 2.41780725177450611770E-1, 2.27238449892691845833E-2, 7.74545014278341407640E-4 };
	
	   double d[] = new double[]{ 2.05319162663775882187E0, 1.67638483018380384940E0, 6.89767334985100004550E-1, 1.48103976427480074590E-1, 1.51986665636164571966E-2, 5.47593808499534494600E-4, 1.05075007164441684324E-9 } ;
	
	  double e[] = new double[]{ 6.65790464350110377720E0, 5.46378491116411436990E0, 1.78482653991729133580E0, 2.96560571828504891230E-1, 2.65321895265761230930E-2, 1.24266094738807843860E-3, 2.71155556874348757815E-5, 2.01033439929228813265E-7 } ;
	
	   double f[] = new double[] { 5.99832206555887937690E-1, 1.36929880922735805310E-1, 1.48753612908506148525E-2, 7.86869131145613259100E-4, 1.84631831751005468180E-5,  1.42151175831644588870E-7, 2.04426310338993978564E-15 } ;
	
	   double q = p - 0.5;
	   double r;
	   double x = 0.0;
	
	   if ( Math.abs( q ) < SPLIT1 )
	   {
	      r = CONST1 - q * q ;
	      return q * ((((((( a[7] * r + a[6] ) * r + a[5] ) * r + a[4] ) * r + a[3] ) * r + a[2] ) * r + a[1] ) * r + a[0] ) /  ((((((( b[6] * r + b[5] ) * r + b[4] ) * r + b[3] ) * r  + b[2] ) * r + b[1] ) * r + b[0] ) * r + 1.0 ) ;
	   }
	   else
	   {
			if ( q < 0 )
				r = p ;
			else
				r = 1.0 - p ;
	
			if ( r < 1e-200)
				throw new RuntimeException("ERROR: p-value " + r + " outside range in Normal Inverse Transformation");
	
			if ( r > 0.0 )
			{
				r = Math.sqrt( -Math.log( r ) ) ;
				if ( r <= SPLIT2 )
				{
					r -= CONST2 ;
					x = ((((((( c[7] * r + c[6] ) * r + c[5] ) * r + c[4] ) * r + c[3] ) * r + c[2] ) * r + c[1] ) * r + c[0] ) /  ((((((( d[6] * r + d[5] ) * r + d[4] ) * r + d[3] ) * r + d[2] ) * r + d[1] ) * r + d[0] ) * r + 1.0 ) ;
				}
				else
				{
					r -= SPLIT2 ;
					x =  ((((((( e[7] * r + e[6] ) * r + e[5] ) * r + e[4] ) * r + e[3] ) * r + e[2] ) * r + e[1] ) * r + e[0] ) /  ((((((( f[6] * r + f[5] ) * r + f[4] ) * r + f[3] ) * r + f[2] ) * r + f[1] ) * r + f[0] ) * r + 1.0 ) ;
				}
			 }
			else
				x = HUGE_VAL;
			
			if ( q < 0 )
				x = -x ;

			return x ;
	   }
	}
	
}