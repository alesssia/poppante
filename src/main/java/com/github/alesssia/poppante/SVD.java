/*
 * 	 SVD.java
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
	Performs the Single Value Decomposition (SVD) 


	References for the implementation can be found at:
		Press, William H., et al. Numerical recipes in C. 
		Vol. 2. Cambridge: Cambridge university press, 1996.

	@author      Alessia Visconti <alessia.visconti@kcl.ac.uk>
	@author      Mario Falchi     <mario.falchi@kcl.ac.uk>
	@author      Hariklia Eleftherohorinou <hariklia.eleftherohorinou06@imperial.ac.uk>
	@version     1.0     
*/

class SVD 
{
	/** Threshold for the zero SVD values */
	public static final double TOL = 1.0e-6; 	
	/** Maximum number of iterations */ 	
	public static final double MAX_ITERATION = 30; 

	/** Number of rows */ 
	private int m;	
	/** Number of parameters */ 
	private int n;	
	
	/** The matrix U */ 
	private double[][] u; 	
	/** The matrix V */
	private double[][] v; 	 // (not the transpose V^T)
							 // is output as v[1..n][1..n]
	/**  The solution vector */
	private double[] x; 
	/**  The diagonal matrix of singular values vector  */	
	private double[] w; 	 //w[1..n]
	
	
	/**
		Returns the value of i-th element of the 
		diagonal matrix of singular values.
		
		@param i the index
		@return the element of the diagonal matrix
	*/ 	
	public double w(int i)
	{
		return w[i];
	}
	
	/**
		Returns the solution vector (after the back-substitution 
		obtained by the fit).
		
		@return the solution vector
	*/
	public double[] x()
	{
		return x;
	}
	
	
	/**
		Evaluates the Single Value Decomposition.
		
		@throws NotConvergencingException if do not convergence in 30 iterations
		@param a the linear system
		@param mp number of parameters
		@param np number of rows
	*/
	public void decompose(double[][] a, int mp, int np) throws NotConvergencingException
	{
		int flag, i, its, j, jj, k, l, nm;
		double c, f, g, h, s, scale, x1, y, z;
		double anorm;
		double[] rv1;

		empty();

		l = nm = 0;
		m = (mp == -1) ? a.length : mp;
		n = (np == -1) ? a[0].length : np;

		u = new double[m][n];
		for (int ir = 0; ir < m; ir++)
			System.arraycopy(a[ir], 0, u[ir], 0, n);

		v = new double[n][n];
		w = new double[n];
		
		rv1 = new double[n];

		g = scale = anorm = 0.0; // Householder reduction to bidiagonal form
		
		for (i = 0; i < n; i++) 
		{
			l = i + 1;
			rv1[i] = scale * g;
			g = s = scale = 0.0;
			
			if (i < m) 
			{
				for (k = i; k < m; k++)
					scale += Math.abs(u[k][i]);
					
				if (scale > Utilities.FPMIN) 
				{
					for (k = i; k < m; k++) 
					{
						u[k][i] /= scale;
						s += u[k][i] * u[k][i];
					}
					f = u[i][i];
					
					g = f < 0.0 ? Math.sqrt(s) : (-Math.sqrt(s));
					h = f * g - s;
					u[i][i] = f - g;
					for (j = l; j < n; j++) 
					{
						s = 0.0;
						for (k = i; k < m; k++)
							s += u[k][i] * u[k][j];
						f = s / h;
						for (k = i; k < m; k++)
							u[k][j] += f * u[k][i];
					}
					for (k = i; k < m; k++)
						u[k][i] *= scale;
				}
			}
			
			w[i] = scale * g;
			g = s = scale = 0.0;
			
			if ((i < m) && (i != n - 1)) 
			{
				for (k = l; k < n; k++)
					scale += Math.abs(u[i][k]);
				if (scale > Utilities.FPMIN) 
				{
					for (k = l; k < n; k++) 
					{
						u[i][k] /= scale;
						s += u[i][k] * u[i][k];
					}
					f = u[i][l];
					g = f < 0.0 ? Math.sqrt(s) : (-Math.sqrt(s));
					h = f * g - s;
					u[i][l] = f - g;
					for (k = l; k < n; k++)
						rv1[k] = u[i][k] / h;
					for (j = l; j < m; j++) 
					{
						s=0.0;
						for (k = l; k < n; k++)
							s += u[j][k] * u[i][k];
						for (k = l; k < n; k++)
							u[j][k] += s * rv1[k];
					}
					for (k = l; k < n; k++)
						u[i][k] *= scale;
				}
			}
			
			x1 = Math.abs(w[i]) + Math.abs(rv1[i]);
			if (anorm < x1)
				anorm = x1;
		}
		
		// accumulation of right-hand transformations
		for (i = n - 1; i >= 0; i--) 
		{
			if (i < n - 1) 
			{
				if (Math.abs(g) > Utilities.FPMIN) 
				{
					// double division to avoid possible underflow
					for (j = l; j < n; j++)
						v[j][i] = (u[i][j] / u[i][l]) / g;
					for (j = l; j < n; j++) 
					{
						s = 0.0;
						for (k = l; k < n; k++)
							s += u[i][k] * v[k][j];
						for (k = l; k < n; k++)
							v[k][j] += s * v[k][i];
					}
				}
				
				for (j = l; j < n; j++)
					v[i][j] = v[j][i] = 0.0;
			}
			
			v[i][i] = 1.0;
			g = rv1[i];
			l = i;
		}

		// accumulation of left-hand transformations
		for (i = (Math.min(n, m)) - 1; i >= 0; i--) // piu' complicata per le covariate hks2009
		{
			l = i + 1;
			g = w[i];
			for (j = l; j < n; j++)
				u[i][j] = 0.0;
			
			if (Math.abs(g) > Utilities.FPMIN) 
			{
				g = 1.0 / g;
				for (j = l; j < n; j++) 
				{
					s = 0.0;
					for (k = l; k < m; k++)
						s += u[k][i] * u[k][j];
					f = (s / u[i][i]) * g;
					
					for (k = i; k < m; k++)
						u[k][j] += f * u[k][i];
				}
				for (j = i; j < m; j++)
					u[j][i] *= g;
			}
			 else
				for (j = i; j < m; j++)
					u[j][i] = 0.0;
			++u[i][i];
		}

		// Diagonalization of the bi-diagonal form:
		// Loop over singular values and over allowed iterations
		for (k = n - 1; k >= 0; k--) 
		{
			for (its = 1; its <= 30; its++) 
			{
				flag = 1;
				for (l = k; l >= 0; l--) 
				{
					// Test for splitting
					// note that rv1[1] is always zero.
					nm = l - 1;
					if ((Math.abs(rv1[l]) + anorm) == anorm) 
					{
						flag = 0;
						break;
					}
					if ( (Math.abs(w[nm]) + anorm) == anorm)
						break;
				}
				
				if (flag == 1) 
				{
					c = 0.0; // cancellation of rv1[l], if l > 1
					s = 1.0;
					for (i = l; i < k; i++) 
					{
						f = s * rv1[i];
						rv1[i] = c * rv1[i];
						if ( (Math.abs(f) + anorm) == anorm)
							break;
					
						g = w[i];
						h = pythag(f, g);
						w[i] = h;
						h = 1.0 / h;
						c = g * h;
						s = (-f * h);
						for (j = 0; j < m; j++) 
						{
							y = u[j][nm];
							z = u[j][i];
							u[j][nm] = y*c + z*s;
							u[j][i] = z*c - y*s;
						}
					}
				}
				
				z = w[k];
				if (l == k) // Convergence
				{
					if (z < 0.0) // Singular value is made nonnegative
					{
						w[k] = -z;
						for (j = 0; j < n; j++)
							v[j][k] = -v[j][k];
					}
					break;
				}

				if (its == MAX_ITERATION)
					throw new NotConvergencingException("SVD: No convergence in " + MAX_ITERATION + " iterations");
					
					
				// shift from bottom 2-by-2 minor
				x1 = w[l];
				nm = k - 1;
				y = w[nm];
				g = rv1[nm];
				h = rv1[k];
				f = ((y-z) * (y+z) + (g-h) * (g+h)) / (2.0*h*y);
				g = pythag(f, 1.0);
				f = ((x1-z)*(x1+z) + h*((y/(f+(f<0?-g:g)))-h))/x1;
				c = s = 1.0;
			
				// next QR transformation
				for (j = l; j <= nm; j++) 
				{
					i = j + 1;
					g = rv1[i];
					y = w[i];
					h = s*g;
					g = c*g;
					z = pythag(f, h);
					rv1[j] = z;
					c = f/z;
					s = h/z;
					f = x1*c + g*s;
					g = g*c - x1*s;
					h = y*s;
					y *= c;
					for (jj = 0; jj < n; jj++) 
					{
						x1 = v[jj][j];
						z = v[jj][i];
						v[jj][j] = x1*c + z*s;
						v[jj][i] = z*c - x1*s;
					}
					z = pythag(f, h);
					w[j] = z; // Rotation can be arbitrary if z=0
					
					if (Math.abs(z) > Utilities.FPMIN) 
					{
						z = 1.0 / z;
						c = f*z;
						s = h*z;
					}
					f = c*g + s*y;
					x1 = c*y - s*g;
					for (jj = 0; jj < m; jj++) 
					{
						y = u[jj][j];
						z = u[jj][i];
						u[jj][j] = y*c + z*s;
						u[jj][i] = z*c - y*s;
					}
				}
				
				rv1[l] = 0.0;
				rv1[k] = f;
				w[k] = x1;
			}
		}
	}

	public void edit() 
	{
		int j;
		double wmax, wmin;
		wmax = 0.0;

		for (j = 0; j < n; j++)
			if (w[j] > wmax)
				wmax = w[j];

		wmin = wmax * TOL;

		for (j = 0; j < n; j++)
			if (w[j] < wmin)
				w[j] = 0.0;
		
	}

	/**
		Solves the linear system of equations Ax=b.
		
		@param b the vector b
	*/	
	public void backSubst(double[] b) 
	{
		int jj, j, i;
		double s;
		double[] tmp = new double[n];
		x = new double[n];

		// calculate U^T * B
		for (j = 0; j < n; j++) 
		{
			s = 0.0;
			if (w[j] != 0) // Nonzero result only if wj is nonzero
			{
				for (i = 0; i < m; i++)
					s += u[i][j] * b[i];
				s /= w[j];
			}
			tmp[j] = s;
		}

		// Matrix multiply by V to get answer
		for (j = 0; j < n; j++) 
		{
			s = 0.0;
			for (jj = 0; jj < n; jj++)
				s += v[j][jj] * tmp[jj];
			x[j] = s;
		}
	}

	public double RSS(double[][] M, double[] b) 
	{
		double rss = 0.0;
		for (int i = 0; i < m; i++) 
		{
			double partial = 0.0;
			for (int j = 0; j < n; j++)
				partial += x[j] * M[i][j];
			partial -= b[i];
			rss += partial * partial;
		}
		return rss;
	}

	/**
		Sets the dimensions to zero.
		
		It corresponds to create an empty object.
	*/
	private void empty() 
	{
		m = 0;
                n = 0;
	}

	/**
		Calculates sqrt(a*a + b*b) safely.
		
		@param a the first value
		@param b the second value
		@return sqrt(a*a + b*b)
	*/
	private double pythag(double a, double b) {
		double absa, absb;
		absa = Math.abs(a);
		absb = Math.abs(b);
		if (absa > absb)
			return absa * Math.sqrt(1.0 + (absb / absa) * (absb / absa));
		else
			return (absb == 0.0 ? 0.0 : absb * Math.sqrt(1.0 + (absa / absb) * (absa / absb)));
	}

}
