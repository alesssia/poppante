/*
 * 	 Utilities.java
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
 *   alessia.visconti@kcl.ac.uk
 */

package com.github.alesssia.poppante;

import java.util.*;
import java.text.*;
import Jama.*;

/**
	Helper class.

	@author      Alessia Visconti <alessia.visconti@kcl.ac.uk>
	@version     1.0     
*/

class Utilities 
{
	//////////////////////////
	// Numerical constants
	//////////////////////////
	
	/** Number close the smallest representable number */	
	public static final double FPMIN   = 1.0e-320;   
 
	/** A tiny number */
    public static final double TOL = 1e-06;
	
 
	//////////////////////////
	//List of non-valid values
	//////////////////////////
	
	/** Invalid double value */	
	public static final double INVALID_D = Double.MAX_VALUE; 
	/** Invalid integer value */		
	public static final int    INVALID_I = Integer.MAX_VALUE;	
	/** Invalid integer value (string format) */
	public static final String INVALID_S = String.valueOf(Integer.MAX_VALUE); 
	
	//////////////////////////
	//Lists of codes
	//////////////////////////
	
	/** Code for unknown gender or affection */
	public static final int MISSING	= 0;	
	/** Code for gender (male) */
	public static final int MALE 	= 1;	
	/** Code for gender (female) */
	public static final int FEMALE 	= 2;	
	/** Code for unaffected individual */
	public static final int UNAFFECTED = 1; 
	/** Code for affected individual */
	public static final int AFFECTED   = 2;
	
	/** Code for analysis mode "pedcheck" */
	public static final int MODE_PEDCHECK 		= 0; 
	/** Code for analysis mode "heritability" */	
	public static final int MODE_HERITABILITY 	= 1;
	/** Code for analysis mode "association" */	
	public static final int MODE_ASSOCIATION 	= 2;	
	
	//////////////////////////
	//List of non-valid values
	//////////////////////////
	
	/** Strings used for representing missing value (array version)*/
	private static final String[] missingValues_a = {"X", "x", "NA", "na", "NaN", "NAN", "nan"}; 
	/** Strings used for representing missing value (Vector version)*/
	public static final Vector<String> MISSING_VALUES = new Vector<String>(Arrays.asList(missingValues_a));
	
	/** Strings used for representing the valid chromosomes (array version)*/
	private static final String[] validChr_a = {"1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22"}; 
	/** Strings used for representing valid chromosomes (Vector version)*/
	public static final Vector<Object> VALID_CHR = new Vector<Object>(Arrays.asList(validChr_a));
	
	
	/**
		Returns whether the code is valid.
		
		\see PedParser
		
		@param code the code to check
		@param validCodes the list of valid code
		@return true if the code is valid, false otherwise
		@see com.github.alesssia.poppante.Family
		@see com.github.alesssia.poppante.PedParser
		@see com.github.alesssia.poppante.Person
	*/
	public static boolean isValid(Object code, Vector<Object> validCodes)
	{
		return validCodes.indexOf(code) != -1;
	}
	
	
	/**
		Returns whether the value is missing.
		
		@param value the value to check
		@return true if the code is missing, false otherwise
		@see com.github.alesssia.poppante.MyFileReader
		@see com.github.alesssia.poppante.PedParser
	*/
	public static boolean isMissing(String value)
	{
		return MISSING_VALUES.indexOf(value) != -1;
	}

	/**
		Evaluates the inner product between two vectors 
		of double.
		
		@precondition the two vectors must have the same size.		
		@param v1 the first vector
		@param v2 the second vector
		
		@return the inner product
	*/
	public static double innerProduct(double[] v1, double[] v2)
	{
		assert (v1.length == v2.length) : "Internal error. Impossible to evaluate inner product.";
		
		double innerProduct = 0.0;
		for (int i = 0; i < v1.length; i++)
			innerProduct += v1[i] * v2[i];
		
		return innerProduct;
	}
	
	
	/** 
		Creates a r by c matrix of double and initialises 
		it with the given value.
		
		@param r the number of rows
		@param c the number of columns
		@param v the initial value
		@return the new matrix
	*/
	public static double[][] set(int r, int c, double v) 
	{
		double[][] m = new double[r][c];
		for(int i=0; i<r; i++)
			for(int j=0; j<c; j++)
				m[i][j] = v;
		return m;
	}
	
	/** 
		Creates a vector of double of length n and initialises 
		it with a given value.
		
		@param n the length of the vector
		@param v the initial value
		@return the new vector
	*/
	public static double[] set(int n, double v)
	{
		double[] r = new double[n];
		for (int i=0; i<n; i++)
			r[i] = v;
		return r;
	}
	
	/** 
		Creates a vector of integer of length n and initialises 
		it with the given value.
		
		@param n the length of the vector
		@param v the initial value
		@return the new vector
	*/
	public static int[] set(int n, int v)
	{
		int[] r = new int[n];
		for (int i=0; i<n; i++)
			r[i] = v;
		return r;
	}
	
	/** 
		Extracts from a vector of Double the values in the 
		positions listed by a vector of Integer
		
		@param v the vector of Double
		@param indices the indices of the values to extract
		@return the extracted values
	*/
	public static double[] getValuesfromVector(Vector<Double> v, Vector<Integer> indices)
	{
		double[] result = new double[indices.size()];
		
		for(int i=0; i< indices.size(); i++)
			result[i] = v.get(indices.get(i));
		
		return result;
	}

	/**
		Given a vector of array of double transform it in a
		Matrix Object. 
	
		@param v a vector of array of double
		@return a Matrix Object
	*/
	public static Matrix getMatrix(Vector<double[]> v)
	{
		double[][] tmp = new double[v.size()][v.get(0).length];
		for(int i=0; i<v.size(); i++)
			tmp[i] = v.get(i);

		return new Matrix(tmp);
	}
	
	/**
		Returns a string codifying an array of double. Values are 
		represented using the provided decimal format.
		
		@param array the array
		@param df the decimal format
		@param ef the esponential format
		@param threshould the threshould under which use the esponential format
		@return a string representing the array
	*/
	public static String array2string(double[] array, DecimalFormat df, DecimalFormat ef, double threshould)
	{
		String s = "(";
		for (int i=0; i<array.length; i++)
			
			s += " " + (array[i] < threshould ? ef.format(array[i]) : df.format(array[i]));

		s += " )";

		return s;
	}
	

	/**
		Evaluates the mean of a vector of Double.
	
		@param v the vector
		@return the mean
	*/
	public static double mean(Vector<Double> v)
	{
	    double s = 0.0;
            for (Double v1 : v) 
                s += v1;
    
	    return s/v.size();
	}

	/**
		Evaluates the variance of a vector of Double.
		
		An unbiased estimates can be calculated.
	
		@param v the vector
		@param correct whether given the unbiased estimate
		@return the variance
	*/
	public static double variance(Vector<Double> v, boolean correct)
	{
	    double mean = mean(v);
	    double s = 0.0;
            for (Double v1 : v) 
            {
                double a = v1;
                s += (mean-a)*(mean-a);
            }
    
	    if (correct)
	        return  s/(v.size()-1);

	    return s/v.size();
	}
	
	/**
		Transforms a vector of Double in an array of 
		double.
	
		@param v the vector of Double
		@return the array of double
	*/

	public static double[] vector2double(Vector<Double> v)
	{
	    double[] d = new double[v.size()];
	    for(int i=0; i<v.size(); i++)
	        d[i] = v.get(i);
	    return d;
	}

	/**
		Shuffles an array of double.
		
		@param ar the array to shuffle
		@return the shuffled array
	*/
	public static double[] shuffle(double[] ar)
	{
		Random rnd = new Random();
		for (int i = ar.length - 1; i > 0; i--)
		{
			int index = rnd.nextInt(i + 1);
			// Simple swap
			double a = ar[index];
			ar[index] = ar[i];
			ar[i] = a;
		}
         
		return ar;
	}
	
	/**
		Returns wheher a string codifies a double
		
		@param str the string to check
		@return whether the string codifies a double
	*/
	
	public static boolean isDouble(String str)
	{
	  return str.matches("-?\\d+(\\.\\d+)?");  //match a number with optional '-' and decimal.
	}
	
}
