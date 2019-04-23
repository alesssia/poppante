/*
 * 	 Result.java
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

import java.text.*;

/**
	Represents the result of an Association or 
	Heritability test

	@author      Alessia Visconti <alessia.visconti@kcl.ac.uk>
	@version     1.0     
*/

class Result 
{
	/**
		Value from which to swithch to decimal format to scientific
		notation */
	private static final double SCIENTIFIC_THRESHOLD = 0.001;
	/** Header of the result file.
		
		According to the set options more columns can be add.*/
	private static final String header = "Nobs\tLnlk_Null\tLnlk_Full\tdf_Null\tdf_Full\tchi^2\tpvalue\tadj_pvalue"; 
	
	/** Decimal format.
		
		It formats a value to have 2 decimal position */
	private static final DecimalFormat df2 = new DecimalFormat("####0.00");	
	/** Decimal format.
		
		It formats a value to have 3 decimal position */
	private static final DecimalFormat df3 = new DecimalFormat("0.000"); 	
	/** Decimal format.
		
		It formats a value to have 4 decimal position */
	private static final DecimalFormat df4 = new DecimalFormat("0.0000");	
	/** Decimal format.
		
		It formats a value to be in scientific notation */
	private static final DecimalFormat dfe = new DecimalFormat("0.00E0");	
	/** Decimal format.
		
		It formats a value to be in percentage notation */
    private static final DecimalFormat dfp = new DecimalFormat("##0.00%");	
	
	/** Phenotype name.
		Not set nor print if analysis mode is "heritability". */			
	private String phenotype;	
	/** Name of the methylation site. */
	private String marker;	
	/** Methylation site  chromosome. */
	private String chr;	
	/** Methylation site  position. */	
	private long position;	

       /** Number of observations. */
	private int nobs;
	/** Log likelihood of the Null model. */
	private double lnlkNull;	
	/** Log likelihood of the Full model. */ 
	private double lnlkFull;	
	/** Degree of freedom of the Null model. */ 
	private int dfNull;     
	/** Degree of freedom of the Full model. */ 	
	private int dfFull;    
	/** Chi-square of the test. */  	
	private double chiq;    
	/** p-value of the test. */ 	
	private double pvalue;  
	/** empirical p-value of the test.
		
		It is evaluated only if the experiment-wise error rate
		is specified as option.	*/ 	
	private double epvalue;  
	/** adjusted p-value of the test. */ 	
	private double adjpvalue;  
	/** adjusted empirical p-value of the test.
		
		It is evaluated only if the experiment-wise error rate
		is specified as option. */
	private double adjepvalue;  
		
	/** Percentage of families showing a positive contribution  
		
		Not set nor print if relc parameter is bot set set 
		nor if no family are supplied, that is we are working with 
		an external kinship matrix and unrelated individuals. 
		@see com.github.alesssia.poppante.Constants
		@see com.github.alesssia.poppante.VC
		*/
	private double posF;	
	/** Gini coefficient to assess family contribution to the chi-square  
		
		Not set nor print if no family are supplied, that is 
		we are working with an external kinship matrix and
		unrelated individuals. */	
	private double giniC;   	
	/** Variance of the Null model */
	private double[] varNull;	
	/** Variance of the Full model */
	private double[] varFull;	
	/** Heritability value 
		
		Not set nor print if analysis mode is "association".*/
	private double heritability;
	
	/** Coefficient of the regression model for the
		methylation site.
		
		Not set nor print if analysis mode is "heritability". */
		private double beta;	
	/** Standard error for the coefficient of the regression 
		model for the methylation site.
	
		Not set nor print if analysis mode is "heritability". */	
	private double se;	
		
	/** Variance explained by the	methylation site.
		
		Not set nor print if analysis mode is "heritability". */		
	private double ve; 		
	/** Error message.
		
		Set if an exception has been raised during the 
		evaluation of the test. */
	private String emessage;
	
	
	/**
		Constructor. 
		
		Crates an empty object.
	*/	
	public Result()
	{
		marker = null;
	}
	
	/**
		Constructor. 
		
		Used when the test has been successful.
		
                @param n the number of observation used
		@param lnlkN log likelihood of the Null model
		@param lnlkF log likelihood of the Full model
		@param dfN degree of freedom of the Null model
		@param dfF degree of freedom of the Full model
		@param chi2 chi-square of the test
		@param pval p-value of the test
		@param epv empirical p-value
		@param pF percentage of families showing a positive contribution
		@param gini Gini coefficient
		@param vNull variance of the Null model
		@param vFull variance of the Full model
		@param h the heritability value
		@param b the beta
		@parma s the standard error
		@param v the variance explained
	*/
	public Result(int n, double lnlkN, double lnlkF, int dfN, int dfF, double chi2, double pval, double epv, double pF, double gini, double[] vNull, double[] vFull, double h, double b, double s, double v)
	{
		phenotype = null;
		marker = null;
		chr = null;
		position = -1;
		
		nobs = n;
		lnlkNull = lnlkN;
		lnlkFull = lnlkF;
		dfNull = dfN;
		dfFull = dfF;
		chiq = chi2;
		pvalue = pval;
		epvalue = epv;
		posF = pF;
		giniC = gini;
		heritability = h;
		beta = b;
		se=s;
		ve = v;
		
		varNull = vNull;
		varFull =vFull;
		
		emessage = null;
	}
	
	/**
		Constructor. 
		
		Used when the test has not been successful and an error has
		been raised.
		
		@param emsg error message
	*/
	public Result(String emsg)
	{
		emessage = emsg;
	}
	
	/**
		Sets the phenotype and methylation site information.
		
		@param p phenotype name
		@param m methylation site name
		@param c methylation site chromosome
		@param pos methylation site position
	*/
	public void setInfo(String p, String m, String c, long pos)
	{
		phenotype = p;
		marker = m;
		chr = c;
		position = pos;
	}
	
	/**
		Returns the pvalue
		
		@return the pvalue
	*/
	public double pvalue()
	{
		return pvalue;
	}
	
	/**
		Returns the empirical p-value
		
		@return the empirical p-value
	*/
	public double epvalue()
	{
		return epvalue;
	}
	
	/**
		Returns the chromosome
		
		@return the chromosome
	*/
	public String chr()
	{
		return chr;
	}
	
	/**
		Returns the position
		
		@return the position
	*/
	public long position()
	{
		return position;
	}
	
	/**
		Returns the marker name
		
		@return the marker name
	*/
	public String marker()
	{
		return marker;
	}
	
	/**
		Sets the adjusted p-value
		
		@param adjpv the adjusted p-value
	*/
	public void setAdjPvalue(double adjpv)
	{
		adjpvalue = adjpv;
	}
	
	
	/**
		Sets the adjusted empirical p-value
		
		@param adjpv the adjusted empirical p-value
	*/
	public void setAdjePvalue(double adjpv)
	{
		adjepvalue = adjpv;
	}
	
	/**
		Updates the header.
		
		The standard header is updated according to the set
		options and returned.
		
		@param methFields the number of fields used to describe a methylation site
		@return the updated header
	*/
	public static String header(int methFields)
	{
		String r = header; //standard header
		
		if (methFields == 1)
			r = "Predictor\t" + r;
		else if (methFields == 3)
			r = "Predictor\tChr\tPosition\t" + r;
		
		if (Constants.alpha != Utilities.INVALID_D)
			r += "\tepvalue\tadj_epvalue";
		
		if (Constants.mode == Utilities.MODE_HERITABILITY)	//associations starts with phenotype
			r += "\tHeritability";
		
		if (Constants.mode == Utilities.MODE_ASSOCIATION)	//associations starts with phenotype
			r = "Response\t" + r + "\tbeta\tse\tVariance.Explained";
		
		if (Constants.relc != Utilities.INVALID_D)	 //extra statistics
			r += "\tPosF\tGiniC";
		
		if (Constants.variance) //variance
			r += "\tvar_Null\tvar_Full";
			
		return r;
	}
	
	
	/** 
		Returns the result.
		
		The result is converted to a String and the field printed 
		vary according to the set options.
		
		@precondition the result should be initialised
		
		@return the result
	*/	
        @Override
	public String toString()
	{
		assert marker != null : "Internal error: not well-formed result";
		
		String result = "";
		
		if (Constants.mode == Utilities.MODE_ASSOCIATION)
			result += phenotype + "\t" ;
		
		result += marker + "\t" ;
		
		if (chr != null)
			result += chr + "\t" + position + "\t";
		
		if (emessage != null) //exception thrown or other error message
		{
			result += emessage;
			return result; 
		}
			
		result += nobs + "\t" + df2.format(lnlkNull) + "\t" + df2.format(lnlkFull) + "\t" + dfNull + "\t" + dfFull + "\t";
		result += df2.format(chiq) + "\t" + (pvalue < SCIENTIFIC_THRESHOLD ? dfe.format(pvalue) : df4.format(pvalue)); 
		result += "\t" + (adjpvalue < SCIENTIFIC_THRESHOLD ? dfe.format(adjpvalue) : df4.format(adjpvalue)); 
		
		if (Constants.alpha != Utilities.INVALID_D)
		{
			result += "\t" + (epvalue < SCIENTIFIC_THRESHOLD ? dfe.format(epvalue) : df4.format(epvalue));
			result += "\t" + (adjepvalue < SCIENTIFIC_THRESHOLD ? dfe.format(adjepvalue) : df4.format(adjepvalue));
		}
		
		if 	(Constants.mode == Utilities.MODE_HERITABILITY)
			result += "\t" + df4.format(heritability);

		if 	(Constants.mode == Utilities.MODE_ASSOCIATION)
		{
			result += "\t" + (Math.abs(beta) < SCIENTIFIC_THRESHOLD ? dfe.format(beta) : df4.format(beta));
			result += "\t" + (Math.abs(se) < SCIENTIFIC_THRESHOLD ? dfe.format(se) : df4.format(se));
			result += "\t" + dfp.format(ve);
		}

		if (Constants.relc != Utilities.INVALID_D)
			result += "\t" + dfp.format(posF) + "\t" + df2.format(giniC);

		if (Constants.variance)
		{
			String vNull = Utilities.array2string(varNull, df4, dfe, SCIENTIFIC_THRESHOLD);
			String vFull = Utilities.array2string(varFull, df4, dfe, SCIENTIFIC_THRESHOLD);
		
			result += "\t" + vNull +  "\t" + vFull;
		}

		return result;
	}
		
}
	
