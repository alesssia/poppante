/*
 * 	 ResultTable.java
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
import com.github.alesssia.statutils.*;
import manhattanplotter.*;

/**
	Represents the result table.

	It is used to adjust the p-values (theoretical and empirical)
	and to plot a summary of the results (Q-Q plot, Manhattan plot).

	To plot the data it uses the manhattanplotter java package 
	(see https://github.com/CarlosMorcilloSuarez/ManhattanPlotter for details)
	
	@author      Alessia Visconti <alessia.visconti@kcl.ac.uk>
	@version     1.0     

	@see com.github.alesssia.arrayutils.FDR
	@see com.github.alesssia.arrayutils.Result
*/
class ResultTable extends Experiment
{
	/** Results */
    private final Vector<Result>  results;		
	/** FDR engine */
    private FDR fdr;		
	/** List of p-values for FDR */
	private Double[] pvalues;
    	
	/**
		Constructor. 
		
		Initialises the data structures.
		It fills the information needed by the ManhattanPlotter, and 
		extracts the pvalues for correction (performed by means of the 
		Benjamini–Hochberg procedure).
	
		It works both with theoretical and empirical p-valus 

		@param res the list of independent results
		@param empirical whether the empirical p-values should be considered instead of the theoretical
	*/
    public ResultTable(Vector<Result> res, boolean empirical) 
    {
        super(MANHATTAN_EXPERIMENT);
		this.results = res;
		
		pvalues = new Double[results.size()];
		for (int i = 0; i < results.size(); i++)
		{
			Result r = results.get(i);
			
			manhattanplotter.Test mytest = convertExperiment(r, empirical);
			if (mytest != null)
				testsArray.add(mytest);
			
			if (empirical)
				pvalues[i] = r.epvalue();
			else
				pvalues[i] = r.pvalue();
		}
		
		terminateExperiment();
        fdr = new BH();
    }
	
	/**
		Converts a Result in a object that can be used by the 
		plotting methods.
		
		@param r a result
		@param empirical whether the empirical p-values should be considered instead of the theoretical
		@return the converted result
	*/
	private manhattanplotter.Test convertExperiment(Result r, boolean empirical)
	{
        manhattanplotter.Test mytest = new manhattanplotter.Test();

		//chr info
		String chr = r.chr();	
		if (!Utilities.isValid(chr, Utilities.VALID_CHR))
			return null;
        
		mytest.chromosome = Integer.valueOf(chr);
		
		//position info
		mytest.position   = r.position();
		mytest.coordinate = genomeInfo.chromosome[mytest.chromosome].initCoordinate + mytest.position;

		if (mytest.position > genomeInfo.chromosome[mytest.chromosome].length)
		{
            genomeInfoOK = false;
            genomeInfo.chromosome[mytest.chromosome].length = mytest.position;
        }
		
		//mrk info
		mytest.markerName = r.marker();
       
		if (empirical)
			mytest.pValue = r.epvalue();
		else
			mytest.pValue = r.pvalue();
		
		mytest.logPValue = -Math.log10(mytest.pValue);
        if (mytest.logPValue > maxLogPValue) 
            maxLogPValue = mytest.logPValue;

      
		return mytest;
	}

	/**
		Terminates the creation of the Experiment object, 
		used by the plotting methods.	
	*/
	private void terminateExperiment()
	{
	    // converts arrayList to array
	    super.tests = new manhattanplotter.Test[testsArray.size()];

	    for (int i = 0; i < tests.length; i++)
	        tests[i] = testsArray.get(i);
	    //testsArray = new ArrayList<manhattanplotter.Test>();

	    // Recreates genomeInfo and recalculates tests coordinates if necesary
	    if(!genomeInfoOK)
		{
	        genomeInfo.rebuild();
	        for(int i = 0; i < tests.length; i++)
	            tests[i].coordinate = genomeInfo.chromosome[tests[i].chromosome].initCoordinate + tests[i].position;
	    }
	}
	
	/**
		Performs the FDR correction and set the adjusted p-values
		(theoretical or empirical)
		
		@param empirical whether the empirical p-values should be considered instead of the theoretical
	*/
    public void adjustPvalues(boolean empirical)
    {
		//theoretical p-values
		fdr = new BH(pvalues);
        fdr.adjust();
		
		if (empirical)
		{
			for (int i = 0; i < results.size(); i++)
				results.get(i).setAdjePvalue(fdr.adjPvalue(i));
		}
		else
		{
			for (int i = 0; i < results.size(); i++)
            	results.get(i).setAdjPvalue(fdr.adjPvalue(i));
		}
    }
	
	/**
		Creates and saves the Q-Q plot. If no valid data are available 
		an error message is returned instead.
		
		@param empirical whether the empirical p-values should be considered instead of the theoretical
		@return the filepath of the created Q-Q plot or an error message
	*/
    public String qqplot(boolean empirical)
    {
		//the marker map or the tests do not contain any entry
		//with a valid position
		String error = validData2Plot(empirical, "Q-Q plot");
		if (error.length() != 0)
			return error;
		
		QQPlot qq = new QQPlot();
		qq.setExperiment(this);
		qq.drawPlot();
		
		String outfile = outputfile("QQ-plot", empirical);
		qq.savePlot(outfile);	
		return outfile + ".png";
	}
	
	/**
		Creates and saves the Manhattan plot. If no valid data are available 
		an error message is returned instead.
		
		@param empirical whether the empirical p-values should be considered instead of the theoretical
		@return the filepath of the created Manhattan plot or an error message
	*/
    public String manhattan(boolean empirical)
    {
		//the marker map or the tests do not contain any entry
		//with a valid position
		String error = validData2Plot(empirical, "Manhattan plot");
		if (error.length() != 0)
			return error;
		
		ManhattanPlot mp = new ManhattanPlot(this);
		mp.setMode(ManhattanPlot.ALL_DATA, 0, 0, 0);
		mp.setDotSize(7);
		mp.drawPlot();
		
		String outfile = outputfile("Manhattan-plot", empirical);
		mp.savePlot(outfile);
		return outfile + ".png";
	}
	
	/**
		Creates the ouitput filename for the plot
		
		@param plot which plot is going to be saved
		@param empirical whether the empirical p-values should be considered instead of the theoretical
		@return the output file path
	 */
	private String outputfile(String plot, boolean empirical)
	{
		String out = "";
		
		if (Constants.output == null)
			out += "./";
		else
			out += Constants.output + "_";
		
		out += plot;
		if (empirical)
			out += "_epvalue";
		
		return out;
	}
	
	/**
		Verifies whether there are valid results to plot, that
		is if there is at least one test in a an autosomal
		chromosome.
		
		@param empirical whether the empirical p-values should be considered instead of the theoretical
		@param plot plot to be created (Q-Q or manhattan plot)
		@return an empty string if there is at least a valid results, a string with a error message otherwise
	*/
	private String validData2Plot(boolean empirical, String plot)
	{
		//the marker map or the tests do not contain any entry
		//with a valid position
		String error = "";
		if (testsArray.size() < 1)
		{
			error += "no " + plot + " created for ";
			if (empirical) 
				error += "empirical ";
			error += "p-values. Please check your methylation data file";
		}
		
		return error;
	}
	
	/**
		Returns the vector of Results where the adjusted p-values 
		have been set.
		
		@return the  the list of independent results
 	 */
	public Vector<Result> results()
	{
		return results;
	}

}
