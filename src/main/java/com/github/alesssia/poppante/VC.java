/*
 * 	 VC.java
 *
 *   PopPAnTE is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   PopPAnTE is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with PopPAnTE.  If not, see <http://www.gnu.org/licenses/>.
 *
 *   For any bugs or problems found, please contact us at
 *   mario.falchi@kcl.ac.uk  
 *	 hariklia.eleftherohorinou06@imperial.ac.uk
 *   alessia.visconti@kcl.ac.uk
 */

package com.github.alesssia.poppante;

import mathutils.*;
import java.util.*;
import com.github.alesssia.probutils.*;

/** 
	Describes a Variant Component (VC) Models.

	Each VC describes a test and is composed by two system of equations, 
	one describing the so called "full model" and one describing a "null
	model" where the tested variable is removed. When the analysis mode is 
	set to "heritability" the null model does not include the effect of 
	individuals’ relationships.	When the analysis mode is set to "association" 
	the null model does not include the value of the tested methylation site.

	References for the implementation can be found at:
		Press, William H., et al. Numerical recipes in C. 
		Vol. 2. Cambridge: Cambridge university press, 1996.

	@author      Alessia Visconti <alessia.visconti@kcl.ac.uk>
	@author      Mario Falchi     <mario.falchi@kcl.ac.uk>
	@author      Hariklia Eleftherohorinou <hariklia.eleftherohorinou06@imperial.ac.uk>
	@version     1.0     
*/

class VC extends Permutable
{
	/** This is where the data are stored */
	private DataManager datamanager;
	/** List of analysable families along with the
		number of members that can be analysed */
	private Analysable analysable;
	/** VC dimension of the Full model */
	public static final int VC_FULL = 2; 	
	/** VC dimension of the Null model */
	public static final int VC_NULL = 2; 	/** VC dimension parameter */
	
	/** Number of Subjects used in the model */
	private int numSubjects;		 // used for the degree of freedom 
	/** Number of equations currently loaded */
	private int counter;			/** number of models currently loaded */
	/** The position of the phenotype in the person's phenotype list */                            	
	private final int phn;				
	/** The position of the methylation site in the person's methylation site list */
	private final int site;				
	/** Sites used in the region-based tests */
	private Vector<Integer> sites;	

	/** Linear equations for Null models */
	private final NormalSet nullSet;	
	/** Linear equations for Full models */
	private final NormalSet fullSet;	
	/** Number of components in the Null models */
	private final int linearNull;
	/** Number of components in the Full models */			
	private final int linearFull;	
	
	/** Percentage of families showing a positive contribution.
		
		It is evaluated only if the relc parameter is set and if 
		family are supplied.
		@see com.github.alesssia.poppante.Constants  
	*/                                 
	private double posF;	
	/** Gini coefficient to assess family contribution to the chi-square.
		
		It is evaluated only if the relc parameter is set and if 
		family are supplied.
		@see com.github.alesssia.poppante.Constants  
	*/         		
	private double giniC;			
	/** Predictor values. 
		
		They are used for variance explained. */
	private final Vector<Double> meths;	
	/** Outcome values.
		
		They are used for variance explained.  */
	private final Vector<Double> ys; 		
	/** Log likelihood of the null model.
		
		It is evaluated only once for tests having 
		the same missingness pattern ant it is
		also used multiple times during the evaluation
		of permutation tests  */
	private double statisticNull; 	
	/** Values for the null model.
		
		If not evalated yet it equal to null. */
	private NullModel nullModel;
	/** Whether then null model should be set */
	private boolean setNull;

	/**
		Constructor. 
		
		@param p position of the phenotype in the person's phenotypes list
		@param s position of the methylation site in the person's methylation sites list
		@param families the list of families to analyse
		@param dm the data manager object
	*/
    public VC(int p, int s, Analysable families, DataManager dm) 
	{
		phn = p;
		site = s;
		sites = null;
		
		datamanager = dm;
		analysable = families;
			
		//Initialises linear models' size and model themselves
		linearNull = datamanager.numCovar() + 1; 
		if (Constants.mode == Utilities.MODE_ASSOCIATION)
		{
			//Gets the null model
			nullModel = datamanager.getNull(phn);
			
			//tries to undestand whether the null model can be used or not
			//If the missingness pattern are the same for the variables 
			//then I can use it and there is not need to set the null model
			//However, this can be the first time I analyse this phenotype, 
			//that is still not set, so I need to evaluate it (and set it)
			//anyway
			setNull = !nullModel.isSet() || !datamanager.missingnessMethPattern(site).equals(nullModel.missingnessPattern());
			
			linearFull = datamanager.numCovar() + 2;  
			if (setNull)
				nullSet = new NormalSet(analysable.numFamilies(), linearNull, VC_NULL);
			else
				nullSet = null;
		}
		else
		{
			setNull = true;
			linearFull = datamanager.numCovar() + 1; //the linear model has not fixed effect
			nullSet = new NormalSet(analysable.numFamilies(), linearNull, VC_NULL-1 ); //vg is not loaded
		}
		fullSet = new NormalSet(analysable.numFamilies(), linearFull, VC_FULL);
		
		//initialises an handful of data structures
		counter = -1;
		numSubjects = 0;	
	
		posF = Utilities.INVALID_D;
		giniC = Utilities.INVALID_D;
		
		meths = new Vector<>();
		ys = new Vector<>();
	}
	
	
	/**
		Sets the list of methylation site values to use in
		the region-base tests. 
		
		The values are represented by their position in the 
		methylation sites list.
		
		@param sites the position of the methylation site values to use.
	*/
	public void setSites(Vector<Integer> sites)
	{
		this.sites = sites;
	}
	
	
	/**
		Returns the set of linear equations for the Null model(s)
		
		@return the the set of linear equations for the Null model(s)
	*/
	public NormalSet nullSet()
	{
		return nullSet;
	}
	
	/**
		Returns the set of linear equations for the Full model(s)
		
		@return the set of linear equations for the Full model(s)
	*/
	public NormalSet fullSet()
	{
		return fullSet;
	}
		

	/**
		Initialises the VC information for each equation (family). 
		
		Each model represents a family. When unrelated individuals
		are used there is only a family.
	*/
	public void fill()
	{
		for (int i=0; i<analysable.numFamilies(); i++)
		{
			prepare(analysable.analysableFamiles(i).numAnalysable());
			AnalysableFamily family = analysable.analysableFamiles(i);
			
			fillKinship(family); 			//setVarComponents
			
			if (datamanager.numCovar() != 0)
				fillCovariates(family); 		//setLinearModel
		
			if (Constants.mode == Utilities.MODE_ASSOCIATION)
			{	
				fillPhenotScore(family); //setScores			
				fillMethSites(family);  	//setLinearModel 
			}
			else
	 			fillMethScore(family); 	//setScores	
			
		}	
	}

	/**
		Crates the equation structure.
		
		It initialises the equation parameter according to the number of 
		analysable individual in the given family.
		
		An individual is analysable if she has a valid value for 
		the methylation site m (when the analysis mode is 
		"heritability") or both the phenotype p and methylation 
		site m (when the analysis mode is "association").
		
		@param size the number of analysable family members
		@see com.github.alesssia.poppante.Normal
	*/
	private void prepare(int size) 
	{
		counter++;		//it's the position of the family in the sets of model

		(fullSet.sets(counter)).prepare(size);
		
		if (setNull)
			(nullSet.sets(counter)).prepare(size);
		
		numSubjects += size;
	}
		
	/**
		Fills the linear models (predictors) with the covariate values.
		
		@param f the family to be modelled
		@see com.github.alesssia.poppante.Normal
	*/	
	private void fillCovariates(AnalysableFamily f)
	{
	    int l = 0;
		Family family = datamanager.families().get(f.famID());
		for (int i=0; i<family.numMembers(); i++)
		{
			int position = f.startPos()+i;
			if (analysable.isAnalysable(position))
			{
				double[] covariates = new double[datamanager.numCovar()];
				for (int c = 0; c < covariates.length; c++) 
				{
					(fullSet.sets(counter)).setLinearModel(l, (linearFull-datamanager.numCovar()+c), datamanager.covariates()[c][position]);
					
					if (setNull)
						(nullSet.sets(counter)).setLinearModel(l, (linearNull-datamanager.numCovar()+c), datamanager.covariates()[c][position]);
				}
				l++;
			}
		}
	}

	/**
		Fills the variant component with the kinship values.
		
		@param f the family to be modelled
		@see com.github.alesssia.poppante.Normal
	*/
	private void fillKinship(AnalysableFamily f) 
	{
		Family family = datamanager.families().get(f.famID());
		Kinship kinship = family.kinship();
		
		int ik = 0;		
		for (int i=0; i<family.numMembers(); i++)
		{
			if (analysable.isAnalysable(f.startPos()+i))
			{
				int jk = ik;
				for (int j=i; j<family.numMembers(); j++)
				{
					if (analysable.isAnalysable(f.startPos()+j))
					{
						double value = kinship.getValue(i, j);

						//vg component should be always filled when association is calculated
						//and only for the full model when heritability is calculated
						if (Constants.mode == Utilities.MODE_ASSOCIATION && setNull)
							(nullSet.sets(counter)).setVarComponents(1, ik, jk, value);
						
						(fullSet.sets(counter)).setVarComponents(1, ik, jk, value);
						
						jk++;
					}	
				}
				ik++;				
			}
		}
	}


	/**
		Fills the model outcomes with the phenotype values. 
		
		This function is used only when the analysis mode is set 
		to "association".
		
		@param f the family to be modelled
		@see com.github.alesssia.poppante.Normal
	*/		
	private void fillPhenotScore(AnalysableFamily f) 
	{
		assert Constants.mode == Utilities.MODE_ASSOCIATION : "Internal error: VC should not use the response values as dependent variable when heritability testing is performed.";
     
		int l = 0;
		Family family = datamanager.families().get(f.famID());
		for (int i=0; i<family.numMembers(); i++)
		{
			int position = f.startPos()+i;
			if (analysable.isAnalysable(position))
			{
				double phenotype = datamanager.phenotypes()[phn][position];
				ys.add(phenotype); //FIXME: why vector?
				
				if (setNull)
					(nullSet.sets(counter)).setScores(l, phenotype);
				
				(fullSet.sets(counter)).setScores(l, phenotype);
				l++;
			}
		}
	}
	
	/**
		Fills the model scores with the methylation values. 
		
		This function is used only when the analysis mode is set 
		to "heritability".
		
		@param f the family to be modelled
		@see com.github.alesssia.poppante.Normal
	*/	
	private void fillMethScore(AnalysableFamily f)
	{
		
		assert Constants.mode == Utilities.MODE_HERITABILITY : "Internal error: VC should not use the predictor values as dependent variable when association testing is performed.";
			
		if (sites == null || sites.size() == 1)
			fillMethScoreStandard(f);
		else
			fillMethScorePC(f);
	}
	
	/**
		Fills the model outcomes with the methylation values. 
		
		This function is used only when the analysis mode is set 
		to "heritability" and we are testing a single methylation
		site.
				
		@param f the family to be modelled
		@see com.github.alesssia.poppante.Normal
	*/	
	private void fillMethScoreStandard(AnalysableFamily f) 
	{     		
		int l = 0;
		Family family = datamanager.families().get(f.famID());
		for (int i=0; i<family.numMembers(); i++)
		{
			int position = f.startPos()+i;
			if (analysable.isAnalysable(position))
			{
				double methylation = datamanager.methylations()[site][position];
				(nullSet.sets(counter)).setScores(l, methylation); 
				(fullSet.sets(counter)).setScores(l, methylation); 
				l++;
			}
		}
	}

	/**
		Fills the model outcomes with the methylation values. 
		
		This function is used only when the analysis mode is set 
		to "heritability" and we are performing a region-based test.
		It fills the outcomes with the PC1 of the methylation
		values included in the given window.

		@param f the family to be modelled
		@see com.github.alesssia.poppante.Marker
		@see com.github.alesssia.poppante.MarkerRegion	
		@see com.github.alesssia.poppante.Normal
	*/
	private void fillMethScorePC(AnalysableFamily f) //FIXME
	{
		//I can use directly the PC values because I have a value per person,
		//and the invalid people have already been marked as such (in fact the
		//PC has valid values for them)
		MarkerRegion markerRegion = new MarkerRegion(sites, datamanager.methylations());
		double[] pc = markerRegion.summarise();
		
		int l = 0;
		Family family = datamanager.families().get(f.famID());
		for (int i=0; i<family.numMembers(); i++)
		{
			int position = f.startPos()+i;
			if (analysable.isAnalysable(position))
			{
				double methylation = pc[position];
				(nullSet.sets(counter)).setScores(l, methylation);
				(fullSet.sets(counter)).setScores(l, methylation);
				l++;
			}
		}	
	}
	
	/**
		Fills the linear model with the methylation values. 
		
		This function is used only when the analysis mode is set 
		to "association".
		
		@param f the family to be modelled	
		@see com.github.alesssia.poppante.Normal
	*/
	private void fillMethSites(AnalysableFamily f) 
	{
		assert Constants.mode == Utilities.MODE_ASSOCIATION : "Internal error: VC should not use the predictor values as independent variable when heritability testing is performed.";
			
		if (sites == null || sites.size() == 1)
			fillMethSitesStandard(f);
		else
			fillMethSitesPC(f);
	}
	
	/**
		Fills the linear model with the methylation values. 
		
		This function is used only when the analysis mode is set 
		to "association" and we are testing a single methylation
		site.
				
		@param f the family to be modelled
		@see com.github.alesssia.poppante.Normal
	*/
	private void fillMethSitesStandard(AnalysableFamily f) 
	{
		int l = 0;
		Family family = datamanager.families().get(f.famID());
		for (int i=0; i<family.numMembers(); i++)
		{
			int position = f.startPos()+i;
			if (analysable.isAnalysable(position))
			{
				double methylation = datamanager.methylations()[site][position];
                meths.add(methylation); 
				(fullSet.sets(counter)).setLinearModel(l, (linearFull-datamanager.numCovar()-1), methylation); 
				l++;
			}	
		}
	}
	
	/**
		Fills the linear model with the methylation values. 
		
		This function is used only when the analysis mode is set 
		to "association" and we are performing a region-based test.
		It fills the outcomes with the PC1 of the methylation
		values included in the given window.

		@param f the family to be modelled
		@see com.github.alesssia.poppante.Marker
		@see com.github.alesssia.poppante.MarkerRegion	
		@see com.github.alesssia.poppante.Normal		
	*/	
	private void fillMethSitesPC(AnalysableFamily f) //FIXME
	{
		MarkerRegion markerRegion = new MarkerRegion(sites, datamanager.methylations());
		
		//I can use directly the PC values because I have a value per person,
		//and the invalid people have already been marked as such (in fact the
		//PC has valid values for them)
		double[] pc = markerRegion.summarise();
		
		int l = 0;
		Family family = datamanager.families().get(f.famID());
		for (int i=0; i<family.numMembers(); i++)
		{
			int position = f.startPos()+i;
			if (analysable.isAnalysable(position))
			{
				(fullSet.sets(counter)).setLinearModel(l, (linearFull-datamanager.numCovar()-1), pc[position]);		
				l++;
			}
		}
	}
	
	/**
		Solves the VC model.
		
		After finding the best model through a optimisation procedure,
		the tests is performed and its significance is evaluated by 
		comparing the likelihood of the full model and the likelihood of 
		a null model. In fact, twice the difference between the logarithms 
		of the full and the null likelihoods asymptotically follows a 
		chi-square distribution with 1 degree of freedom and can be used 
		to evaluate the difference in goodness-of-fit between the two models
		and to derive a p-value for the full model.
		
		If the relc parameter is set and if family are supplied the 
		Percentage of families showing a positive contribution and the
		Gini coefficient to assess family contribution to the chi-square
		are evaluated.
                 
		If the permutation option is set and the analysis mode is 
		"association" it evaluates an empirical p-value by means of
		a set of permutation tests.    
		
		@return the result of the test
		@see com.github.alesssia.poppante.Constants  
		@see com.github.alesssia.poppante.Normal
		@see com.github.alesssia.poppante.Result
	*/                                 	
	public Result solve() 
	{
		statisticNull = Utilities.INVALID_D;
		int dfNull = Utilities.INVALID_I;
		double[] variancesNull = null;
		double[] detailsNull = null;
		
		try 
		{
			//If not done beforehands the Null Model should be
			//evaluated
			if (setNull)
			{
				nullSet.disableConstant();
				nullSet.solve();
				nullSet.enableConstant();	
				
				//Permute must not take any parameter as required by the
				//adapte function. To speed up the computation I can,
				//however, avoid the evaluation of the "statisticNull"
				//and I am thus storing it in an istance variable (while the
				//"statisticFull" is stored in a local variable)
				//This value is also set globally, if this test has a 
				//common missingness pattern
				statisticNull = nullSet.evaluate();
				dfNull = numSubjects - nullSet.countParameters();
				variancesNull = nullSet.variances();
				detailsNull = nullSet.evaluateDetails();
				
				//now that I have calculate it, I store it (but only
				//if this is my missingness pattern), otherwise who cares?
				//Also, I need to store it only for association testing
				if (Constants.mode == Utilities.MODE_ASSOCIATION && datamanager.missingnessMethPattern(site).equals(nullModel.missingnessPattern()))
				{
					nullModel.setlogLikelihood(statisticNull);
					nullModel.setDf(dfNull);
					nullModel.setVariances(variancesNull.clone());
					nullModel.setDetails(detailsNull.clone());
				}
			}
			else
			{
				statisticNull = nullModel.logLikelihood();
				dfNull = nullModel.df();
				variancesNull = nullModel.variances();
				detailsNull = nullModel.details();
			}
			
			fullSet.solve();
			fullSet.disableConstant();
			fullSet.enableConstant();
		} 
		catch (Exception e) 
		{
			return new Result(e.getMessage());
		}
		
				
		double statisticFull = fullSet.evaluate();
		int dfFull = numSubjects - fullSet.countParameters();
		
                if (dfFull < 1 || dfNull < 1)
                    return new Result("Warning : not enough observations to estimate the model parameters");
                
		double chi2 = 2 * (statisticNull - statisticFull);
		if (chi2 < 0)
			chi2 = 0;
		double pvalue = ProbDist.gammq(0.5 * 1, 0.5 * chi2);

		if (Constants.relc != Utilities.INVALID_D && pvalue <= Constants.relc) 
			LKstats(fullSet.evaluateDetails(), detailsNull);

		double heritability = fullSet.variances(1)/(fullSet.variances(0) + fullSet.variances(1));
		
		//sets the variance explained
		double ve = Utilities.INVALID_D;
		double beta = Utilities.INVALID_D;
		double se = Utilities.INVALID_D;
		if (Constants.mode == Utilities.MODE_ASSOCIATION)
		{
			beta = fullSet.means(1);
			double var = Utilities.variance(meths, true);
			ve = beta * beta *  var / Utilities.variance(ys, true);
			se = Math.abs(beta/Math.sqrt(chi2)); 
		}
		//performs the permutation tests
		double epvalue = Utilities.INVALID_D;
		if (Constants.alpha != Utilities.INVALID_D)
		{
			Adaptive<VC> adaptive = new Adaptive<>(Constants.alpha, Constants.c, this);
			epvalue = adaptive.adapt(pvalue);
		}
		
		return new Result(numSubjects, statisticNull, statisticFull, dfNull, dfFull, chi2, pvalue, epvalue, posF, giniC, variancesNull, fullSet.variances(), heritability, beta, se, ve);
	}

	/**
		Performs a permutation test and returns the p-value
		obtained by the permuted model.
		
		If the analysis mode is "association" the evaluation of
		the null model is not needed, since the permuted values
		are the methylation sites in the predictor. The only 
		values of the null model needed to evaluate the new 
		p-value is its log-likelihood, that is given as parameter.
		If the analysis mode is "heritability" the permutation
		is not performed.

		@return the p-value of the permuted model
	*/
        @Override
	public double permute()
	{
		try
		{
			//permutes the data of each family
			//after this step the full models are changed forever
			for (int i=0; i<counter; i++)
			{
				double[] meth = (fullSet.sets(i)).predictor(1); //Meths			
				meth = Utilities.shuffle(meth);
				(fullSet.sets(i)).setPredictor(meth, 1);
			}
			
			//solves the new full model
			fullSet.disableConstant();
			fullSet.solve();
			fullSet.enableConstant();
			double statisticFull = fullSet.evaluate();
			
			double testValue = 2 * (statisticNull - statisticFull);
			testValue = (testValue < 0) ? 0 : testValue;

			return ProbDist.gammq(0.5 * 1, 0.5 * testValue);
		}
		catch (Exception e)
		{
			//doing nothing, is the system can't be solved I assume
			//that the p-value of the permutation is higher than that
			//of the real model
			return 1.1;
		}
	}


	/**
		Percentage of families showing a positive contribution and the
		Gini coefficient.
		
		@param famLKfull the -(log-likelihood) of the equations belonging to the Full model
		@param famLKnull the -(log-likelihood) of the equations belonging to the Null model
		@see com.github.alesssia.poppante.Normal
	*/	
	private void LKstats(double[] famLKfull, double[] famLKnull) 
	{
		double sum = 0;
		double gsum = 0;
		int posit = 0;
		double[] valuesAll = new double[famLKfull.length];

		for (int i = 0; i < famLKfull.length; i++) 
		{
			valuesAll[i] = famLKnull[i] - famLKfull[i];
			if (valuesAll[i] > 0)
				posit++;
		}
		
		double[] values = new double[posit];
		
		int idx = 0;
		for (int i = 0; i < famLKfull.length; i++)
			if (valuesAll[i] > 0)
				values[idx++] = valuesAll[i];

		Arrays.sort(values);

		for (int i = 0; i < posit; i++) 
		{
			sum += values[i];
			gsum += 2 * values[i] * (i + 1);
		}

		sum *= values.length;

		giniC = (gsum/sum) - ((double)(posit+1)/posit);
		posF = (double) posit / valuesAll.length;
	}


}

