/*
 *   PopPAnTe: Population and Pedigree Association Testing of the Epigenome
 *   Copyright (C) 2015 Dr Alessia Visconti, Dr. Mario Falchi                                  
 *
 *   PopPAnTe is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.

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
 *   alessia.visconti@kcl.ac.uk
 */

package com.github.alesssia.poppante;

import java.util.*;
import java.io.*;

/**
	PopPAnTe evaluates the heritability of DNA methylation]
	sites and their association with phenotypic traits as 
	well as with gene expression levels. PopPAnTe works with
	family data of any size and complexity, including highly 
	inbred individuals, also when the family structure is unknown.
	It allows also region-based tests by collapsing sites within 
	a region of interest. 

	It can be highly configured with a set of command line options.

	References:
		TBA
*/
public class Main 
{
	/** Version number and release data */
	private static final String year = "2018";
	private static final String version = "Version 1.0.2 -- Mar 20th, " + year;

	/**
		Return the tool information
		
		@return the tool information
	*/	
	private static String info()
	{
		String s = "\nPopPAnTe: Population and Pedigree Association Testing of Quantitative Data\n"; 
		s += version + "\n";
		s += "Copyright (C) " + year + " Dr. Alessia Visconti\t<alessia.visconti@kcl.ac.uk>\n\t\t   Dr. Mario Falchi\t<mario.falchi@kcl.ac.uk>\n";
		
		if (Constants.verbose)
		{
			s += "\nPopPAnTe is distributed in the hope that it will be useful,\n";
			s += "but WITHOUT ANY WARRANTY. See the GNU GPL v3.0 for more details.\n";
			
			s += "\nDetails about PopPAnTe can be found at:\n";
			s += "\tVisconti, A., et al., 'PopPAnTe: population and pedigree\n";
			s += "\tassociation testing for quantitative data'. BMC genomics\n";
			s += "\t18.150 (2017) https://doi.org/10.1186/s12864-017-3527-7\n"; 
			
			s += "\nPlease report comments and bugs to:\n";
			s += "\t- alessia.visconti@kcl.ac.uk\n";
		}
		
		return  s;
	}
	
	/**
		Return a help message
		
		@return a help message
	*/
	private static String help()
	{
		return  "\n" + Constants.optionList() + "\n";
	}
	

	/**
		PopPAnTe evaluates the heritability of DNA methylation]
		sites and their association with phenotypic traits as 
		well as with gene expression levels. PopPAnTe works with
		family data of any size and complexity, including highly 
		inbred individuals, also when the family structure is unknown.
		It allows also region-based tests by collapsing sites within 
		a region of interest. 

		It can be highly configured with a set of command line options.
		
		@param args command line
		@see com.github.alesssia.poppante.Constants
		@see com.github.alesssia.poppante.Result
	*/
	public static void main(String[] args) 
	{
		try 
		{
		   //  _____                               _                       _   _   _
		   // |  __ \                             | |                     | | | | (_)
		   // | |__) |_ _ _ __ __ _ _ __ ___   ___| |_ ___ _ __   ___  ___| |_| |_ _ _ __   __ _
		   // |  ___/ _` | '__/ _` | '_ ` _ \ / _ \ __/ _ \ '__| / __|/ _ \ __| __| | '_ \ / _` |
		   // | |  | (_| | | | (_| | | | | | |  __/ ||  __/ |    \__ \  __/ |_| |_| | | | | (_| |
		   // |_|   \__,_|_|  \__,_|_| |_| |_|\___|\__\___|_|    |___/\___|\__|\__|_|_| |_|\__, |
		   //                                                                               __/ |
		   //                                                                              |___/
     	
			//creates a stopwatch and starts it
			Stopwatch stopwatch = new Stopwatch();
			stopwatch.start();
			
			//thrads manager
			Thread t = new Thread();
			
			//reads options and prints those in use
			Constants.parse(args);
			System.out.println(info());
			if (Constants.help)
			{
				System.out.println(help());
				System.exit(0);
			}
			
			
			if (Constants.verbose)
			{
				System.out.println(Constants.printConstants());
				
				//Prints suhhestions on how to preceed if the system does not converge
				if (Constants.kinship != null && Constants.decomposition == null)
					System.out.println("\n\tIn the unluckily case PopPAnTe raises exception during the\n\tvariant component evaluation please try with the option\n\t --decomposition LU or --decomposition QR\n");
				else if (Constants.kinship != null && !Constants.decomposition.equals("QR"))
					System.out.println("\n\tIn the unluckily case PopPAnTe raises exception during the\n\tvariant component evaluation please try with the option\n\t --decomposition QR\n");
				
				if (t.threads() > 1)
					System.out.println("Multithread version used. " + t.threads() + " threads running.\n");
				System.out.println("Loading data...\n");
			}
		
						
			//creating our workhorse and the file manager
			DataManager datamanager = new DataManager();
			MyFileReader filereader = new MyFileReader(datamanager);
			

			
 //			    _____  ______ _____         _               _
 // 		   |  __ \|  ____|  __ \       | |             | |
 // 		   | |__) | |__  | |  | |   ___| |__   ___  ___| | _____ _ __
 // 		   |  ___/|  __| | |  | |  / __| '_ \ / _ \/ __| |/ / _ \ '__|
 // 		   | |    | |____| |__| | | (__| | | |  __/ (__|   <  __/ |
 // 		   |_|    |______|_____/   \___|_| |_|\___|\___|_|\_\___|_|
 //
                                                          
			if(Constants.mode == Utilities.MODE_PEDCHECK)
			{
				try
				{
					System.out.println("Checking pedigree file");
					
					//I do not need to set the phenotypes number, also
					//because phenotypes values aren't read in this setting 
					FamilySet familyset = filereader.readFamilies();
					String message = datamanager.checkFamilies(familyset.message());
					if (message.equals(""))
						System.out.println("Pedigree correct");
					else
					{
						System.out.println("The pedigree checker has found the following problems: ");
						System.out.println(message);
					}
				}
				catch(IOException | NotWellFormedLineException e)
				{
					System.out.println(e.getMessage());
					System.out.println();
				}
				
				stopwatch.stop();
				if (Constants.verbose)
					System.out.println("Pedigree assessed in " + stopwatch.getTime());
				
				return;
			}
	
			
			
		// 			_____        _          _                 _ _
	    // 		   |  __ \      | |        | |               | (_)
	    // 		   | |  | | __ _| |_ __ _  | | ___   __ _  __| |_ _ __   __ _
	    // 		   | |  | |/ _` | __/ _` | | |/ _ \ / _` |/ _` | | '_ \ / _` |
	    // 		   | |__| | (_| | || (_| | | | (_) | (_| | (_| | | | | | (_| |
	    // 		   |_____/ \__,_|\__\__,_| |_|\___/ \__,_|\__,_|_|_| |_|\__, |
	    // 		                                                         __/ |
	    // 		                                                        |___/
			
			//Time to read the first set of data, that is the PED file and
			//the phenotypes values. Since readFamilies chechs whether the
			//number of loaded phenotypes is correct, the phenotype list
			//should be loaded beforehand.
			//The phenotype values are loaded only if the mode is association.
			//In fact they are useless otherwise. In this case the phenotype
			//number, used to chech the validity of the PED file, is set to 
			//zero.
			//The phenotype values are temporary loaded within the Person
			//object, that represents one of population member. This is due 
			//to the fact that the individuals are not sorted yet and the
			//phenotype matrix in datamanager can't be loaded properly.
			//In fact, the order of the data loaded in the phenotype matrix
			//(as well as that of the data loaded in the methylation matrix)
			//depends on the order of the object of the population member.
			//The loaded IDs correspond to the real ones when the family 
			//structure is available and to mock ones when the family structure
			//is missing.
			int numPhenotypes = 0;
			int numPhenos2analyse = 0;
			if (Constants.filter != null) 
				numPhenos2analyse = filereader.readIncludedPhenos();
			if (Constants.mode == Utilities.MODE_ASSOCIATION)
				numPhenotypes = filereader.readPhenotypesInformation();
			int people =  filereader.readFamilies().people();
			
			//Once the families are read, if a family structure is available
			//(that is I am estimating the kinship matrix from the genealogical
			//data), I need to sort the individual within families BEFORE
			//initialising the positionTable, and evaluating the kinship matrix.
			//If the theoretical kinship values are used it is also time to
			//evaluate them. If the external matrix is used, this job is deferred 
			//so the kinship matrix does not need to be re-organised once 
			//individuals are removed (see afterwards)
			if (Constants.kinship == null)
			{
				
				datamanager.sort();
				datamanager.initialiseKinships();
				datamanager.evaluateKinship();
			}
			

			//Since the mock parents are no loger used (I know who the mock 
			//individuals are, beacuse they have been marked as such during
			//the reading of the pheonotypes -- they do not have any valid
			//phenotype value). I can remove them from the dataset.
			//This will save both space and time for the analysis.
			//Removing individuals may create empty families. If they exist
			//they are removed as well (family composed only by mock individuals
			//are empty families).
			datamanager.removeMock();
			datamanager.removeEmptyFamilies();
			
			//I then read the covariates and I store them in a temporary structure
			//witin each person. I then remove those individuals having
			//a missing value for at least one covariate. In fact these individuals
			//can't be analysed. This will again save both space and time for the 
			//analysis. Also in this case the empty families will be removed.
			//To correctly read the covariate files I assume that their number
			//is equal to the first line of the covariate file. I can also 
			//create the position table, that will be used by the kinship
			//and methylation reading. 
			if (Constants.covariate != null) 
			{
				filereader.setNumCovar();
				filereader.readCovariates();
				datamanager.removeMissingCovariates();
				
				datamanager.removeEmptyFamilies();
				datamanager.updatePositionTable();
				datamanager.resetCovariates();
			}
			else	
				datamanager.updatePositionTable();

			int numcov = datamanager.numCovar();
			//Now that I removed all the individuals I no loger need, I can
			//update the phenotypes and covariate structure, and remove the data 
			//temporary loaded within the single people (actually covariates have 
			//been updatete beforehand). 
			datamanager.resetPhenotypes();
			
			//time to read the external kinship values. The kinship matrix should be 
			//initialised and the diagonal values set. The actual values can be read
			//afterwards. 
			//The user may also have been chosen to use a bending procedure to make 
			//the external kinship matrix positive semidefinite (and thus solvable 
			//with the Cholesky decomposition, instead of using the more computational
			//expensive LU/QR decomposition. No actual check on positive definiteness is
			//acctually performed
			if (Constants.kinship != null)
			{
				datamanager.initialiseKinships();
				datamanager.adjustKinship();
				filereader.readKinship();	
				
				if (Constants.decomposition == null)
					datamanager.bending();
			}

			//Time to read the second set of data, that is the methylation
			//values and information. Since readMethylationData chechs whether 
			//the number of loaded methylation values is correct, the list of
			//methylation sites should be loaded beforehand.
			//Also, if the user decide of testing only a set of methylation 
			//sets, only these values are loaded in both the list of methylation
			//sites and in the matrix of methylation values. 
			//This is why the include file is read before any other file and the
			//dimension of the methylation array is set immedetialy after.
			//Please note that the list of methylation sites to analyse may
			//include also sites that are not in the dataset.
			int numSites2analyse = 0;
			if (Constants.include != null) 
				numSites2analyse = filereader.readIncludedMeths();
			int numSites = filereader.readMethylationInformation();
			
			datamanager.resetMethylations();
			filereader.readMethylationData();	
			
			//If the user provides also covariates for the correction
			//of the methylation values, they must be read as well.
			int correctionValues = 0;
			if (Constants.correct != null && !Utilities.isDouble(Constants.correct))
			{
				correctionValues = filereader.getNumCorrectionCov();
				datamanager.resetCorrectionCovs(correctionValues);
				filereader.readCorrectionCovariates();
			}
			stopwatch.stop();
		
			//prints data information
			if (Constants.verbose)
			{
				System.out.println("\t" + people + "\tindividuals\tread from [ " + Constants.ped + " ]");
				if (Constants.kinship == null)
					System.out.println("\t" + datamanager.numFamilies() + "\tfamilies\tread from [ " + Constants.ped + " ]");
				else 
				{
					System.out.println("\t \trelationships\tread from [ " + Constants.kinship + " ]");
					System.out.println("\t\t\t\t(unrelated individuals)");
				}
			
				if (Constants.mode == Utilities.MODE_ASSOCIATION)
		    		System.out.println("\t" + numPhenotypes + "\tresponse\tread from [ " + Constants.response + " ]");
				if (Constants.mode == Utilities.MODE_ASSOCIATION && Constants.filter != null)
					 System.out.println("\t" + numPhenos2analyse + "\tresponses\thave been selected from [ " + Constants.filter + " ]");
				if (Constants.include != null)
					 System.out.println("\t" + numSites2analyse + "\tpredictors\thave been selected from [ " + Constants.include + " ]");
				System.out.println("\t" + numSites + "\tpredictors\tloaded from [ " + Constants.predictor + " ]");
				
		
				if (Constants.covariate != null) 
					System.out.println("\t" + datamanager.numCovar() + "\tcovariates\tread from [ " + Constants.covariate + " ]");
			
				if (Constants.correct != null && !Utilities.isDouble(Constants.correct))
					System.out.println("\t" + correctionValues + "\tcovariates (correction)\tread from [ " + Constants.correct + " ]");
			
				System.out.println("\nData loaded in " + stopwatch.getTime());
			}

 		   //  _   _                            _ _           _   _
 		   // | \ | |                          | (_)         | | (_)
 		   // |  \| | ___  _ __ _ __ ___   __ _| |_ ___  __ _| |_ _  ___  _ __
 		   // | . ` |/ _ \| '__| '_ ` _ \ / _` | | / __|/ _` | __| |/ _ \| '_ \
 		   // | |\  | (_) | |  | | | | | | (_| | | \__ \ (_| | |_| | (_) | | | |
 		   // |_| \_|\___/|_|  |_| |_| |_|\__,_|_|_|___/\__,_|\__|_|\___/|_| |_|
               
			  	
 			//normalises phenotypes
 			if (Constants.normalise != null && (Constants.normalise.equals("response") || Constants.normalise.equals("both")))
 			{
 				stopwatch.reset();
 				stopwatch.start();
 				if (Constants.verbose)
 					System.out.println("Normalising response data...");
 				try
 				{
 					datamanager.inverseNormalTransform("response");
 				}
 				catch(Exception e)
 				{
 					System.out.println("ERROR during data normalisation.\nNormalisation can't be performed");
 					System.out.println();
 					return;
 				}
 				stopwatch.stop();
 				if (Constants.verbose)
 					System.out.println("Data normalised in " + stopwatch.getTime());
 			}
			
			//normalises methylation values
			if (Constants.normalise != null && (Constants.normalise.equals("predictor") || Constants.normalise.equals("both")))
			{
				stopwatch.reset();
				stopwatch.start();
				if (Constants.verbose)
					System.out.println("Normalising predictor data...");
				try
				{
					datamanager.inverseNormalTransform("predictor");
				}
				catch(Exception e)
				{
					System.out.println("ERROR during data normalisation.\nNormalisation can't be performed");
					System.out.println();
					return;
				}
				stopwatch.stop();
				if (Constants.verbose)
					System.out.println("Data normalised in " + stopwatch.getTime());
			}
			
			
  		  //   _____                         _   _
  		  //  / ____|                       | | (_)
  		  // | |     ___  _ __ _ __ ___  ___| |_ _  ___  _ __
  		  // | |    / _ \| '__| '__/ _ \/ __| __| |/ _ \| '_ \
  		  // | |___| (_) | |  | | |  __/ (__| |_| | (_) | | | |
  		  //  \_____\___/|_|  |_|  \___|\___|\__|_|\___/|_| |_|
           	
  			//correct the methylation values
  			if (Constants.correct != null)
  			{
  				stopwatch.reset();
  				stopwatch.start();
  				if (Constants.verbose)
  					System.out.println("Correcting predictors values...");
				
  				int n = datamanager.correct();
  				if (Constants.verbose)
  					System.out.println("Correction was performed using " + n + " variables");
				
  				stopwatch.stop();
				if (Constants.verbose)
	  				System.out.println("Data corrected in " + stopwatch.getTime());
  			}
		
 		   //  _______        _   _
 		   // |__   __|      | | (_)
 		   //    | | ___  ___| |_ _ _ __   __ _
 		   //    | |/ _ \/ __| __| | '_ \ / _` |
 		   //    | |  __/\__ \ |_| | | | | (_| |
 		   //    |_|\___||___/\__|_|_| |_|\__, |
 		   //                              __/ |
 		   //                             |___/
		
			//prints info
			if (Constants.verbose)
			{
				System.out.println();
				if (Constants.mode == Utilities.MODE_HERITABILITY)
					System.out.print("Heritability analysis selected, ");
				else
					System.out.print("Association analysis selected, ");

				if (Constants.region != Utilities.INVALID_I)
					System.out.println("Region-based testing chosen.");
				
			    if (Constants.alpha != Utilities.INVALID_D)
			   		System.out.println("An adaptive permutations procedure will be performed to assess\nthe significance of the tests.");
					
				System.out.println("\nAnalysis started. Each dot corresponds to 10 tests.\n") ;
			}
				                                     
			//does analysis
			stopwatch.reset();
			stopwatch.start();
			
			//I need to initialise the missingness pattern for both phenotypes
			//and methylations. I need also to create an empty null model for each 
			// phenotype. In fact each phenotype has is own missingness pattern
 			//and some methylation sites may have it as well, and the same
 			//null model can be used to assess multiple sites.
			//This is done to speed up the evaluation of the variance component
			//system, because the null model will be common among these "datasets".
			if (Constants.mode == Utilities.MODE_ASSOCIATION)
			{
				datamanager.setMissingnessMethPattern();
				datamanager.setMissingnessPhenoPattern();
			}

			
			MyTest test = new MyTest(datamanager);
			Vector<Result> res = test.analyse();
			stopwatch.stop();
		
			if (Constants.verbose)
				System.out.println("\n\nAnalysis ended.\n\t" + test.tests() + " tests performed in " + stopwatch.getTime() + ".\n") ;
			
					
		   //  _____          _                                            _
		   // |  __ \        | |                                          (_)
		   // | |__) |__  ___| |_ ______ _ __  _ __ ___   ___ ___  ___ ___ _ _ __   __ _
		   // |  ___/ _ \/ __| __|______| '_ \| '__/ _ \ / __/ _ \/ __/ __| | '_ \ / _` |
		   // | |  | (_) \__ \ |_       | |_) | | | (_) | (_|  __/\__ \__ \ | | | | (_| |
		   // |_|   \___/|___/\__|      | .__/|_|  \___/ \___\___||___/___/_|_| |_|\__, |
		   //                           | |                                         __/ |
		   //                           |_|                                        |___/
			
			//does FDR correction and plots manhattan and qqplot
			if (Constants.verbose)
				System.out.println("Correcting for multiple comparisons...");
			
			String plotNames = "";
			
			stopwatch.reset();
			stopwatch.start();
			ResultTable rt = new ResultTable(res, false);
			rt.adjustPvalues(false);
			
			if (Constants.plot)
			{
				plotNames += rt.qqplot(false) + "\n";
				plotNames += "\t\t\t" + rt.manhattan(false) + "\n";
			}
			
			//empirical p-values are corrcted and plotted only if necessary
			if (Constants.alpha != Utilities.INVALID_D)
			{
				rt = new ResultTable(res, true);
				rt.adjustPvalues(true);
				if (Constants.plot)
				{
					plotNames += "\t\t\t" + rt.qqplot(true) + "\n";
					plotNames += "\t\t\t" + rt.manhattan(true) + "\n";
				}
			}
			stopwatch.stop();
			
			if (Constants.verbose)
			{
				System.out.println("Correction performed in " + stopwatch.getTime() + ".") ;
				if (Constants.plot)
					System.out.println("\nPlots saved as\t" + plotNames + ".") ;
			}
			
			
		   //  _____        _                        _ _   _
		   // |  __ \      | |                      (_) | (_)
		   // | |  | | __ _| |_ __ _  __      ___ __ _| |_ _ _ __   __ _
		   // | |  | |/ _` | __/ _` | \ \ /\ / / '__| | __| | '_ \ / _` |
		   // | |__| | (_| | || (_| |  \ V  V /| |  | | |_| | | | | (_| |
		   // |_____/ \__,_|\__\__,_|   \_/\_/ |_|  |_|\__|_|_| |_|\__, |
		   //                                                       __/ |
		   //                                                      |___/
			
			//prints results
			if (Constants.verbose)
				System.out.println("Writing results...");
			
			if (Constants.output != null)
				System.setOut(new PrintStream(new FileOutputStream(new File(Constants.output + ".tsv"))));
			else
				System.out.println();
		
			if (Constants.header)
				System.out.println(Result.header(datamanager.methFields()));
			
			stopwatch.reset();
			stopwatch.start();
			for (Enumeration<Result> elements = res.elements(); elements.hasMoreElements();)
				System.out.println(elements.nextElement());
			stopwatch.stop();
		
			if (Constants.output != null)
			{
				System.setOut(new PrintStream(new FileOutputStream(FileDescriptor.out)));
				System.out.print("Results written in " + Constants.output + ".tsv");
				if (Constants.verbose)
					System.out.println(" in " + stopwatch.getTime());
				else
					System.out.println("\n");
			}
			else if (Constants.verbose)
				System.out.println("\nResults written in " + stopwatch.getTime() + ".") ;
			
			if (Constants.verbose)
				System.out.println("PopPAnTe ends.\n");	
				
		} 
		catch (Exception e)  //something bad happened
		{
			//e.printStackTrace();
			System.out.println(e.getMessage());
			System.out.println();
		}
	}

}
