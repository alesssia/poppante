/*
 *   MyTest.java
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
 *   alessia.visconti@kcl.ac.uk
 */

package com.github.alesssia.poppante;

import java.util.*;
import java.util.concurrent.*;

/**
	Performs the association/heritability tests

	@author      Alessia Visconti <alessia.visconti@kcl.ac.uk>
	@version     1.0               
	@see		 com.github.alesssia.poppante.MyTest
	@see		 com.github.alesssia.poppante.MyFileReader
*/

public class MyTest
{
	/** This is where the data are stored */
	private DataManager datamanager;
	/**  Counts the number of tests performed. */
	private int tests;   
	/**  Regulates the (concurrent) printing on the standard output */
	private final Semaphore semaphore; 
	
	/**
		Constructor. 
		
		Initialises the data structure. 
		It uses the data stored in the DataManager object.
		
		@param dm the DataManager object
	*/
	public MyTest(DataManager dm)
	{
		datamanager = dm;
		tests = 0;
		semaphore = new Semaphore(1);
	}
        
    /**
		Returns the number of tests performed.
		
		@return the number of tests
	*/
	public int tests()
	{
		return tests;
	}

	/**
		Does the analysis. 
	
		In details, it:
		- creates a vector of tests to be performed;
		- runs the tests;
		- collects and returns the results;
		
		The phenotypes to be tested are filtered according tho the information
		in the FILTER file, if present
	
		@return the results of the analysis
		@throws Exception if an error occurs during the analysis
		@see com.github.alesssia.poppante.InputTest
	*/
	public Vector<Result> analyse() throws Exception
	{	
		Vector<InputTest> inputs = new Vector<>();

		for (int m = 0; m < datamanager.listMeths().size(); m++)
		{
			//mock phenotype, it is not used in the heritability test
			if (Constants.mode == Utilities.MODE_HERITABILITY)
				inputs.add(new InputTest(Utilities.INVALID_I, m));
			else 
				for (int p = 0; p < datamanager.phenotypeNames().size(); p++) 
						inputs.add(new InputTest(p, m));
		}
	
		return processInputs(inputs);		
	}


	/**
		Prints information about the number of tests performed.
	
		The information is printed on the standard output.
	
		@param n the number of tests performed
	*/
	private void printStatistics(int n)
	{
		if (n % 10 == 0)
			System.out.print(".");
		if (n % 500 == 0)
			System.out.print(" [" + tests + "]\n\t");
	}

	/**
		Extracts the families that can be analysed.

		A family is analysable if at least one of its member
		is analysable, that is if at least one of its member
		has a valid value for the methylation site m (when the
		analysis mode is "heritability") or both the phenotype p
		and methylation site m (when the analysis mode is "association").

		@param p the position of the phenotype in the person's phenotype list
		@param m the position of the methylation site in the person's methylation site list
		@return the list of family that can be analysed
		@see com.github.alesssia.poppante.Family
		@see com.github.alesssia.poppante.Person
	*/
	private Analysable family2analise(int p, int m)
	{
		Analysable analysable = new Analysable();
		
		//used to store whether the individuals are analysable
		boolean[] isAnalysable = new boolean[datamanager.people()];
	
		//first family starts from position 0
		int start = 0;
        String[] familyKeys = datamanager.familyKeys();
		for (int f=0; f<familyKeys.length; f++)
		{
			int numAnalysable = 0;
			int offset = datamanager.families().get(familyKeys[f]).numMembers();
			
			//counts analysable family members and sets the vector of 
			//analysable people
			for(int i=start; i<start+offset; i++)
				if ((Constants.mode == Utilities.MODE_ASSOCIATION && datamanager.methylations()[m][i] != Utilities.INVALID_D && datamanager.phenotypes()[p][i] != Utilities.INVALID_D) 
					||  (Constants.mode == Utilities.MODE_HERITABILITY && datamanager.methylations()[m][i] != Utilities.INVALID_D))
				{
					isAnalysable[i] = true;
					numAnalysable++;
				}
				
			//This family has at least one analysable member
			if (numAnalysable > 0)
				analysable.addFamily(new AnalysableFamily(familyKeys[f], numAnalysable, start));
			
			//where the next family starts	
			start += offset;
		}
		
		analysable.setIsAnalysable(isAnalysable);
		
		return analysable;
	}

	/**
		Evaluates a single test.
	
		It evaluates the actual test in the VC environment. 
		If a test raises an exception, for any reason, the exception 
		is caught and added to the results.
	
		@param vc the variant component object
		@param p the position of the phenotype in the person's phenotype list
		@param mrk the methylation site
		@return the analysis results
		@see com.github.alesssia.poppante.VC
		@see com.github.alesssia.poppante.Result
	*/	
	private Result evaluate(VC vc, int p, Marker mrk)
	{
		Result result = new Result();	
	
		//Fill throws an Exception in PCA evaluation.
		//It's caught here, with those thrown by VC.solve
		try
		{
			vc.fill();
			result = vc.solve();
		}
		catch(Exception e)
		{
			result = new Result(e.getMessage());
		}

		//Information about the test that has been performed
		String phn = (p == Utilities.INVALID_I) ?  "" : datamanager.phenotypeNames().get(p).name();
		result.setInfo(phn, mrk.name(), mrk.chromosome(), mrk.position());

		return result;
	}

	/**
		Prepares and runs the tests. 
	
		PopPAnTe uses a multi-thread system -- it is possible because 
		tests are independent from one another. 		
		If the verbose mode is selected it prints a summary of
		the number of tests at run-time.
		If a list of methylation sites to analyse has been included,
		only these methylation sites are tested.
	
		@param inputs the list of tests
		@return the list of results
	
		@throws InterruptedException if an error occurs during the analysis.
		@throws ExecutionException if an error occurs during the analysis.
		@see com.github.alesssia.poppante.InputTests
	*/
	private Vector<Result> processInputs(Vector<InputTest> inputs) throws InterruptedException, ExecutionException 
	{		
		//sets threads
		Thread t = new Thread();
		ExecutorService service = Executors.newFixedThreadPool(t.threads());
		Vector<Future<Result>> futures = new Vector<>();
	
		if (Constants.verbose)
			System.out.print("\t");

		for (final InputTest input : inputs)
		{
			Callable<Result> callable = new Callable<Result>()
			{
				@Override
				public Result call() throws Exception
				{
                                        Result result;
                                        
					//extract information				
					int p = input.p();
					int m = input.m();
					
					Analysable analysable = family2analise(p, m);
										
					//no family to analise
					if (analysable.numFamilies() == 0)
                                                result = new Result("Warning : no observations");
					
					Marker mrk = datamanager.listMeths().get(m);
					VC vc = new VC(p, m, analysable, datamanager);
					
					//If the regin-based testing is selected I need to extract the sites
					//within the region
					if (Constants.region != Utilities.INVALID_I)		
						vc.setSites(mrk.getMarkersInWindow(datamanager.listMeths()));
					
					//does the test
					result = evaluate(vc, p, mrk);
					
					//critical session to count the done tests
					semaphore.acquire();
					tests++;
					int n = tests;
					semaphore.release();
				
					if (Constants.verbose)
						printStatistics(n);
				 			
					return result;
				}
			};

			futures.add(service.submit(callable));
		}

		service.shutdown();
		
		Vector<Result> outputs = new Vector<>();
		for (Future<Result> future : futures)
		{
			Result result = future.get();
			if (result != null)
				outputs.add(result);
		}
		return outputs;
	}

}


/**
	Helper class for the multi-thread test session.

	It describes a single test by means of the position of 
	the phenotype and of methylation site in the phenotypes 
	and methylation sites list, respectively
*/
class InputTest
{
	/**  The position of the phenotype in the person's phenotype list */
	int p; 
	/**  The position of the methylation site in the person's methylation site list */
	int m; 

	/** 
		Constructor. 
		
		Defines a test.
		
		@param p the position of the phenotype in the person's phenotype list
		@param m the position of the methylation site in the person's methylation site list
	*/
	public InputTest(int p, int m)
	{
		this.p = p;
		this.m = m;
	}
	
	/**
		Returns the position of the phenotype in the person's phenotype list 
		
		@return the position of the phenotype
	*/
	public int p()
	{
		return p;
	}
	
	/**
		Returns the position of the methylation site in the person's methylation site list 
		
		@return the position of the methylation site
	*/
	public int m()
	{
		return m;
	}

}

