/*
 *   DataManager.java
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

import Jama.Matrix;
import com.github.alesssia.algebrautils.MyLinearRegression;
import java.io.IOException;
import java.util.*;

import org.apache.commons.math3.stat.regression.MultipleLinearRegression;
import org.hashids.*;

import com.github.alesssia.arrayutils.*;
import com.github.alesssia.algebrautils.MyPCA;


/**
	Stores and manages the data.
	
	It is responsible of storing the data (loaded using the MyFileReader class), 
	using mainly three objects:
		- families, that stores the family structures and the individual 
		  information, along with individuals' phenotypes, methylation values,
		  and covariates
		- phenotypeNames, that stores the list of phenotype values as they
	      stored in the Person objects
	 	- listMeths, that stores the list of methylation values as they
	      stored in the Person objects.

	It performs data transformation, that is:
		- quantile normalisation of the phenotype values
		- correction of the methylation values

	It is used as parameter for the association/heritability tests, that 
	are delegated to the class Test.

	@author      Alessia Visconti <alessia.visconti@kcl.ac.uk>
	@version     1.0               
	@see		 com.github.alesssia.poppante.Test
	@see		 com.github.alesssia.poppante.MyFileReader
 */

class DataManager 
{
	/**  Family and subjects information. */
	private Hashtable<String, Family> families; 
	/**  Mapping between individuals' and their position.
		
		 It is used to store the absolute position of every values
		 loaded and it is initialised after the updatePositionTable is called.  */	
	private Hashtable<String, Integer> positionTable;
	/** Order in which the families should be analysed */
	private String[] familyKeys;
	
    /**  Phenotypes to analyses. */		 
	private Vector<String> includedPheno;
	/**  Phenotype information. */
	private Vector<Phenotype> phenotypeNames;  
	/** Phenotypes values.
		
		Each row represents a phenotype and each column represents 
		an individual. Order of the phenotypes is done by the 
		phenotypeNames list, while order of the individuals is done
		according to the position in positionTable. */
	private double[][] phenotypes;
		 
		 
	/**  Methylation sites to analyses. */		 
	private Vector<String> includedMeth;
	/**  Methylation site map. */
	private Vector<Marker> listMeths; 
	/** Methylation values.
		
		Each row represents a methylation site and each column represents 
		an individual. Order of the sites is done by the 
		listMeths list, while order of the individuals is done
		according to the position in positionTable. */
	private double[][] methylations;
	
	
	/** Covariate values.
		
		Each row represents a covariate and each column represents 
		an individual. Order of the individuals is done
		according to the position in positionTable. */
	private double[][] covariates;
	
	
	/** Number of covariates/ */
	private int numCovar;

	/** Covariates (for correction of methylation sites).
		
		Each row represent a variable used for correction and 
		each column represents an individual. Order of the 
		individuals is done according to the position in 
		positionTable. */
	private double[][] correctionCovs;
	
	/** Unique hash used to define pattern of missingness
		of methylation data.
		
		It converts the list of missing sites for each
		person in a unique key, that will be used to
		analyse the people having the same pattern of
		missingness together. */
	private Hashids hashids;
	
	/** Assignes each methylation site to the set that 
		identifies its missingness pattern */
	private String[] missingnessMethPattern;
		
	/**	Set of missingness pattern identified 
		for the phenotypes and the methylation
		sites and their null log likelihood.
		They hash key represent the phenotype (each
		phenotype will have its null model) and the
		value is the null model itself */
	private Hashtable<Integer, NullModel>  missingnessPattern;
	
	/** Number of fields used to describe the
		methylation values */
	private int methField;

	
	/**
		Constructor. 
		
		Initialises the data structure.
	*/
	public DataManager()
	{
		families = new Hashtable<>();
		positionTable = new Hashtable<>();
		familyKeys = null;
		
		phenotypeNames = new Vector<>();
                includedPheno = new Vector<>();
		phenotypes = null;
		
		listMeths = new Vector<>();
		includedMeth = new Vector<>();
		methylations = null;
		methField = 1;
		
		covariates = null;
		numCovar = 0;
		correctionCovs = null;
		
		//Variable used for the missingness patterns
		hashids = new Hashids("42"); //set salt
		missingnessMethPattern = new String[0];
		missingnessPattern = new Hashtable();
	}
	
	/**
		Returns all the information loaded for family members,
		along with the family structure.
		
		@return the family structure
	*/
	public Hashtable<String, Family> families()
	{
		return families;
	}
	
	
	/**
		Returns the number of families that have been loaded.
		
		@return the number of families
	*/
	public int numFamilies()
	{
		return families.size();
	}
	
	/**
		Returns the number of individual that have been loaded.
		
		@return the number of individuals
	*/
	public int people()
	{
		return positionTable.size();
	}
	
	/**
		Returns the sorted list of family keys 
		
		@return the list of family keys
	*/
	public String[] familyKeys()
	{
		return familyKeys;
	}
        
    /**
		Returns the list of loaded phenotypes
		
		@return the the list of loaded phenotypes
	*/
	public Vector<Phenotype> phenotypeNames()
	{
		return phenotypeNames;
	}
        
    /**
		Returns the list of phenotypes sites to be analysed.
		
		@return the list of phenotypes sites to be analysed.
	*/
	public Vector<String> includedPheno()
	{
		return includedPheno;
	}
	
    /**
		Returns the list of loaded methylation sites.
		
		@return the the list of loaded methylation sites.
	*/
	public Vector<Marker> listMeths()
	{
		return listMeths;
	}

    /**
		Returns the list of methylation sites to be analysed.
		
		@return the list of methylation sites to be analysed.
	*/
	public Vector<String> includedMeth()
	{
		return includedMeth;
	}
	
    /**
		Returns the number of fields used to describe a methylation
		site.
		
		@return number of fields
	*/
	public int methFields()
	{
		return methField;
	}
	
    /**
		Sets the number of fields used to describe a methylation
		site.
		
		@param n number of fields
	*/
	public void setMethFields(int n)
	{
		methField=n;
	}
	
    /**
		Returns whether the correctionCovs table is initialised
		
		@return whether the correctionCovs table is initialised
	*/
	public boolean isCorrectionCovsInitialised()
	{
		return correctionCovs != null && correctionCovs.length != 0 && correctionCovs[0].length != 0;
	}
	
    /**
		Returns whether the methylations table is initialised
		
		@return whether the methylations table is initialised
	*/
	public boolean isMethylationsInitialised()
	{
		return methylations != null && methylations.length != 0 && methylations[0].length != 0;
	}
	
    /**
		Returns the table of methylation values
		
		@return the table of methylation values
	*/
	public double[][] methylations()
	{
		if (methylations == null)
			return null;
		
		return methylations.clone();
	}
	
	
    /**
		Returns the table of phenotype values
		
		@return the table of phenotype values
	*/
	public double[][] phenotypes()
	{
		if (phenotypes == null)
			return null;
		
		return phenotypes.clone();
	}
	
    /**
		Returns the table of covariates values
		
		@return the table of covariates values
	*/
	public double[][] covariates()
	{
		if (covariates == null)
			return null;
		
		return covariates.clone();
	}
	
    /**
		Returns the table of phenotype values
		
		@return the table of phenotype values
	*/
	public double[][] correctionCovs()
	{
		if (correctionCovs == null)
			return null;
		
		return correctionCovs.clone();
	}
	
	
	/**
		Sets the specified methylation site to the given value
		
		@precondition the request position is valid
		
		@param row the row in the methylation matrix
		@param col the column in the methylation matrix
        @param value the value to set
	*/
	public void setMeth(int row, int col, double value)
	{
		assert row >= 0 && row < methylations.length : "Internal error: invalid predictor position";
		assert col >= 0 && col < methylations[0].length : "Internal error: invalid predictor position";
		
		methylations[row][col] = value;
	}
	
	
	/**
		Sets the specified correction covariate value to the given value
		
		@precondition the request position is valid
		
		@param row the row in the methylation matrix
		@param col the column in the methylation matrix
        @param value the value to set
	*/
	public void setCor(int row, int col, double value)
	{
		assert row >= 0 && row < correctionCovs.length : "Internal error: invalid correctionCovs position";
		assert col >= 0 && col < correctionCovs[0].length : "Internal error: invalid correctionCovs position";
		
		correctionCovs[row][col] = value;
	}
	
	
	/**
		Sets the number of covariates that have been loaded.
		
		@param n the number of covariates
	*/
	public void setNumCovar(int n)
	{
		numCovar = n;
	}
	
	/**
		Returns the number of covariates that have been loaded.
		
		@return the number of covariates
	*/
	public int numCovar()
	{
		return numCovar;
	}
	

	/**
		Returns the position of the given individual
		in the positionTable.
		
		Individuals are codified as the concatenation
		of their famID with their ID.
                
        @param individual the individual 		
		@return the position of the individual
	*/
	public Integer getPosition(String individual)
	{
		return positionTable.get(individual);
	}
	
	
	/**
		Sorts the individuals within familes
		
		The sorting of the individuals is a precodintion of 
		the kinship evaluation	
	*/
	public void sort()
	{
		for (Enumeration<Family> elements = families.elements(); elements.hasMoreElements();)
			elements.nextElement().sort();
	}
	
	
	/**
		Updates (or creates) the position table object.
		
		The position table determines the position of any individual value  
		within the data structures used by Poppante (i.e., covariates, 
		phenotypes and methylation values, and kinship matrix).
		For instance the phenotype (methylation) values for the individual in 
		the i-th position according to the position table will be in the i-th 
		column of the phenotypes (methylations) matrix.
		
		It is also necessary to map the alfanumeric IDs to numeric one.
		
		It must be created once the family are read from the data files,
		and the updated again every time an individual is removed from 
		the population (e.g., when she misses one of her covariates).
	*/
	
	public void updatePositionTable()
	{
		
		positionTable = new Hashtable<>();
		familyKeys = new String[families.size()];
		
		//Sorting is important when real family structures are available.
		//In this case, indeed, the family members should have a specific
		//order, as request by the Kennet Lange's method for kinship estimation
		if (Constants.kinship == null)
		{
			//I can't updatePositionTable the HashTable, but I can create a sorted vector of keys.
			Object [] tmp = families.keySet().toArray();
			familyKeys = Arrays.copyOf(tmp, tmp.length, String[].class);
			Arrays.sort(familyKeys);
		}
		//When individuals are unrelated (and only one family is available)
		//a real "updatePositionTable" is not necessary, and any order of the memberIDs list
		//is fine.
		else
		{
			familyKeys = new String[1];
			familyKeys[0] = Utilities.INVALID_S;
		}
		
		//Time to determine the position of every individual in the
		//analysis. The position is determined both by the family order
		//in the family keys and, within each family, by the order of 
		//the member vector array.
		int position = 0;
		for (int f=0; f<familyKeys.length; f++)
		{
			String famID = familyKeys[f];
			Vector<Person> members = families.get(famID).members();
			for (Enumeration<Person> elements = members.elements(); elements.hasMoreElements();)
			{
				Person p = elements.nextElement();
				positionTable.put(new String(famID + "" + p.id()), position);
				position++;
			}
		}
	}
	
	/**
		Removes individuals with missing covariates.
		
		This can be done safely because individuals with 
		a missing value for a covariate are not analysed.
		Mock individuals have been set to have all covariates, 
		and they are retained (they are used for the evaluation
		of the kinship matrix).
		
		Please note that this methods does not remove 
		family that may not have any (real) individual left
		after this test.
		
		@precondtion the number of covariates should be greater than zero
	*/
	public void removeMissingCovariates()
	{
		assert numCovar > 0 : "Internal error: removeMissingCovariates() called on individual with no covariates read";
		
		for (Enumeration<Family> elements = families.elements(); elements.hasMoreElements();)
		{
            Family family = elements.nextElement();
			Vector<Person> members = family.members();
			
			//select person to be removed
			Vector<Person> marked = new Vector<>();
			Vector<Integer> position = new Vector<>();
                        
			for (int i=0; i<members.size(); i++)
				if (members.get(i).missesCovariates())
				{
					marked.add(members.get(i));
					position.add(i);
				}
			
			//remove them
			for (Enumeration<Person> people = marked.elements(); people.hasMoreElements();)
				members.remove(people.nextElement());
			
			//If the kinship matrix has been evaluated I need to remove the 
			//people also from here (if any has been removed)
			if (Constants.kinship == null && position.size() > 0)
				family.kinship().reset(position, members.size());
		}
	}
	
	/**
		Removes individuals with missing covariates.
		
		This can be done safely because individuals with 
		a missing value for a covariate are not analysed.
		Mock individuals have been set to have all covariates, 
		and they are retained (they are used for the evaluation
		of the kinship matrix).
		
		Please note that this methods does not remove 
		family that may not have any (real) individual left
		after this test.
	*/
	public void removeMock()
	{
		for (Enumeration<Family> elements = families.elements(); elements.hasMoreElements();)
		{
			Family family = elements.nextElement();
			Vector<Person> members = family.members();
			
			//select person to be removed
			Vector<Person> marked = new Vector<>();
			Vector<Integer> position = new Vector<>();
			for (int i=0; i<members.size(); i++)
				if (members.get(i).isMock())
				{
					marked.add(members.get(i));
					position.add(i);
				}
			//remove them
			for (Enumeration<Person> people = marked.elements(); people.hasMoreElements();)
				members.remove(people.nextElement());
			
			//If the kinship matrix has been evaluated I need to remove the 
			//people also from here (if any has been removed)
			if (Constants.kinship == null && position.size() > 0)
				family.kinship().reset(position, members.size());
		}
	}
	
	
	
	/**
		Removes empty families.
		
		A family is considered empty either if it has no member or
		if it includes only mock individuals.	
	*/
	public void removeEmptyFamilies()
	{
		Vector<String> marked = new Vector<>();
		
		Object [] tmp = families.keySet().toArray();
		String[] keys = Arrays.copyOf(tmp, tmp.length, String[].class);
		
		for (int f=0; f<keys.length; f++)
		{
			//check if there is at least one memeber that is not mock
			Vector<Person> memebers = families.get(keys[f]).members();
			int mocks = 0;
			Enumeration<Person> elements = memebers.elements();
			while (elements.hasMoreElements() && elements.nextElement().isMock() )
				mocks++;
			
			//select family to remove (it removes also familues with zero members
			//(in that case mock is zero and member.size() too).
			if (memebers.size() == mocks)
				marked.add(keys[f]);
		}
		
		//remove them
		for (Enumeration<String> elements = marked.elements(); elements.hasMoreElements();)
			families.remove(elements.nextElement());
		
	}
	
	
	/**
		Sets the phenotype values that have been stored 
		within each individual in the correct data structure.
		
		The temporary data are also removed to free memory space.
		
		@precondition the phenotypes table should have been initialised
		@precondition the position table should have been initialised
	*/
	public void resetPhenotypes()
	{
		assert phenotypeNames.size() > 0 : "Internal error: response list not initialised";
		assert positionTable.size() > 0 : "Internal error: position table not initialised";
		
		phenotypes = new double[phenotypeNames.size()][positionTable.size()];
		
		//any order is fine. The real order depends on the position table
		for (Enumeration<Family> elements = families.elements(); elements.hasMoreElements();)
		{
			Vector<Person> members = elements.nextElement().members();
			for (Enumeration<Person> people = members.elements(); people.hasMoreElements();)
			{
				Person person = people.nextElement();			
				int position = getPosition(new String(person.famID() + "" + person.id()));
				double [] phenos = person.phenotypes();
				for (int i=0; i<phenos.length; i++) 
				{
					phenotypes[i][position] = phenos[i];
					if (phenos[i] == Utilities.INVALID_D)
						phenotypeNames().get(i).addMissing((long)position);
				}
				person.resetPhenotypes();	
			}
		}
	}
	
	
	
	/**
		Sets the dimension of the methylation values matrix,
		and initialises its values to not valid (missing).
		In fact missing values are those that are not read/written
		and need to be set beforehand.
		
		@precondtion the number of methylation sites should be greater than zero
		@precondition the position table should have been initialised
	*/
	public void resetMethylations()
	{
		assert listMeths.size() > 0 : "Internal error: no predictor information available";
		assert positionTable.size() > 0 : "Internal error: position table not initialised";
		
		//FIXME: Bottleneck
		methylations = Utilities.set(listMeths.size(), positionTable.size(), Utilities.INVALID_D);
	}
	
	/**
		Sets the dimension of the methylation values matrix,
		and initialises its values to not valid (missing).
		In fact missing values are those that are not read/written
		and need to be set beforehand.
		
		@precondition the position table should have been initialised
		
		@param nCorrCov the number of correction covariates
	*/
	public void resetCorrectionCovs(int nCorrCov)
	{
		assert positionTable.size() > 0 : "Internal error: position table not initialised";
		
		//FIXME: Bottleneck
		correctionCovs = Utilities.set(nCorrCov, positionTable.size(), Utilities.INVALID_D);
	}
	
	/**
		Sets the covariate values that have been stored 
		within each individual in the correct data structure.
		
		The temporary data are also removed to free memory space.
		@precondtion the number of covariates should be greater than zero
		@precondition the position table should have been initialised
		@precondition individuals with missing covariates should have been removed
		
        @throws java.lang.Exception if there are individuals with covariates that haven't removed beforehands
	*/
	public void resetCovariates() throws Exception
	{
		assert numCovar > 0 : "Internal error: no covariates available";
		assert positionTable.size() > 0 : "Internal error: position table not initialised";
			
		covariates = new double[numCovar][positionTable.size()];

		//any order is fine. The real order depends on the position table
		for (Enumeration<Family> elements = families.elements(); elements.hasMoreElements();)
		{
			Vector<Person> members = elements.nextElement().members();
			for (Enumeration<Person> people = members.elements(); people.hasMoreElements();)
			{
				Person person = people.nextElement();
                double [] covars = person.covariates();
                if (covars == null)
					throw new Exception("Internal error: individuals with missing covariates should be removed before covariate resetting.");
				                
				int position = getPosition(new String(person.famID() + "" + person.id()));
				for (int i=0; i<covars.length; i++) 
					covariates[i][position] = covars[i];
				
				person.resetCovariates();	
			}
		}
	}
			

	
	/**
		Checks the pedigree information.
		
		It verifies whether the first seven columns are in 
		PLINK format and whether the information are consistent. 
		A message listing all the error is returned.
		Used when mode is "pedcheck".
		
		@param m the error message created during the loading of the ped file
		@return the error message including issues detected within families
		@throws IOException If the pedigree file can't be read
		@throws NotWellFormedLineException If the file has a line that is not well-formed
		@see com.github.alesssia.poppante.MyFileReader
	*/
	public String checkFamilies(String m) throws IOException, NotWellFormedLineException
	{		
        String message = m;
		for (Enumeration<Family> elements = families.elements(); elements.hasMoreElements();)
			message += elements.nextElement().checkFamily();
		
		return message;
	}
	
	
	/**
		Evaluates the kinship matrix using the Kenneth Lange's 
		approach.
		
		If no kinship matrix is supplied the algorithm evaluates 
		a kinship matrix is evaluated for each family.
	
		@precondition the kinship matrix should have been initialised
		@precondition the family data must have been sorted
		
	*/
	public void evaluateKinship()
	{
		assert Constants.kinship == null : "Internal error: kinship should be evaluated only when not provided.";
		
		for (Enumeration<Family> elements = families.elements(); elements.hasMoreElements();)
			elements.nextElement().evaluateKinship();
	}	
	
	
	/**
		Adjusts the self-kinship values.
		
		Since self kinship values are not necessary always available
		in the kinship matrix provided (for instance PLINK --genome 
		does not evaluate them), they are then set now to 1. 
		If self-kinship values are available, they will be overwritten.
	*/
	public void adjustKinship()
	{
		families.get(Utilities.INVALID_S).kinship().adjustKinship();
	}
	
	/**
		Initialises the families' kinship matrices.
		
		This is done to:
		- satisfy the precondition of kinship evaluation according
		  to the Kenneth Lange's approach.
		- to set up the kinship structure needed when the kinship
		  is read from an external file.
	*/
	public void initialiseKinships()
	{	
		for (Enumeration<Family> elements = families.elements(); elements.hasMoreElements();)
		{
			Family family= elements.nextElement();
			family.kinship.init(family.members().size());
		}
	}
	
	
	/**
		Applies the 'bending' procedure to modify the eigenvalues of non
		positive defined (external) kinship matrix.
	
		See http://www.aps.uoguelph.ca/~lrs/Summer2012Full/PDforce.pdf for
		the R function that inspired this code
	*/
	public void bending()
	{	
		assert Constants.kinship != null : "Internal error. Bending is available only for external kinship matrix.";
		
		Family family = families().get(Utilities.INVALID_S);
		family.kinship().bending();
	}

	/**
		Performs the quantile normalisation (inverse normal transformation) 
		on the phenotype or methylation values
		
		@param which which values normalise 
		@throws RuntimeException if the transformation can't be performed
	*/
	public void inverseNormalTransform(String which) throws RuntimeException 
	{		
		if (which.equals("response"))
			inverseNormalTransformPheno();
		else 
			inverseNormalTransformMeth();
	}
	
	/**
		Performs the quantile normalisation (inverse normal transformation) 
		on the phenotype values
		 
		@throws RuntimeException if the transformation can't be performed
	*/
	private void inverseNormalTransformPheno() throws RuntimeException 
	{
		for (int p=0; p<phenotypes.length; p++)
		{
			//extracts the valid phenotypes
			Vector<Double> v = new Vector<>();
			for (int i = 0; i < phenotypes[p].length; i++)
				if (phenotypes[p][i] != Utilities.INVALID_D)
				    v.add(phenotypes[p][i]);
			
			//Creates the rank vector
		    RankVector rankVector = new RankVector(v.toArray(new Double[]{}));
			rankVector.rank();
						
			//Does the transformation
			double[] zs = rankVector.transform();
		
			//conts the position inside zs
			int c = 0;
			for (int i = 0; i < phenotypes[p].length; i++)
				if (phenotypes[p][i] != Utilities.INVALID_D)
				{
					phenotypes[p][i] = zs[c];
					c++;
				}					
		}
	}
	
	/**
		Performs the quantile normalisation (inverse normal transformation) 
		on the methylation values
		 
		@throws RuntimeException if the transformation can't be performed
	*/
	private void inverseNormalTransformMeth() throws RuntimeException 
	{
		for (int p=0; p<methylations.length; p++)
		{
			//extracts the valid phenotypes
			Vector<Double> v = new Vector<>();
			for (int i = 0; i < methylations[p].length; i++)
				if (methylations[p][i] != Utilities.INVALID_D)
				    v.add(methylations[p][i]);
			
			//Creates the rank vector
		    RankVector rankVector = new RankVector(v.toArray(new Double[]{}));
			rankVector.rank();
						
			//Does the transformation
			double[] zs = rankVector.transform();
		
			//conts the position inside zs
			int c = 0;
			for (int i = 0; i < methylations[p].length; i++)
				if (methylations[p][i] != Utilities.INVALID_D)
				{
					methylations[p][i] = zs[c];
					c++;
				}					
		}
	}
	
	
	/**
		Corrects the methylation values either by the given 
		covariates (if Constants.corrects describes a file path)
		or by a number of Principal Components (PCs) sufficient to account
		for a given overall variability threshold (if Constants.corrects
		is a double value).
		
		In the former case the returned value is the number of read
		covariates, while in the latter case represents the number of 
		PCs used for the correction.
		
		If a covariate file is used the values for correction have already
		been loaded by the MyFileReader object.
		
		@precondition Fileloader loaded the correction values from file, if necessary
		
		@return the number of variable used for correction
		@see com.github.alesssia.poppante.MyFileReader
	*/
	public int correct()
	{
		//Decides which correction approach use according to the
		//Constants.correct value. If the user decided to use her own
		//covariates, they have been already loaded in the individuals'
		//data structure. If the user decided to correct the methylation
		//values using the PCs, then the PCA should be performed, the
		//correct number of PCs evaluated, and the selected values
		//loaded in the individuals's data structure.
		//The number of methylation covariates is the number of row in
		//the methylation covariates  matrix
		if (Utilities.isDouble(Constants.correct))
			calculateMethCovariates();
		
		int numCorrectionCovs = correctionCovs.length;
		
		//for each sites
		for (int m=0; m<listMeths.size(); m++)
		{
			//initialises the data structure, that are vectors instead
			//of arrays because I don't know how many people I have
			//with missing values
			Vector<Double> tmpY = new Vector<>();	
			Vector<Vector<Double>> tmpX = new Vector<>();
			for (int j = 0; j<numCorrectionCovs; j++)
				tmpX.add(new Vector());
					
			//extracts the valid methylations values and the
			//correction covariates for that person
			for (int i = 0; i < methylations[m].length; i++)
				if (methylations[m][i] != Utilities.INVALID_D)
				{
				    tmpY.add(methylations[m][i]);
					for (int j = 0; j<numCorrectionCovs; j++)
						tmpX.get(j).add(correctionCovs[j][i]);
				}
		
			//creates the Y vector with the read methylation values
			double[] data = new double[tmpY.size()+tmpY.size()*numCorrectionCovs];
			for (int i=0; i<tmpY.size(); i++)
			{
				data[i*(numCorrectionCovs+1)] = tmpY.get(i);
				for(int j=0; j<numCorrectionCovs; j++)
					data[i*(numCorrectionCovs+1)+j+1] = tmpX.get(j).get(i);
			}
						
			//Does the correction
			MyLinearRegression regressor = new MyLinearRegression(data, tmpY.size(), numCorrectionCovs); 
			double[] zs = regressor.residuals();

			//Sets the value back
			//c conts the position inside zs
			int c = 0;
			for (int i = 0; i < methylations[m].length; i++)
				if (methylations[m][i] != Utilities.INVALID_D)
				{
					methylations[m][i] = zs[c];
					c++;
				}		
		}

		//once the data structure have been populated the correction is
		//performed and then the data structure cleaned to free space.
		correctionCovs = null;
		return numCorrectionCovs;
	}
	
	
	/**
		Evaluates the PCA of the methylation sites, evaluates the 
		number of PCs necessary to account for the total variability,
		and extract them. 
		
		It substitutes the missing values with zero. However, a matrix
		line (that is one of the PCA variables) composed only by zero 
		changes the PCA results. It is then necessary to remove all the
		mock individuals (that are those having only missing values) before
		calling this function.
		
		@precondition the correction should be performed by means of PCs
		@precondition removeMock should have been called.
	*/
	private void calculateMethCovariates()
	{
		assert Utilities.isDouble(Constants.correct) : "Internal error : --correct requires a number.";
	
		//Variables should stay in the columns (I already creted the matrix
		//that way)
		MyPCA pca = new MyPCA(methylations, Utilities.INVALID_D);
		pca.evaluatePCA();
		
	    // selects the number of PC to use and initialises
		//the data structure accordingly 
		pca.evaluateProportionOfVariance();
		int numMethCovar = pca.howMany(Double.parseDouble(Constants.correct));
		resetCorrectionCovs(numMethCovar);
		
	    //loads correction values (they need to be transposed)
		Matrix tmp = new Matrix(pca.getPCs(numMethCovar));
		correctionCovs = tmp.transpose().getArray();
	}
	
	/**
		Sets the missingness code for each methylation site.
		
		Each site is assigned to a missingness pattern, that
		is a unique code identifying the people having
		missing values, as recorded when the methylation data 
		were read.
		
		This code will be used to identify group of poeple/
		methylation sites to be analysed together (that is:
		using the same null model)
		
		@precondition the request position is valid
		@see com.github.alesssia.poppante.Marker
		@see com.github.alesssia.poppante.MyFileReader.readMethylationData
	*/
	public void setMissingnessMethPattern()
	{
		assert Constants.mode == Utilities.MODE_ASSOCIATION : "Internal error: missingness pattern should not be set when predictor heritability is assessed.";
		assert listMeths.size() > 0 : "Internal error: no predictor information available";
		
		missingnessMethPattern = new String[listMeths.size()];
		for (int i=0; i<listMeths.size(); i++)
		{
			Marker marker = listMeths.get(i);
			missingnessMethPattern[i] = hashids.encode(marker.missing());
			//once the code is set, the list of person can be deleted
			//from the marker, in order to save space
			marker.resetMissing();
		}		
	}
	
	/**
		Sets the missingness code for each phenotype.
		
		Each phenotype is assigned to a missingness pattern, that
		is a unique code identifying the people having
		missing values, as recorded when the phenotype data 
		were loaded in the DataManager structures.
		This code will be used to identify group of poeple/
		methylation sites to be analysed together (that is:
		using the same null model).
		
		This function also initialises a null model for each
		phenotype, that will be used if there is a methylation
		site sharing the same missingness pattern.
		
		@precondition the request position is valid
		@see com.github.alesssia.poppante.Phenotypes
	*/
	public void setMissingnessPhenoPattern()
	{
		assert Constants.mode == Utilities.MODE_ASSOCIATION : "Internal error: missingness pattern should not be set when predictor heritability is assessed.";	
		assert phenotypeNames.size() > 0 : "Internal error: no response variable available";
		
		for (int i=0; i<phenotypeNames.size(); i++)
		{
			Phenotype pheno = phenotypeNames.get(i);
			//the pattern is set in order to compare it with the 
			//methylation site's one
			//No missing values are an empty string
			missingnessPattern.put(i, new NullModel(hashids.encode(pheno.missing())));

			//once the code is set, the list of person can be deleted
			//from the marker, in order to save space
			pheno.resetMissing();
		}		
	}
	
	
	/**
		Returns the null model associated to the missingness pattern
		of the given phenotype.
		
		It is the VC that decide whether used it or not
		
		@param phenotypePosition the phenotype for which we want to extract the null model
		@return the null model
	*/	
	public NullModel getNull(int phenotypePosition)
	{
		assert missingnessMethPattern.length > 0 : "Internal error: missingness pattern not initialised.";
		return missingnessPattern.get(phenotypePosition);
	}
	
	/**
		Returns the missingness pattern of the given methylation site.
		
		@param methylationPosition the methylation site for which we want to know the pattern
		@return the missingness pattern (code)
	*/
	public String missingnessMethPattern(int methylationPosition)
	{
		assert Constants.mode == Utilities.MODE_ASSOCIATION : "Internal error: missingness pattern should not be required when predictor heritability is assessed.";
		
		return missingnessMethPattern[methylationPosition];
	}
}
