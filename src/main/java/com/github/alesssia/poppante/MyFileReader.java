/*
* 	 MyFileReader.java
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
import java.io.*;


/**
	Reads the information from the input files.

	@author      Alessia Visconti <alessia.visconti@kcl.ac.uk>
	@version     1.0  
*/

class MyFileReader 
{
	/** This is where the data are stored */
	private DataManager datamanager;
	
	/** Position of the methylation sites whose value should
		be read. */
	private Vector<Integer> position2Read;
	/** Position of the phenotypes whose value should
		be read. */
	private Vector<Integer> position2ReadP;
	
	/** Total number of methylation sites that should be
		recorded in the METH file */
	private int totMethSite;
	/** Total number of phenotypes that should be
		recorded in the PHENO file */
	private int totPhenos;
	
		
	
	/**
		Constructor.
		
		Initalise the data structure.
		@param dm the DataManager object
	*/
	public MyFileReader(DataManager dm)
	{
		datamanager = dm;
		position2Read = new Vector<>();
		position2ReadP = new Vector<>();
		totMethSite = Utilities.INVALID_I;
		totPhenos = Utilities.INVALID_I;
	}
	
	
	/**
		Reads the family structure and the individuals' pedigree 
		information along with the phenotype values. 
		
		It also verifies the validity of the pedigree file. 
		When the mode "pedcheck" is selected instead of throwing an
		exception the errors in the pedigree file are saved in a 
		string and returned in the FamilySet object.
		
		@precondition the datamanager object should know the number of phenotypes
		to be read when the mode is not PEDCHECK
		
		@throws Exception if the phenotypes list is not available when mode is not PEDCHECK
		@throws IOException if the pedigree file can't be read.
		@see com.github.alesssia.poppante.FamilySet
		@return the family information
	*/
	public FamilySet readFamilies()  throws Exception, IOException
	{	
		if (Constants.mode == Utilities.MODE_ASSOCIATION && datamanager.phenotypeNames().size() == 0)
			throw new Exception("Internal Error: response list not initialised.");
			
		int people = 0;
		String message = "";
		
		//creates a fake family for unrelated people analisys
		if (Constants.kinship != null)
			datamanager.families().put(Utilities.INVALID_S, new Family(Utilities.INVALID_S));

		try 
		{	
			RandomAccessFile file = new RandomAccessFile(Constants.ped, "r");
			String s;
			while ((s = file.readLine()) != null) 
			{
				people++;
				Person person = readIndividual(s, people);
	
				//assign the person to her own family
				if (Constants.kinship == null)
				{
					String personFamily = person.famID();
					if (datamanager.families().get(personFamily) == null)
						datamanager.families().put(personFamily, new Family(personFamily));
						
					datamanager.families().get(personFamily).addMember(person);
				}
				//all the individuals are supposed to be unrealed an thus assigned to 
				//the same family (with id == INVALID_I), and the individual IDs 
				//becames the concatenation of famID and ID
				else 
				{
					//resets the IDs to have a mock values
					person.setID(new String(person.famID() + "" + person.id())); 	
					
					//founders remain founders (that is their parent's ID will be always zero)
					if (!person.isFounder())
					{
						person.setMotherID(new String(person.famID() + "" + person.motherID())); 	
						person.setFatherID(new String(person.famID() + "" + person.fatherID())); 	
					}
					person.setfamID(Utilities.INVALID_S);				
					
					datamanager.families().get(Utilities.INVALID_S).addMember(person);
				}	
			}
		}
		catch(Exception e)
		{
			if (Constants.mode == Utilities.MODE_PEDCHECK)
				message += "\t- " + e.getMessage() + "\n";
			else
				throw new IOException("ERROR: " + e.getMessage());
		}
		
		return new FamilySet(people, message); 
	}
	
	/**
		Parses the information of an individual and checks errors
		in the format.
		
	
		@precondition the datamanager object should know the number of phenotypes
		to be read when the mode is not PED CHECK
	
		@param s the string codifying a person
		@param counter the position of the string in the file
		@return the new Person	
		@throws NotWellFormedLineException if the file has a line that is not well-formed
		@see com.github.alesssia.poppante.Person
		@see com.github.alesssia.poppante.Utilities.MISSING_VALUES
	*/
	private Person readIndividual(String s, int counter) throws NotWellFormedLineException
	{
		int INFO_SIZE = 7;	
		String famID, id, fatherId, motherId;
		int sex, affection, twin;
		double[] phenotypes = new double[datamanager.phenotypeNames().size()];
		boolean isMock = false;
		
		try 
		{
			String[] info = s.split("\\s+");
			
			//the line is not valid: some information is missing
			if (info.length < INFO_SIZE) //affections and twins aren't int!
				throw new NotWellFormedLineException("Not well-formed PED file at line " + counter);
			
			//when no value is specified a missing value is used
			sex = 0;
			//Lets check and transforms some values.
			if (info[4].equals("1") || info[4].equals("M"))
				sex = 1;
			else if (info[4].equals("2") || info[4].equals("F"))
				sex = 2;
			else if (info[4].equals("0") || Utilities.isMissing(info[4]))
				sex = 0;
			//This error is reported only if the PED CHECK mode is selected,
			//because it is not used into the analysis
			else if (Constants.mode == Utilities.MODE_PEDCHECK)
				throw new NotWellFormedLineException("Sex information is not valid at line " + counter);
			
			//Affection is used for case/control studies.
			//It must be [0,1]. Missing values are not modified
			affection = Utilities.MISSING;
			if (info[5].equals("0"))
				affection = 0;
			else if (info[5].equals("1"))
				affection = 1;
			else if (info[5].equals("-9"))
				isMock = true;
			
			//This error is reported only if the PED CHECK mode is selected,
			//because it is not used into the analysis
			else if (Constants.mode == Utilities.MODE_PEDCHECK)
				throw new NotWellFormedLineException("Affection information is not valid at line " + counter);
			
			
			if (info[6].equals("1") || info[6].equals("MZ"))
				twin = 1;
			else if (info[6].equals("2") || info[6].equals("DZ"))
				twin = 2;
			else if (info[6].equals("0"))
				twin = 0;
			else
				throw new NotWellFormedLineException("Twin information is not valid at line " + counter);
                        
			//real IDs are read, they will be changed afterwards if necessary
			famID 	 = info[0]; 
			id 		 = info[1];
			fatherId = info[2]; 
			motherId = info[3]; 
			
			
			//phenotype values are loaded only if the analysis mode is Association
			if (Constants.mode == Utilities.MODE_ASSOCIATION)
			{
				if (info.length != INFO_SIZE + totPhenos)
					throw new NotWellFormedLineException("The number of response variables does not match the expected value at line " + counter);
			
				//Continue t parsing,otherwise
				int skippedPheno = 0;
				for(int i=INFO_SIZE; i<info.length; i++) 
				{
					int phenoPosition = i-INFO_SIZE;
					
					//if this position is not in the list of included
					//position is not read
					if (!position2ReadP.contains(phenoPosition))
					{
						skippedPheno++;
						continue;
					}
					
					//since not all the position are read, I can't simply use:
					//	int methPosition = i-INFO_SIZE;
					//to identify the position of the methylation value in the methylation
					//matrix, but I also need to remove those that have been skipped.
					phenoPosition -= skippedPheno;
															
					String phenotype = info[i];
					try 
					{
						if (Utilities.isMissing(phenotype)) 
							phenotypes[phenoPosition] = Utilities.INVALID_D;
						else
							phenotypes[phenoPosition] = Double.parseDouble(phenotype); 
					} 
					catch (Exception exc) 
					{
						throw new NotWellFormedLineException("Wrong response values for individual " + id + "\n\tCheck PED file at line " + counter);
					}
				}
			}
			
			return new Person(famID, id, fatherId, motherId, sex, affection, twin, phenotypes, isMock);
		}
		catch (NumberFormatException e)
		{
			throw new NotWellFormedLineException("Error at line " + counter + " of the PED file. The format is not correct.");
		}
		catch (NotWellFormedLineException exception)
		{
			throw exception;
		}	
	}
	
	
	
	
	/**
		Reads the map of the methylation sites and returns the
		number of sites.
		
		The list of methylation sites is stored in the datamanager object.
		
		If the file format is not correct (e.g., the number of read 
		fields is smaller than three) it throws an exception specifying 
		which row raises the problem. A format problem is usually
		a not numeric values for one of the data (methylation site
		coordinate) 

		@return the number of loaded methylation sites.
		@throws IOException if the methylation file can't be read
		@throws NotWellFormedLineException if the file has a line that is not well-formed
		
		@see com.github.alesssia.poppante.Marker
	*/
	public int readMethylationInformation() throws IOException, NotWellFormedLineException
	{
		//vector of the position to read, according to the list of sites to
		//include in the analysis
		try 
		{
			RandomAccessFile file = new RandomAccessFile(Constants.map, "r");

			totMethSite = 0;
			
			String s = file.readLine();
			String[] info = s.split("\\s+");
			
			//The first line decided whether the information includes
			//only the name or name and position
			int INFO_SIZE = (info.length >= 3) ? 3 : 1;
			datamanager.setMethFields(INFO_SIZE);
			
			do
			{	
				totMethSite++;
				info = s.split("\\s+");
				
				//If we are supposed to use only the name but some marker 
				//has also position is fine, otherwise it will throw an excpetion             
				if (info.length < INFO_SIZE)
					throw new NotWellFormedLineException("ERROR: line " + totMethSite + " is not a valid predictor description: it should include " + INFO_SIZE + "columns. Please check your MAP file");		
				//it is loaded only if the used asked for it to be analysed
				else if (Constants.include == null || datamanager.includedMeth().contains(info[0]))
				{
					//totMethSite counts the lines, that starts from 1, but the 
					//positions starts from 0 
					position2Read.add((totMethSite-1));
					if(INFO_SIZE == 3)
						datamanager.listMeths().add(new Marker(info[0], info[1], Integer.parseInt(info[2])));
					else 
						datamanager.listMeths().add(new Marker(info[0]));
				}
			} while ((s = file.readLine()) != null);

			file.close();
			
			return datamanager.listMeths().size();	
		}
		catch (NumberFormatException e)
		{
			throw new NotWellFormedLineException("ERROR: line " + totMethSite + " does not describe a valid marker site (not numeric position). Please check your MAP file");
		}
		catch(IOException e) 
		{
			throw new IOException("ERROR: Check --map parameter or MAP file");
		}	
	}

	/**
		Reads the phenotype names and returns the number of 
		phenotypes that have been loaded. 
		
		The list of phenotypes is stored in the datamanager object.
		
		It reads only the phenotype names (that is the first 
		column of the file), whilst the other information are 
		discarded.
		If an error occurs an exception is thrown. The exception 
		messages specifies the line showing the issues. 
		A format problem is usually a not numeric values for one
		of the data.
	
		@return the number of phenotypes in the dataset
		@throws IOException if the phenotype file can't be read
		@throws NotWellFormedLineException if the file has a line that is not well-formed
		@see com.github.alesssia.poppante.Phenotype
	*/
	public int readPhenotypesInformation() throws IOException, NotWellFormedLineException
	{
		try 
		{
			RandomAccessFile file = new RandomAccessFile(Constants.response, "r");
			totPhenos = 0;
			
			String s = file.readLine();
			String[] info = s.split("\\s+");
						
			do
			{	
				totPhenos++;
				info = s.split("\\s+");
				if (info.length < 1)
					throw new NotWellFormedLineException("ERROR: line " + totPhenos + " does not describe a valid response variables. Please check your RESPONSE file");
				
				//it is loaded only if the used asked for it to be analysed
				else if (Constants.filter == null || datamanager.includedPheno().contains(info[0]))
				{
					//totPhenos counts the lines, that starts from 1, but the 
					//positions starts from 0 
					position2ReadP.add((totPhenos-1));
					datamanager.phenotypeNames().add(new Phenotype(info[0]));
				}
			} while ((s = file.readLine()) != null);
			
			file.close();
			
			return datamanager.phenotypeNames().size();	
		} 
		catch (IOException exception) 
		{
			throw new IOException("ERROR: Check --response parameter or RESPONSE file.");
		}		
	}

	/**
		Reads the first line of the covariate file and sets the 
		number of covariate in the datamanager object.
		
		This method assumes that every line has the same number of 
		covariates (as the readCovariates method does as well).

		@throws IOException If the covariate file can't be read
		@throws NotWellFormedLineException If the file has a line that is not well-formed
	*/
	public void setNumCovar() throws IOException, NotWellFormedLineException
	{
		try 
		{
			RandomAccessFile file = new RandomAccessFile(Constants.covariate, "r");
			int INFO_SIZE = 2;
			String s = file.readLine();
			file.close();
			
			if (s != null)
			{
				String[] info = s.split("\\s+");
				//family and individual IDs are mandatory, and at least a covariate should be available
				if (info.length < INFO_SIZE) 
					throw new NotWellFormedLineException("ERROR: covariate file does not describe valid covariate values. Please check your COVARIATE file");
				datamanager.setNumCovar(info.length - INFO_SIZE);
			}
			else
				throw new NotWellFormedLineException("ERROR: covariate file does not describe a valid covariate values. Please check your COVARIATE file");
		}
		catch (IOException  exception) 
		{
			throw new IOException("ERROR: Check --covariate parameter or COVARIATE file.");
		}		
	}
	
	/**
		Reads the first line of the correction file and gets the 
		number of correction covariates, used by datamanager to 
		initialise the data structure.
		
		This method assumes that every line has the same number of 
		covariates (as the readCorrectionCovariates() method does as well).

		@throws IOException If the covariate file can't be read
		@throws NotWellFormedLineException If the file has a line that is not well-formed
		
		@return the number of correction covariates
	*/
	public int getNumCorrectionCov() throws IOException, NotWellFormedLineException
	{
		try 
		{
			RandomAccessFile file = new RandomAccessFile(Constants.correct, "r");
			int INFO_SIZE = 2;
			String s = file.readLine();
			file.close();
			
			if (s != null)
			{
				String[] info = s.split("\\s+");
				//family and individual IDs are mandatory, and at least a covariate should be available
				if (info.length < INFO_SIZE) 
					throw new NotWellFormedLineException("ERROR: COVARIATE file does not describe valid covariate values for predictors correction. Please check your file");
				return (info.length - INFO_SIZE);
			}
			else
				throw new NotWellFormedLineException("ERROR: COVARIATE file does not describe valid covariate values for predictors correction. Please check your file");
		}
		catch (IOException  exception) 
		{
			throw new IOException("ERROR: Check --correct parameter or COVARIATE file.");
		}		
	}


	/**
		Reads the covariate values.
		
		It verifies the the file is well-formed, that is the family and 
		the subject ID are specified. 
		If a set of covariates refers to a subject that is not listed in
		the pedigree file the information is discarded. 
		A format problem is usually	a not numeric values for one of 
		the data.
		
		@precodintion PED file should have been read

		@throws IOException If the covariate file can't be read
		@throws NotWellFormedLineException If the file has a line that is not well-formed
		@see com.github.alesssia.poppante.datamanager
		@see com.github.alesssia.poppante.Utilities.MISSING_VALUES
	*/
	public void readCovariates() throws IOException, NotWellFormedLineException
	{	
		assert datamanager.numFamilies() != 0 : "Internal error. No family loaded.";
			
		int counter = 0;
		try 
		{
			RandomAccessFile file = new RandomAccessFile(Constants.covariate, "r");
			int INFO_SIZE = 2;
			String s;
	
			while ((s = file.readLine()) != null) 
			{
				counter++;
				String[] info = s.split("\\s+");
		
				//family and individual IDs are mandatory, and at least a covariate should be available
				if (info.length < INFO_SIZE) 
					throw new NotWellFormedLineException("ERROR: line " + counter + " does not describe a valid COVARIATE. Please check your COVARIATE file");
			
				//The check of the correct number of covariates is done only for
				//individual for whom we have information in the PED file			
				if (info.length - INFO_SIZE != datamanager.numCovar())
					throw new NotWellFormedLineException("ERROR: line " + counter + " has an unexpected number of covariates.");
					
				//gets the correct IDs	
				Family family;
				if (Constants.kinship == null)
					family = datamanager.families().get(info[0]);
				else
				{
					family = datamanager.families().get(Utilities.INVALID_S);
					info[1] = new String(info[0] + "" + info[1]);
				}
			
				//this family does not exist
				if (family == null)
					continue;
				
				//this person does not exist
				Person person = family.getMemberByID(info[1]);
				if (person == null)
					continue;
								
				//read the actual covariates and check whether there are missing 
				//covariate values
				double[] covariates = new double[datamanager.numCovar()];
				for (int i=INFO_SIZE; i<info.length; i++)
				{
					int covPosition = i-INFO_SIZE;
										
					// if a covariate is missing the individual is set to have missing 
					// covariates, the covariate read so far removed and the line parsing stopped
					if (Utilities.isMissing(info[i]))
					{
						covariates = null;
						break;
					}
					else 
						covariates[covPosition] = Double.parseDouble(info[i]);
				}
			
				person.setCovariates(covariates);

			}
	
			file.close();
		}
		catch (NumberFormatException e)
		{
			throw new NotWellFormedLineException("ERROR: line " + counter + " does not describe a valid covariates (not numeric value). Please check your COVARIATE data file");
		}
		catch (IOException  exception) 
		{
			throw new IOException("ERROR: Check --covariate parameter or COVARIATE file.");
		}	
	}

	/**
		Reads the methylation value data.
		
		It also verifies the the file is well-formed, that is the 
		family and the subject ID are specified. If a set of 
		methylation values refers to a subject that is not listed
		in the pedigree file the information is discarded.
		A format problem is usually a not numeric values for one 
		of the data.
	
		@precodintion the positionTable should have been initialised
		@precodintion methylation information should have been read
		@precodintion the methylation table should have been initialised
	
		@throws IOException If the methylation data file can't be read
		@throws NotWellFormedLineException if the file has a line that is not well-formed
		
		@see com.github.alesssia.poppante.Utilities.MISSING_VALUES
	*/
	public void readMethylationData() throws IOException, NotWellFormedLineException
	{
		assert datamanager.people() != 0 : "Internal error: position table has not ben initialised.";
		assert datamanager.listMeths().size() != 0 : "Internal error: no predictor information."; 
		assert datamanager.isMethylationsInitialised() : "Internal error: predictor table not initialised";
		
		int counter = 0;
		try 
		{
			int methNumber = datamanager.listMeths().size();
			RandomAccessFile file = new RandomAccessFile(Constants.predictor, "r");
		
			int INFO_SIZE = 2;
			String s;
		
			while ((s = file.readLine()) != null) 
			{
				counter++;
				String[] info = s.split("\\s+");
		
				//family and individual IDs are mandatory, and the number of read methylatin sites
				//should match the expected one
				if (info.length != INFO_SIZE + totMethSite)
					throw new NotWellFormedLineException("The number of predictor values does not match the expected value at line " + counter);
			
				//Extract the position of that individual in the methylation matrix
				String key = "";
				if (Constants.kinship == null)
					key = new String(info[0] + "" + info[1]);
				else
					key = new String(Utilities.INVALID_S + "" + info[0] + "" + info[1]);
				
				Integer personPosition = datamanager.getPosition(key);
				
				//If the individual is not present in the dataset (PED file) 
				//the methylation data will be discharged 
				if (personPosition == null)
					continue;
				
				int skippedMeth = 0;
				//read the actual values
				for (int i=INFO_SIZE; i<info.length; i++)
				{
					int methPosition = i-INFO_SIZE;
					
					//if this position is not in the list of included
					//position is not read
					if (!position2Read.contains(methPosition))
					{
						skippedMeth++;
						continue;
					}
					
					//since not all the position are read, I can't simply use:
					//	int methPosition = i-INFO_SIZE;
					//to identify the position of the methylation value in the methylation
					//matrix, but I also need to remove those that have been skipped.
					methPosition -= skippedMeth;

					// add to each marker information the position of the person
					// for whom a value is missing. This list of person will be
					// used in the main program to identify missingness pattern,
					// that is group of people with the same missing values.
					if (Utilities.isMissing(info[i]))
					{
						datamanager.setMeth(methPosition, personPosition, Utilities.INVALID_D);
						datamanager.listMeths().get(methPosition).addMissing((long)personPosition);
					}
					else 
						datamanager.setMeth(methPosition, personPosition, Double.parseDouble(info[i]));
				}
				
			}
	
			file.close();		
		
		}
		catch (NumberFormatException e)
		{
			throw new NotWellFormedLineException("ERROR: line " + counter + " does not describe a valid predictor (not numeric value). Please check your PREDICTOR data file");
		}
		catch (IOException  exception) 
		{
			throw new IOException("ERROR: Check --predictor argument or PREDICTOR data file");
		}
		
		
	}
	
	/**
		Reads the kinship data.
		
		It also verifies whether the the file is well-formed,
		If a kinship value refers to a subject that is not listed
		in the pedigree file the information is discarded.
		
		@precondition a mock family should be created 
		@precodintion the kinship matrix initialised beforehand
		@precodintion the positionTable should have been initialised
		
		@throws Exception if the mock family has not been created, or the kinship/positionTable not initialised
		@throws IOException if the kinship file can't be read
		@throws NotWellFormedLineException if the file has a line that is not well-formed
		
		@see com.github.alesssia.poppante.Kinship
	*/
	public void readKinship() throws Exception, IOException, NotWellFormedLineException
	{
		//This function does specify the family since only one
		//family (with ID == Utilities.INVALID_S)) is stored when there
		//are unrelated individuals
		Family family = datamanager.families().get(Utilities.INVALID_S);
		if (family == null)
			throw new Exception("Internal error: mock family for unrelated individuals has not ben initialised."); 
		
		if (family.kinship().matrix() == null)
			throw new Exception("Internal error: kinship matrix for unrelated individuals has not ben initialised."); 
        
		if (datamanager.people() == 0)
			throw new Exception("Internal error: position table has not ben initialised."); 
		
		int counter = 0;   
		try
		{
			int INFO_SIZE = 5;

			RandomAccessFile file = new RandomAccessFile(Constants.kinship, "r");
			String s;	
			while ((s = file.readLine()) != null)
			{
				counter++;
				String[] info = s.trim().split("\\s+");

				if (info.length < INFO_SIZE)
					throw new NotWellFormedLineException("ERROR: line " + counter + " does not describe a valid kinship entry. Please check your KINSHIP file");

				//If one of the individuals is not present in the dataset (PED file)
				//their kinship information will be discharged
				//I am looking into the datamanager map system, that uses the 
				//mock family ID to get the position into the array, so the real family ID
				//(info [0] and [2], should be not considered.
				Integer p1 = datamanager.getPosition(new String(Utilities.INVALID_S + "" + info[0] + "" + info[1]));
				if (p1 == null)
					continue;
				
				Integer p2 = datamanager.getPosition(new String(Utilities.INVALID_S + "" + info[2] + "" + info[3]));
				if (p2 == null)
					continue;
				
				double genomicKinship = Double.parseDouble(info[4]); 
				if ( genomicKinship < Constants.mink) 
					genomicKinship = 0.0;   
				family.kinship().setValue(genomicKinship, p1, p2);	
			}
		
			file.close();			
		}
		catch(NumberFormatException exception)
		{
			throw new NotWellFormedLineException("ERROR: line " + counter + " does not describe a valid kinship entry (not numeric value). Please check your KINSHIP data file");
		}
		catch (IOException  exception)
		{
			throw new IOException("ERROR: Check --kinship argument or KINSHIP file.");
		}
	}
	
	/**
		Reads the list of the methylation sites to be analysed and
		returns their number.
		
		It also verifies whether the the file is well-formed,
		
		@return the number of methylation sites to be analysed.
		@throws IOException if the methylation file file can't be read
		@throws NotWellFormedLineException if the file has a line that is not well-formed
		@see com.github.alesssia.poppante.DataManager
	*/
	public int readIncludedMeths() throws IOException
	{
		try
		{
			RandomAccessFile file = new RandomAccessFile(Constants.include, "r");
			
			int counter = 0;
			String s;
	        
			while ((s = file.readLine()) != null)
			{
				counter++;
				String[] info = s.trim().split("\\s+");
	
				if (info.length < 1)
					throw new NotWellFormedLineException("ERROR: line " + counter + " does not describe a valid predictor name. Please check your INCLUDE file");
			
				datamanager.includedMeth().add( info[0] );
			}
	            
			return datamanager.includedMeth().size();			
		}
		catch (IOException exception)
		{
			throw new IOException("ERROR: Check --include argument or INCLUDE file.");
		}
	}

	/**
		Reads the list of the phenotypes to be analysed and returns
		their number.
		
		It also verifies whether the the file is well-formed,
		
		@return the number of phenotypes to be analysed.
		@throws IOException if the phenotypes file file can't be read
		@throws NotWellFormedLineException if the file has a line that is not well-formed
		@see com.github.alesssia.poppante.DataManager
	*/
	public int readIncludedPhenos() throws IOException
	{
		try
		{
			RandomAccessFile file = new RandomAccessFile(Constants.filter, "r");
			
			int counter = 0;
			String s;
	        
			while ((s = file.readLine()) != null)
			{
				counter++;
				String[] info = s.trim().split("\\s+");
	
				if (info.length < 1)
					throw new NotWellFormedLineException("ERROR: line " + counter + " does not describe a valid response name. Please check your FILTER file");
			
				datamanager.includedPheno().add( info[0] );
			}
	            
			return datamanager.includedPheno().size();			
		}
		catch (IOException exception)
		{
			throw new IOException("ERROR: Check --filter argument or FILTER file.");
		}
	}
	
	
	/**
		Reads the covariates to use for correcting the methylation values 
		and returns their number. 
		
		It verifies the the file is well-formed, that is the family and 
		the subject ID are specified. 
		If a set of covariates refers to a subject that is not listed in
		the pedigree file the information is discarded. 
		A format problem is usually	a not numeric values for one of 
		the data.

		@precondition  the file should have no missing value 
		
		@throws Exception If the covariate file includes a missing value
		@throws IOException If the covariate file can't be read
		@throws NotWellFormedLineException If the file has a line that is not well-formed
		@see com.github.alesssia.poppante.Person
		@see com.github.alesssia.poppante.Utilities.MISSING_VALUES
	*/
	public void readCorrectionCovariates() throws Exception, IOException, NotWellFormedLineException
	{
		assert datamanager.people() != 0 : "Internal error: position table has not ben initialised.";
		assert datamanager.isCorrectionCovsInitialised() : "Internal error: correction table not initialised";
		
		//number of covariates read in the previous line
		int nCor = -1;
		try
		{
			RandomAccessFile file = new RandomAccessFile(Constants.correct, "r");
			int INFO_SIZE = 2;
			String s;

			int counter = 0;
			while ((s = file.readLine()) != null)
			{
				counter++;
				String[] info = s.split("\\s+");

				//family and individual IDs are mandatory, and at least a covariate should be available
				if (info.length < INFO_SIZE)
					throw new NotWellFormedLineException("ERROR: line " + counter + " does not describe a valid covariate. Please check your CORRECTION file");

				///Extract the position of that individual in the methylation matrix
				String key = "";
				if (Constants.kinship == null)
					key = new String(info[0] + "" + info[1]);
				else
					key = new String(Utilities.INVALID_S + "" + info[0] + "" + info[1]);
				
				Integer personPosition = datamanager.getPosition(key);
				
				//If the individual is not present in the dataset (PED file) 
				//the methylation data will be discharged 
				if (personPosition == null)
					continue;
				
				//The check of the correct number of covariates is done only for
				//individual for whom we have information in the PED file
				if (nCor > 0 && nCor != info.length - INFO_SIZE)
					throw new NotWellFormedLineException("ERROR: line " + counter + " has a number of covariates that does not match the previous lines.");

				nCor = info.length - INFO_SIZE;
				
				//read the actual values
				for (int i=INFO_SIZE; i<info.length; i++)
				{
					// if methvalue is missing set affection to missing
					if (Utilities.isMissing(info[i]))
						throw new Exception("ERROR: The -correct file includes a missing value at line " + counter + "\nPlease use a file with no missing value.");
					else 
						datamanager.setCor(i-INFO_SIZE, personPosition, Double.parseDouble(info[i]));
				}
			}
		}
		catch (IOException exception)
		{
			throw new IOException("ERROR: Check --correct parameter or the CORRECTION file.");
		}
	}
	
	
	
}


/**
	Runtime exception. 

	Raised when a line of the input file is not well formed.

	@author      Alessia Visconti <alessia.visconti@kcl.ac.uk>
	@version     1.0               
*/
class NotWellFormedLineException extends RuntimeException
{
	/**  serial number */
	public static final long serialVersionUID = 4242L; 

	/**
		Constructor. 
		
		Initialises the exception message.
		
		@param msg the exception message
	*/
	public NotWellFormedLineException(String msg) 
	{
		super(msg);
	}
}