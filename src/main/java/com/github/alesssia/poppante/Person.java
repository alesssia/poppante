/*
 * 	 Person.java
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
	Represents a person. 

	It includes information about the person herself (such as her 
	sex), the population (such as the family the person belong to) 
	and the structure of the family the person belong to (such has 
	the parents' IDs).

	When the family structure is available the stored IDs are the 
	real ones, while when individuals are not related the stored IDs
	are mock ones.

	It includes also data structure to load phenotype and covariate 
	values, that are used to temporarily store these information. In
	fact, once the data are fully loaded, they will be copied in the 
	main class and removed.
	
	@author      Alessia Visconti <alessia.visconti@kcl.ac.uk>
	@author      Mario Falchi     <mario.falchi@kcl.ac.uk>
	@author      Hariklia Eleftherohorinou <hariklia.eleftherohorinou06@imperial.ac.uk>
	@version     1.0     
*/

class Person 
{
	/** Family's ID. */
    private String famID;   
	/** Individual's ID. */		
    private String id;		 	
	/** Father's ID. 
		
		This value is set to zero if the individual
		is a founder. */	
    private String fatherID; 
	/** Mother's ID. 
		
		This value is set to zero if the individual
		is a founder. */		
    private String motherID;		
	/** Individual's gender.
		
		The values should be in [0,2], 0 being missing, 
		1 male, and 2 female */
    private final int sex;		
	/** Affection status.
		
		This value is not used in the analysis. */	
    private final int affection;  	
	/** Twin status.
		
		The values should be in [0,2], 0 being not-twin/unknown, 
		1 monozygotic twin, and 2 dizygotic twin */
    private final int twin;		
		
	/** Values of this individual's phenotypes */
    private double[] phenotypes; 
	/** Whether is a mock individual. */		
    private final boolean isMock;	
	/** Values of this individual's covariates */
    private double[] covariates; 
	/**
		Constructor. 
		
		Describes the person by means of her pedigree information, and
		loads her phenotypes values. It also initialises the other data 
		structures.
		
		@precondition sex is represented by a valid value
		@precondition twin status is represented by a valid value
		@precondition affection status is represented by a valid value
		@param famID the family's ID
		@param id the individual's ID
		@param fatherID the individual's father's ID
		@param motherID the individual's mother's ID
		@param sex the individual's gender
		@param affection the individual's affection status
		@param twin the individual's twin status
		@param phenotypes the vector of phenotypes
		@param isMock whether it represents a mock individual
	*/
	public Person (String famID, String id, String fatherID, String motherID, int sex, int affection, int twin, double[] phenotypes, boolean isMock) 
    {
		assert sex >= 0 && sex <= 2 : "Not valid sex";
		assert twin >= 0 && twin <=2 : "Not valid twin code";
		assert affection == 0 || affection == 1 || affection == Utilities.MISSING : "Not valid affection status";
		
		this.famID = famID;
		this.id = id;
		this.fatherID = fatherID;
		this.motherID = motherID;
		this.sex = sex;
		this.affection = affection;
		this.twin = twin;
		this.phenotypes = phenotypes;
		this.isMock = isMock;
		
		//Since covariates is set to NULL all Person
		//missesCovariates == true
		covariates = null;
    }

	/**
		Returns the family ID.
		
		@return the family ID
	*/
    public String famID() 
    {
		return famID;
    }
	
	/**
		Returns the individual's ID.
		
		@return the individual's ID
	*/	
    public String id()
    {
		return id;
    }
	
	
	/**
		Sets the individual's ID.
		
		@param id the individual's ID
	*/
	public void setID(String id)
    {
		this.id = id;
    }
	
	
	/**
		Sets the individual's family ID.
		
		@param fid the individual's family ID
	*/	
    public void setfamID(String fid)
    {
		famID = fid;
    }
    
	/**
		Returns whether the individual is a monozygotic twin.
		
		@return true if the individual is a monozygotic twin, false otherwise
	*/    
    public boolean isMZ()
    {
		return twin == 1;
    }
	
	/**
		Returns whether the individual is founder. 
		
		@return true if the individual is founder, false otherwise
	*/		
    public boolean isFounder()
    {
		return fatherID.equals("0") && motherID.equals("0");
    }
	
	/**
		Returns whether the individual's gender and that
		given as parameter are the same.
		
		@param supposed the supposed gender
		@return true if the gender are the same, false otherwise
	*/		
    public boolean isValidParent(int supposed)
    {
		return sex == supposed;
    }
    
	/**
		Returns the individual's father's ID.
		
		@return the individual's father's ID
	*/    
    public String fatherID()
    {
		return fatherID;
    }
	
	/**
		Sets the individual's father's ID.
		
		@param id the new father's ID
	*/
	public void setFatherID(String id)
    {
		fatherID = id;
    }
	
    
	/**
		Returns the individual's mother's ID.
		
		@return the individual's mother's ID
	*/   
    public String motherID()
    {
		return motherID;
    }
	
	/**
		Sets the individual's mother's ID.
		
		@param id the new mother's ID
	*/
	public void setMotherID(String id)
    {
		motherID = id;
    }
	
	
	/**
		Returns the individual's gender as a numeric
		code.
		            
		@return the individual's gender
	*/		
    public int sexCode()
    {
		return sex;
    }
	
	/**
		Returns the individual's gender as a String
		            
		@return the individual's gender
	*/		
    public String sexString()
    {
		switch (sex)
		{
			case 0 : return "U";
			case 1 : return "M";
			case 2 : return "F";
		}
        return "U";
    }
	
	/**
		Returns the individual's affection status.
		            
		@return the individual's affection status
	*/		
    public int affection()
    {
		return affection;
    }
	
	/**
		Returns the individual's twin status in a String
		format.
		            
		@return the individual's twin status
	*/		
    public String twinString()
    {
		switch (twin) 
		{
            case 0 : return "0";
            case 1 : return "MZ";
            case 2 : return "DZ";
		}
        return "0";
	}
	
	/**
		Returns the individual's twin status as a numeric code
		            
		@return the individual's twin status
	*/		
    public int twinCode()
    {
		return twin;
	}
	
	/**
		Returns whether the individual is mock
		
		@return true if the individual is mock, false otherwise
	*/		
    public boolean isMock()
    {
		return isMock;
    }

	/**
		Whether a covariate is missing
	
		@return true if covariates are missing, false otherwise
	*/
	public boolean missesCovariates()
	{
		return covariates == null;
	}
  
	
	/**
		Returns the vector of phenotypes.
		
		@return  the vector of phenotypes
	*/    
    public double[] phenotypes()
    {
		if (phenotypes == null)
			return null;
		
		return phenotypes.clone();
    }
	
	/**
		Sets the i-th phenotype to the given value.
		
		@param i the position of phenotype
		@param value the value to set
	*/
    public void setPhenotypes(int i, double value)
    {
		phenotypes[i] = value;
    }
	
	/**
		Deletes the phenotype values 
	*/
    public void resetPhenotypes()
    {
		phenotypes = null;
    }
	
	
	/**
		Returns the vector of covariates.
		
		@return  the vector of covariates
	*/    
    public double[] covariates()
    {
		if (covariates == null)
			return null; 
		
		return covariates.clone();
    }
	
	/**
		Sets the covariate values 
		
		@param c the vector of covariates
	*/
    public void setCovariates(double[] c)
    {
		if (c == null)
			covariates = null;
		else
		 	covariates = c.clone();
    }
	
	/**
		Deletes the covariate values 
	*/
    public void resetCovariates()
    {
		covariates = null;
    }
	
}
