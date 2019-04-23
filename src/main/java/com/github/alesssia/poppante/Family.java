/*
 * 	 Family.java
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

/** 
	Represents a family.

	@author      Alessia Visconti <alessia.visconti@kcl.ac.uk>
	@author      Mario Falchi     <mario.falchi@kcl.ac.uk>
	@version     1.0    
    @see com.github.alesssia.poppante.Person
	@see com.github.alesssia.poppante.Kinship
 */
class Family 
{
	/** Family ID */
	private final String id; 	
	/** Family members */		   
	private Vector<Person> members;   
	/** Members' ID order */ 
	private final Vector<String> memberIDs;  
	/** Kinship matrix */ 
	Kinship kinship;        
	/** Whether a monozygotic twin belongs to the family */
	private boolean hasMZTwin; 		 
	
	
	/**
		Constructor. 
		
		Initialises the data structure and sets 
		the family's ID.
		
		@param id the family's ID
	*/
	public Family(String id) 
	{
		this.id = id;
		members = new Vector<>();
		memberIDs = new Vector<>();
			
		kinship = new Kinship();
		hasMZTwin = false;
	}
	
	/**
		Returns the family ID.
		
		@return the family ID
	*/
	public String id()
	{
		return id;
	}
	
	
	/**
		Returns the list of people belonging to the family.
		
		@return the family's member
	*/	
	public Vector<Person> members()
	{
		return members;
	}
	
	
	/**
		Returns the list of members' ID.
		
		@return the list of members' ID
	*/
	public Vector<String> memberIDs()
	{
		return memberIDs;
	}
	
	/**
		Returns a family member given her position.
		
		@precondition the position is valid
		@param pos the position in the family members' list
		@return the family member
	*/
	public Person getMember(int pos) 
	{
		assert pos >= 0 && pos < members.size() : "Internal error: asking for member in a wrong position";
		return members.get(pos);
	}
	
	
	/**
		Returns a family member given her ID. 
		
		If the ID is not present, null is returned.
		
		@param id the ID of the family member
		@return the family member or null
	*/
	public Person getMemberByID(String id) 
	{
		Person p = null;
		for (Person member : members) 
		{
			if (member.id().equals(id)) 
			{
				p = member;
				break;
			}
		}
		return p;
	}

	/**
		Returns the number of members.
		
		@return the number of members
	*/
	public int numMembers() 
	{
		return members.size();
	}
	
	/**
		Returns whether the family includes a monozygotic twin.
		
		@return true if the family includes a monozygotic twin, false otherwise
	*/
	public boolean hasMZTwin() 
	{
		return hasMZTwin;
	}

	/**
		Returns the kinship matrix.
		
		@return the kinship matrix
	*/
	public Kinship kinship()
	{
		return kinship;
	}

	/**
		Adds a new member to the family.
		
		If the new member is a monozygotic twin (and if she belongs 
		to the offspring) records her presence.
		It also stores the position of this family member in the family. 
		To do so it uses the memberID list. Storing the member order 
		is a key aspect since ancestors must precede offsprings in 
		the kinship evaluations.
		
		It checks whether the family ID is correct and whether the
		individual is not already present (that is there is no one with
		her ID).
		
		@precondition  the individual's ID should match the family ID
		@throws java.lang.Exception if the person is already present in the dataset
		@param p the member to add
     * 
	*/
	public void addMember(Person p) throws Exception
	{
		assert p.famID().equals(this.id) : "Internal error: assigning individual to the wrong family.";
		
		if (getMemberByID(p.id()) != null)
			throw new Exception("An individual with ID " + p.id() + " is already present in family with ID " + this.id);
			
		if (p.isMZ() && !p.isFounder())
			hasMZTwin = true;
			
		memberIDs.add(p.id());
		members.add(p);	
	}
	
		
	/**
		Rearranges the members.
		
		It arranges the members in the members list following
		the order given by the memberIDs list.	
	*/
	private void rearrange()
	{
		Vector<Person> tmp = new Vector<>();
		for (String memberID : memberIDs) 
			tmp.add(getMemberByID(memberID));
	    
		members = tmp;
	}
	

	/**
		Swaps two member ID in the memberID list given
		their position.
		
		@param a the first position
		@param b the second position
	*/
	private void swapMembers(int a, int b)
	{
		String aID = memberIDs.get(a);
		String bID = memberIDs.get(b);
		memberIDs.set(b, aID);
		memberIDs.set(a, bID);
	}


	/**
		Sorts the family members.
		
		The sorting is performed because when the Kenneth 
		Lange's approach is used the parents must preceed 
		their offsprings. The Kenneth Lange's approach is 
		used to evaluate the theoretical kinship matrix.
		The sorting is not performed when the kinship 
		information is supplied by means of an external
		file since the Kenneth Lange approach will not be
		used.
		
		@see com.github.alesssia.poppante.DataManager
	*/	
	public void sort() 
	{
		int current = 0;	
		
		while(current < memberIDs.size())
		{			
			Person person = getMemberByID(memberIDs.get(current));
			
			//ancestors are ok
			if (person.isFounder())
			{
				current++;
				continue;
			}
			
			//checks the position of the father, and swap the two member if necessary
			int fIndex = memberIDs.indexOf(person.fatherID());
			if (current < fIndex)
			{	
				swapMembers(current, fIndex);
				continue;
			}
			
			//if the father is ok, checks the position of the father, 
			//and swap the two member if necessary
			int mIndex = memberIDs.indexOf(person.motherID());
			if (current < mIndex)
			{
				swapMembers(current, mIndex);
				continue;
			}

			//done for this position
			current++;
		}
		
		//mambers will have the same order of memberID
		rearrange();
	}
	
	
	/**
		Evaluates the kinship matrix (Kenneth Lange's method).
		
		It requires the family members to be sorted, in order to have 
		the ancestors before the offsprings (as done by the method]
		sort of this class) and the kinship initialised.
		
		If the family includes monozygotic twins the 
		kinship matrix is corrected to consider this information.
		
		@precondition  family member must be sort
		@precondition  kinship matrix must have been initialised
		
		@see com.github.alesssia.poppante.Family.sort
		@see com.github.alesssia.poppante.Family.rearrange
		@see com.github.alesssia.poppante.Kinship
	*/	
	public void evaluateKinship()
	{
		assert kinship.matrix() != null : "Internal error: kinship matrix for family " + id + " has not been initialised.";
		
		kinship.evaluate(memberIDs, members, hasMZTwin);	
	}
	
	
	/**
		Checks the member's parents.
		
		Return a message listing whether the member's parents
		contain errors such as:
		- missing information;
		- wrong gender.
		It is used when the analysis mode is "pedcheck".
	
		@param p the family member to check
		@return a message listing the errors
	*/
	private String checkParents(Person p)
	{
		if (p.isFounder())
			return "";
		
		String message = "";
		Person father = getMemberByID(p.fatherID());
		Person mother = getMemberByID(p.motherID());
		
		if (father == null)
			message += "\t- In family " + id + " subject " + p.id() + " has no father information.\n";
		else if (!father.isValidParent(Utilities.MALE))
			message += "\t- In family " + id + " subject " + father.id() + " has been assigned to be subject " + p.id() + "'s father, but he is not male.\n";
		
		if (mother == null)
			message += "\t- In family " + id + " subject " + p.id() + " has no mother information.\n";
		else if (!mother.isValidParent(Utilities.FEMALE))
			message += "\t- In family " + id + " subject " + mother.id() + " has been assigned to be subject " + p.id() + "'s mother, but she is not female.\n";
		
		return message;
	}
	
	
	/**
		Checks the member's siblings.
		
		Return a message listing whether the member's siblings
		contain errors such as:
		- missing twins;
		- monozygotic twins have different gender.
		
		Two members are considered to be a pair of twins if their 
		twin status is the same (either monozygotic or dizygotic) 
		and if they share the same parents. 
		Founders are not considered.
		It is used when the analysis mode is "pedcheck".
		
		@param p the family member to check
		@return a message listing the errors
	*/
	private String checkSiblings(Person p)
	{
		//problems arise when the subject is a twin
		String twinStatus = p.twinString();
		if (p.isFounder() || twinStatus.equals("0"))
			return "";
		
		boolean found = false;
		String message = "";
		for (Enumeration<Person> elements = members.elements(); elements.hasMoreElements();)
		{
			Person tmp = elements.nextElement();
			//skip myself
			if (tmp.id().equals(p.id()))
				continue;
			
			//we don't care about founder information, and we consider two people a twin pair if
			//they have the same twin status and the same pair of parents
			if (!tmp.isFounder() && tmp.twinString().equals(twinStatus) && p.fatherID().equals(tmp.fatherID()) && p.motherID().equals(tmp.motherID()))
			{
				//I've found the other twin, let's check geneder if necessaty (MZ twin)
				found = true;
				if (twinStatus.equals("MZ") && p.sexCode() != tmp.sexCode())
					message += "\t- In family " + id + " subjects " + p.id() + " and " + tmp.id() + " should represent a monozygotic twin pair, but they have different sex.\n";
				
			}
		}
		
		if (!found)
			message +=  "\t- In family " + id + " subjects " + p.id() + " has not information about her " + twinStatus + " twin.\n";
			
		return message;
	}
	
	
	/**
		Checks the family's structure.
		
		Returns a message listing whether the family structure 
		contains errors such as:
		- missing parents;
		- wrong gender assigned to parents;
		- missing twins;
		- monozygotic twins' gender.
		It is used when the analysis mode is "pedcheck".
		
		@return a message listing the errors
	*/
	public String checkFamily()
	{
		String message = "";
		for (Enumeration<Person> elements = members.elements(); elements.hasMoreElements();)
		{
			Person p = elements.nextElement();
			message += checkParents(p);
			message += checkSiblings(p);
		}	
		
		return message;
	}
}