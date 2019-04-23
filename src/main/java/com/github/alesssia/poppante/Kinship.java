/*
 * 	 Kinship.java
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
import Jama.*;

/**
	Represents a kinship matrix as an upper 
	triangular matrix.

	@author      Alessia Visconti <alessia.visconti@kcl.ac.uk>
	@version     1.0     
*/

class Kinship
{
	/** Kinship matrix */
	private double[][] matrix; 	   		
	
	/**
		Constructor. 
		
		Builds an empty kinship matrix
	*/	
	public Kinship()
	{
		matrix = null;
	}
	
	/**
		Constructor. 
		
		Builds a matrix of the given size
		
		@param n the number of family members 
	*/	
	public Kinship(int n)
	{
		init(n);
	}
	
	/** 
		Returns the kinship matrix
		
		@return the kinship matrix
	*/
	public double[][] matrix()
	{
		if (matrix == null)
			return null;
		
		return matrix.clone();
	}
	
	/**
		Initialises the kinship (triangular) matrix.
		
		All values are equals to zero.
		
		@param n the number of family members 
	*/
	public void init(int n)
	{
		matrix = new double[n][];
		for (int i = 0; i < n; i++)
			matrix[i] = new double[n - i];
	}
	
	
	/**
		Evaluates the kinship matrix (Kenneth Lange's method).
		
		It requires the family member to be sorted, in order to have 
		the ancestors before the offsprings. 
		If the family includes monozygotic twins the kinship matrix is 
		corrected to consider this information.
		
		@param memberIDs the sorted list of family member IDs
		@param members the list of family members
		@param hasMZTwin whether the family has a monozygotic twin pair
		
		@see com.github.alesssia.poppante.Family
	*/
	public void evaluate(Vector<String> memberIDs, Vector<Person> members, boolean hasMZTwin) 
	{
		double value;
		for (int current = 0; current < memberIDs.size(); current++)
		{
			Person person = members.get(current);
		
			//sets founder to be 0.5 fixing part of the matrix
			//diagonal, since each subject is identical to herself
			if(person.isFounder())
			{
				setValue(0.5, current, current);
				continue;
			}
			
			//sets the rest of the diagonal taking also in account the 
			//parents inbreed relationship. If parents are unrelated 
			//then the value is set to be 0.5 (subjects are considered 
			//identival to themselves. Otherwise the inbreed relationship
			//is weigthed and summed. 
			int fIndex = memberIDs.indexOf(person.fatherID());
			int mIndex = memberIDs.indexOf(person.motherID());
			value = 0.5 + (0.5 *  getValue(fIndex, mIndex));
			setValue(value, current, current);
			
			//adjusts the genetic relationships between the current subject 
			//and the other family members condidering the genetic relatioship
			//between them and the current subject's parent
			//When the family member is the current subject MZ twin, the value 
			//is set to be 0.5, considering the two twins the same individual 
			for (int i=0; i < current; i++)
			{
				if (hasMZTwin && isSheMe(members, current, i))
					setValue(0.5, current, i);
				else
				{
					double fValue = getValue(fIndex, i);
					double mValue = getValue(mIndex, i);
					value = 0.5 * (fValue + mValue);
					setValue(value, current, i);
				}
			}
		}	
		
		//we need to multiply all the values by two, as definite
		//in the Kenneth Lange's method
		for (int i = 0; i < matrix.length; i++)
			for (int j = 0; j < matrix[i].length; j++)
				matrix[i][j] *= 2;
	}
	
	/**
		Sets the self kinship value of each individual to 1.
		
		It is used when the kinship matrix is provided as 
		external file. It assures that the self-relatedness is
		set also when it is not reported in the file
		(as done, for instance, by PLINK).
		If the file reports the self-relatedness this information
		will be overwritten.
		
		@see com.github.alesssia.poppante.DataManager
	*/
	public void adjustKinship()
	{
		// (The external values are already multiplied by two)
		for (int i = 0; i < matrix.length; i++)
			setValue(1.0, i, i);
	}
	
	
	/**
		Returns whether two subjects are monozygotic twins. 
		
		Two subjects are considered monozygotic twins if:
		- they are labelled as such;
		- they have the same parents; and
		- they have the same gender.
		Please note that pairs of monozygotic twins are considered identical twins of
		one another.
		
		@param members the list of family members
		@param myIndex the index of the first individual
		@param herIndex the index of the second individual
		@return true if the two subjects are monozygotic twins, false otherwise
	*/
	private boolean isSheMe(Vector<Person> members, int myIndex, int herIndex)
	{		
		Person me = members.get(myIndex);
		
		//If I'm not MZ, she can't be me
		if(!me.isMZ())
			return false;
		
		Person she = members.get(herIndex);
		return (she.isMZ() && me.fatherID().equals(she.fatherID()) && me.motherID().equals(she.motherID()) && me.sexCode() == she.sexCode());
	}
	
	
	/**
		Sets an entry of the kinship matrix to a value.
		
		@param value the value
		@param row the row index 
		@param col the col index
	*/
	public void setValue(double value, int row, int col) 
	{
		if (row <= col)
			matrix[row][col-row] = value;
		else
			matrix[col][row-col] = value;
	}
	
	/**
		Returns the value of an kinship matrix entry.
		
		@param row the row index 
		@param col the col index
		@return the value
	*/
	public double getValue(int row, int col) 
	{
		if (row <= col)
			return matrix[row][col-row];
		else
			return matrix[col][row-col];
	}
	
	
	/**
		Removes individuals from the kinship matrix
		
		@param position the position to remove
		@param n the number of individuals in the new matrix
	*/
	public void reset(Vector<Integer> position, int n)
	{
		double[][] oldK = matrix;
		init(n);
		
		int k=0;
		for (int i=0; i<oldK.length; i++)
		{
			//is this line to be removed?
			if (position.contains(i))
				continue;
			
			int l=0;
			for (int j=0; j<oldK.length; j++)
			{
				//is this coloum to ve removed?
				if (position.contains(j))
					continue;
				
				double value = (i <= j) ? oldK[i][j-i] : oldK[j][i-j];
				setValue(value, k, l);
				l++;	
			}
			k++;
		}
	}
	
	
	
	/**
		Applies the 'bending' procedure to modify the eigenvalues of non
		positive defined (external) kinship matrix.
		
		Please note that no check is performed to verify whether the
		matrix is or not positive definite. 
				
		See http://www.aps.uoguelph.ca/~lrs/Summer2012Full/PDforce.pdf for
		the R function that inspired this code
	*/
	public void bending()
	{	
		int n = matrix.length;
		
		//creates a symmetric matrix
		double[][] tmp = new double[n][n];
		for (int i=0; i<n; i++)
			for (int j=i; j<n; j++)
				tmp[i][j] = tmp[j][i] = getValue(i, j);
	
		//	Eigenvalues and eigenvectors		
		Matrix a = new Matrix(tmp);
		EigenvalueDecomposition eig = new EigenvalueDecomposition(a);
		double[] eigenvalues = eig.getRealEigenvalues();
		
		// tolerance level
		double rtol = Utilities.TOL * eigenvalues[0];
		boolean positivedefinite = true;
		for (int i=0; i<n; i++)
			if (eigenvalues[i] < rtol)
			{
				positivedefinite = false;
				eigenvalues[i] = rtol;
			}
				
		//no eigenvalues are smaller than the tolerance level,
		//and the matrix should be positive definite
		if (positivedefinite)
			return; 
	
		//does the trick
		Matrix eigenvectors = eig.getV();
		tmp = new double[n][n];
		for (int i=0; i<n; i++)
			tmp[i] = eigenvalues;
        a = new Matrix(tmp);
		a = a.transpose();
		a = a.arrayTimes(eigenvectors.transpose());
		tmp = eigenvectors.times(a).getArray();
		
		//update the kinship matrix
		for (int i=0; i<n; i++)
			for (int j=i; j<n; j++)
				setValue(tmp[i][j], i, j);
	}
	
}