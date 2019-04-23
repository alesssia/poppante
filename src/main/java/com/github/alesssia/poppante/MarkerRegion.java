/*
 * 	 MarkerRegion.java
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
import com.github.alesssia.algebrautils.MyPCA;

/**
	Collapses multiple markers by means of their first 
	Principal Component.

	It calculates the PCA only over a specific list of
	values out of those given as input. 
	Also it uses only the values that will be used in the
	analysis that is only those that belong to a family 
	member having a valid value for the methylation site m 
    (when the analysis mode is "heritability") or both 
	the phenotype p and methylation site m (when the analysis
    mode is "association").

	@author      Alessia Visconti <alessia.visconti@kcl.ac.uk>
	@version     1.0     
	@see com.github.alesssia.poppante.VC	
	@see com.github.alesssia.poppante.Marker	
*/ 
class MarkerRegion
{
	/** Values to be summarised */
	private double[][] data;				
	
	/**
		Constructor. 
		
		Initialises the data structure.
		
		@param sites methylation sites to be summarised
		@param methylations the matrix of methylation values
		
		@see com.github.alesssia.poppante.VC	
		@see com.github.alesssia.poppante.Marker	
	*/
	public MarkerRegion(Vector<Integer> sites, double[][] methylations)
	{
		data = new double[sites.size()][methylations[0].length];
		for (int i=0; i<sites.size(); i++)
			data[i] = Arrays.copyOf(methylations[sites.get(i)], methylations[sites.get(i)].length);		
	}
	
	/**
		Evaluates the PCA of the selected methylation sites and returns 
		the first principal component.
		
		If there is only one maker in the list, it returns the
		original methylation values.
		
		It substitutes the missing values with zero. However, a matrix
		line (that is one of the PCA variables) composed only by zero 
		changes the PCA results. It is then necessary to remove all the
		mock individuals (that are those having only missing values) before
		calling this function.

		@precondition no mock individuals should be available in the main matrix
		@return the summary of the given methylation sites
	*/
	public double[] summarise()
	{
		if (data.length == 1)
			return data[0];
		
		//Variables should stay in the columns (I already creted the matrix
		//that way)
		MyPCA pca = new MyPCA(data, Utilities.INVALID_D);
		pca.evaluatePCA();
		
		return pca.getPC(1);
	}
}