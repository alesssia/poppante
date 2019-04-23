/*
 * 	 Thread.java
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

/**
	Manages the number of threads.

	@author      Alessia Visconti <alessia.visconti@kcl.ac.uk>
	@version     1.0     
*/


class Thread
{
	/** The number of threads */
	private final int threads;
	
	/**
		Constructor. 
		
		Decides the number of threads to use. The number is the smallest
		values between the requested number of thread and the number of 
		available processors. 
		
		If no threads are requested, a single thread is assigned
		
		@see com.github.alesssia.poppante.Constants
	*/	
	public Thread()
	{
		threads = Math.min(Runtime.getRuntime().availableProcessors(), Constants.threads);
	}
	
	/** 
		Returns the number of threads to use
		
		@return the number of threads to use
	 */		
	public int threads()
	{
		return threads;
	}
	
		
}
	
	