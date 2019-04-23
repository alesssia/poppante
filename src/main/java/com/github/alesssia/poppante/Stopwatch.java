/*
 * 	 Stopwatch.java
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
	Implements a simple stopwatch to measure elapsed time

	@author      Alessia Visconti <alessia.visconti@kcl.ac.uk>
	@version     1.0     
*/

class Stopwatch
{
	/**  Start time (nanoseconds) */
	private long startTime;		
	/**  Stop time (nanoseconds) */
	private long stopTime;		
	
	/**
		Constructor. 
		
		Initialises the timers.
	*/
	public Stopwatch()
	{
		startTime = 0;
		stopTime = 0;
	}
	
	/**
		Starts the stopwatch
	*/	
	public void start()
	{
		startTime = System.nanoTime();
	}
	
	/**
		Stops the stopwatch
	*/
	public void stop()
	{
		stopTime = System.nanoTime();
	}
	
	/**
		Resets the stopwatch
	*/
	public void reset()
	{
		startTime = 0;
		stopTime = 0;
	}
	
	/**
		Returns the elapsed time in milliseconds
		
		@return elapsed time
	*/
	private long getMillisecond()
	{
		return (stopTime - startTime)/1000000;
	}
	
	
	/**
		Returns the elapsed time in the format dd:hh:mm:ss
		
		@return elapsed time
	*/
	public String getTime()
	{
		long milliseconds = getMillisecond();
		
		long seconds = milliseconds / 1000;
		milliseconds -= seconds*1000;
		
		long minutes = seconds/60;
		seconds %= 60;
		
		long hours = minutes / 60;
		minutes %= 60;
		
		if (minutes == 0)
			return  minutes + "m:" + seconds + "s:" + milliseconds + "ms"; 

		return hours + "h:" + minutes + "m:" + seconds + "s"; 
	}
	
}
