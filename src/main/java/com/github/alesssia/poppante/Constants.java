/*
 * 	 Constants.java
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

import java.lang.reflect.Field;
import java.lang.reflect.Modifier;
import org.apache.commons.cli.*;

/**
	Defines and manages the command line options.
	These values are used through all the project. 

	@author      Alessia Visconti <alessia.visconti@kcl.ac.uk>
	@author      Mario Falchi     <mario.falchi@kcl.ac.uk>
	@author      Hariklia Eleftherohorinou <hariklia.eleftherohorinou06@imperial.ac.uk>
	@version     1.0               
 */
class Constants 
{
	/** Default relatedness for mink-2 parameter
		
		It corresponds to the kinship between second cousins*/	
	private static final double SECOND_COUSIN_KINSHIP = 0.0315;
	
	/** Default relatedness for mink-3 parameter
		
		It corresponds to the kinship between third cousins*/	
	private static final double THIRD_COUSIN_KINSHIP = 0.0078;
	
	//Static variables and option names *must* have the same name
    
	/** Whether a help message should be printed on the standard output. */
	public static boolean help 	 = false; 
	/** Analysis to be performed.
      
		Allowed values are: "heritability", "association", "pedcheck". */
	public static int	  mode   = -1 ;   
	
	/**  File path of the pedigree file. 
		
		 The file should be in PLINK format along with
		 twin information. Each column after the 7th 
		 should describe a quantitative traits to be 
		 analysed.
		 If the family structure is unknown any 
		 value can be used to codify the familyID,
		 but the pair <familyID, individualID> should
		 univocally identify each subject.
		*/
	public static String  ped 	 = null; 
	/**  File path of the phenotype information file.
		
		 The file should contains the phenotype name. */
	public static String  response = null;  
	/**  File path of the methylation data file.
		
		 The file should contains the individuals' IDs as
		 described in the pedigree file. Each column after 
		 the 3rd should describe a methylation site
		 value. */
	public static String  predictor 	 = null;  
	/**  File path of the methylation site map file.
		
	     The file should contains the methylation
		 site coordinate. */
	public static String  map    = null;  
	
	/** File path of the output file. */                                     
	public static String  output = null;
	/** File path of the covariate file.
		
		The file should contains the individuals' IDs as
		described in the pedigree file. Each column after 
		the 3rd should describe a covariate	value.
		All the values are used in the analysis. */  
	public static String  covariate 	 = null;
	/** File path of the external kinship.
		
		The file should be formed by five columns: the first 
		four column describe a pairs of individuals, the 
		fifth their kinship value. */  
	public static String  kinship= null;  

	/** File path of the phenotypes to include in the analysis.
		
		The file should be formed by a column that 
		lists the names of the phenotypes to analyse */ 
	public static String filter = null;
	/** File path of the methylation sites to include in the analysis.
		
		The file should be formed by a column that 
		lists the names of the methylation sites to
		analyse */ 
	public static String include = null;
	/** Either the file path of the covariates used to correct
		the methylation values or the threshold of the overall variability
		to be removed using a Principal Component Analysis.
		
		The file should contains the individuals' IDs as
		described in the pedigree file. Each column after 
		the 3rd should describe a covariate	value.
		All the values are used in the analysis. */ 
	public static String correct = null;
	
	/** Whether verbose. */	
	public static boolean verbose  = false;
	/** Whether print Manhattan and Q-Q plot. */	
	public static boolean plot  = false;
	/** Whether to print the header in the output file. */ 	
	public static boolean header   = true; 
	/** Whether to print the variance in the output file. */ 
	public static boolean variance = false;
	/** Whether to perform the quantile normalisation on the quantitative traits. */ 
	public static String normalise = null; 
	/** Whether to force the usage the QR/LU decomposition
		
		It works only when the external kinship matrix is 
		supplied. The default is to use a bending procedure to
		make the matrix positive definite followed by the 
		Cholesky decomposition. */
	public static String decomposition = null;
 
	/** Window size (in bp) for the region-based tests. 
		
		When set it collapses the sites within the window 
		by means of their first principal component. */
	public static int     region = Utilities.INVALID_I; 
	
	/** Controls the experiment-wise error rate (EWER)
		
		In can be estimated using a Bonferroni adjustment
		with the estimate number of informative marker in the
		genome. Common thresholds for establishing significance 
		in EWAS are considered to be 5E-8 and 1E-8 (corresponding
		to a Bonferroni alpha of 0.05 and 0.01,	respectively). 
		Used in the adaptive permutation test. */
	public static double alpha = Utilities.INVALID_D;
	
	/** Controls the desired precision.
		
		Significance of the threshold that equals the 
		standard error in estimation, when the desired
		p-value is equal to alpha.  
		Used in the adaptive permutation test. */
	public static double c = Utilities.INVALID_D;
	
	/** Threshold for the extra statistics evaluation. 
		
		If the p-value of the test is smaller than relc 
		the percentage of families showing a positive 
		contribution and the Gini coefficient to assess 
		family contribution to the chi-square are evaluated
		and printed in the output file. */
	public static double  relc   = Utilities.INVALID_D; 
	/** Threshold for the kinship value. 
        
		Pairs of individuals with a kinship value smaller than mink 
		are considered unrelated and their kinship value set to zero. */
	public static double  mink   = 0.0; 		
	/** Number of threads to use.
         
		This value adjusted a run time to meet the 
		machine availability. */                                          
	public static int     threads = 1;  
	
	 /** List of avaliable options. */
	private static Option[] options;   
	

	/**
		Parses the command line. 
		
		It verifies that the analysis mode is set and correct and
		that the mandatory parameter are set as well.
		
		@postcondition the selected mode is valid
		
		@param args the command line
		@throws Exception if an error occurs when parsing the command line 
		@throws MandatoryParameterException if a mandatory parameter is not specified
		@throws IllegalModeException if the mode is not valid
		@see com.github.alesssia.poppante.MandatoryParameterException 
		@see com.github.alesssia.poppante.IllegalModeException 	
	*/
	public static void parse(String[] args) throws Exception, IllegalModeException, MandatoryParameterException
	{
		try 
		{
			Parser parser = new PosixParser();
			OptionSetter opts = new OptionSetter();
			CommandLine cli = parser.parse(opts.options, args);
			options = cli.getOptions();
			
			init(cli);
			
			assert mode != -1 : "Internal error: mode not assigned.";
		}
		catch(MissingArgumentException e) 
		{
			String[] tokens = e.getMessage().split(" ");
			String missingArgument = tokens[tokens.length-1];
			if (missingArgument.equals("help")) 			
			{
				verbose = true;
				help = true;
			}
			else
				throw new MandatoryParameterException("ERROR: Check the " + missingArgument + " parameter: its value might be missing or in the wrong format.\nUse the option --help for details about PopPAnTe usage.");
		}
		catch(Exception e)
		{
			throw e;
		}	
	}
	
	/**
		Reads the parameter values from the command line and
		sets the class static variable (options).
		
		If parameters are meaningless with respect to one another
		it also sets them accordingly.
		
		@param cli the command line
		@throws Exception if an error occurs when parsing the command line 
		@throws MandatoryParameterException if a mandatory parameter is not specified
		@throws IllegalModeException if the mode is not valid
		
		@see com.github.alesssia.poppante.MandatoryParameterException 
		@see com.github.alesssia.poppante.IllegalModeException 	
	*/
	private static void init(CommandLine cli) throws MandatoryParameterException, IllegalModeException, Exception
	{
		if (cli.hasOption("help"))
		{
			help = Boolean.parseBoolean(cli.getOptionValue("help"));
			if (help)
				return;
		}
			
		readMandatoryParameters(cli);
		readOptionalParameters(cli);
		
		// The evaluation of family contribution has no sense
		// when no families are supplied
		if (kinship != null)
			relc   = Utilities.INVALID_D;
		
		//Both alpha and c must be specified
		if (alpha != Utilities.INVALID_D && c == Utilities.INVALID_D || 
			alpha == Utilities.INVALID_D && c != Utilities.INVALID_D)
			throw new IllegalModeException("ERROR: both alpha and c must be specified in the adaptive permutation procedure.\nUse the option --help for details about PopPAnTe usage.");
                
        // The permutation tests is not performed for
		// heritability analysis.
		if (mode == Utilities.MODE_HERITABILITY)
			alpha = Utilities.INVALID_D;
		
		
		//If correct is a numeric value, then it must be 
		//included in 0-1
		if (correct != null && Utilities.isDouble(correct))
		{
			double v = Double.parseDouble(correct);
			if (v < 0 || v > 1)
				throw new IllegalModeException("ERROR: the value used for correction should be included in [0,1].\nUse the option --help for details about PopPAnTe usage.");	
		}
		
		
		//if the normalise mode is specified it should
		//be one of the following: methylation, phenotype, both
		if (normalise != null && !normalise.equals("response") && !normalise.equals("predictor") && !normalise.equals("both"))
			throw new IllegalModeException("ERROR: the normalisation option is not valid.\nUse the option --help for details about PopPAnTe usage.");	
		
		//if a decomposition is specified it should be one of the
		//following: QR, bending
		if (decomposition != null && !decomposition.equals("QR") && !decomposition.equals("LU"))
			throw new IllegalModeException("ERROR: the decomposition option is not valid.\nUse the option --help for details about PopPAnTe usage.");	
		
		//check if the mink values is valid (that is if it is <= 1)
		if (mink > 1)
			throw new IllegalModeException("ERROR: the minimum genomic relationship coefficient is not valid.\nUse the option --help for details about PopPAnTe usage.");
	}
	
	
	/**
		Reads mandatory parameters. 
		
		When the selected mode is "heritability" or "association"
		the mandatory parameters are: mode, pedigree file, the map of 
		the methylation sites, the methylation values and the
		information about the tested quantitative traits.
		If the selected mode  is "pedcheck", only the mode and the 
		pedigree file are considered mandatory and read.
		
		@param cli the command line
		@throws Exception if an error occurs when parsing the command line 
		@throws MandatoryParameterException if a mandatory parameter is not specified
		@throws IllegalModeException if the mode is not valid
		
		@see com.github.alesssia.poppante.MandatoryParameterException 
		@see com.github.alesssia.poppante.IllegalModeException 	
	*/
	private static void readMandatoryParameters(CommandLine cli) throws MandatoryParameterException, IllegalModeException, Exception
	{
		
		
		if (cli.hasOption("mode")) 
			readMode(cli.getOptionValue("mode"));
		else
			throw new MandatoryParameterException("ERROR: the analysis mode has not been specified as input.\nUse the option --help for details about PopPAnTe usage.");
				
		if (cli.hasOption("ped")) 
			ped = cli.getOptionValue("ped");
		else
			throw new MandatoryParameterException("ERROR: the pedigree file was not specified as input.\nUse the option --help for details about PopPAnTe usage.");

		//I don't need the other files
		if (Constants.mode == Utilities.MODE_PEDCHECK)
			return;	

		if (cli.hasOption("predictor")) 
			predictor = cli.getOptionValue("predictor");
		else
			throw new MandatoryParameterException("ERROR:the  predictor data file was not specified as input.\nUse the option --help for details about PopPAnTe usage.");


		if (cli.hasOption("map")) 
			map = cli.getOptionValue("map");
		else
			throw new MandatoryParameterException("ERROR: the map file for the predictors was not specified as input.\nUse the option --help for details about PopPAnTe usage.");
		
		//I don't need the phenotypes value
		if (Constants.mode == Utilities.MODE_HERITABILITY)
			return;	
		
		if (cli.hasOption("response")) 
			response = cli.getOptionValue("response");
		else
			throw new MandatoryParameterException("ERROR: the description file for the response variable was not specified as input.\nUse the option --help for details about PopPAnTe usage.");
		
		
		
	}
	
	/**
		Reads optional parameters.
		
		@param cli the command line
	*/
	private static void readOptionalParameters(CommandLine cli)
	{
		if (cli.hasOption("output"))
			output = cli.getOptionValue("output");

		if (cli.hasOption("covariate")) 
			covariate = cli.getOptionValue("covariate");
		
		if (cli.hasOption("correct")) 
			correct = cli.getOptionValue("correct");
		
		if (cli.hasOption("kinship")) 
			kinship = cli.getOptionValue("kinship");

		if (cli.hasOption("mink"))
		{
			String s = cli.getOptionValue("mink");
			
			if (s.toLowerCase().equals("c2"))
				mink = SECOND_COUSIN_KINSHIP;
			else if (s.toLowerCase().equals("c3"))
				mink = THIRD_COUSIN_KINSHIP;
			else 
				mink = readNumericParameter(cli, "mink");
		}

	
		if (cli.hasOption("verbose"))
		    verbose = readBooleanParameter(cli, "verbose");
		
		if (cli.hasOption("plot"))
		    plot = readBooleanParameter(cli, "plot");
			
		if (cli.hasOption("header"))
			header = readBooleanParameter(cli, "header");
		
		if (cli.hasOption("normalise"))
			normalise = cli.getOptionValue("normalise");
			
		if (cli.hasOption("variance"))
			variance = readBooleanParameter(cli, "variance");
		
		if (cli.hasOption("decomposition"))
			decomposition =  cli.getOptionValue("decomposition");
		
		if (cli.hasOption("region"))
			region = (int)readNumericParameter(cli, "region");
		
		if (cli.hasOption("threads"))
		    threads = (int)readNumericParameter(cli, "threads");
		
		if (cli.hasOption("relc")) 
			relc = readNumericParameter(cli, "relc");
		
	    if (cli.hasOption("filter"))
			filter = cli.getOptionValue("filter");
		
	    if (cli.hasOption("include"))
			include = cli.getOptionValue("include");	
     
		if (cli.hasOption("alpha"))
			alpha = readNumericParameter(cli, "alpha");
		
		if (cli.hasOption("c"))
			c = readNumericParameter(cli, "c");
	}
	
	
	/**
		Sets the analysis mode and verifies its validity.
			
		@param smode the string mode
		@throws IllegalModeException if the mode is not valid 
		@see com.github.alesssia.poppante.IllegalModeException
	*/
	private static void readMode(String smode) throws IllegalModeException
	{		
		switch (smode) 
		{
			case "pedcheck" : 
				mode = Utilities.MODE_PEDCHECK;
				break;
			case "heritability" :
				mode = Utilities.MODE_HERITABILITY;
				break;
			case "association" :
				mode = Utilities.MODE_ASSOCIATION;
				break;
			default:
				throw new IllegalModeException("ERROR: the selected mode is not available.\nUse the option --help for details about PopPAnTe usage.");
		}
					
	}
	
	/**
		Reads numeric parameters and verifies their validity.
		It is used also for integer parameters, that need to be casted 
		afterwards.
			
		@param option the option to read
		@throws IllegalModeException if the value is not numeric 
		@see com.github.alesssia.poppante.IllegalModeException
	*/
	private static double readNumericParameter(CommandLine cli, String option) throws IllegalModeException
	{
		double value;
		
		try
		{
			return Double.parseDouble(cli.getOptionValue(option));
		} 
		catch(NumberFormatException e)
		{
			throw new IllegalModeException("ERROR: the " + option + " parameter is not numeric.\nUse the option --help for details about PopPAnTe usage.");
		}

	}	
	
	/**
		Reads numeric parameters and verifies their validity.
		It is used also for integer parameters, that needs to be cast 
		afterwards.
			
		@param option the option to read
		@param toInt whether the value should be cast to integer
		@throws IllegalModeException if the value is not numeric 
		@see com.github.alesssia.poppante.IllegalModeException
	*/
	private static Boolean readBooleanParameter(CommandLine cli, String option) throws IllegalModeException
	{
		String value = cli.getOptionValue(option);
				
		if (!value.toLowerCase().equals("true") && !value.toLowerCase().equals("false"))
			throw new IllegalModeException("ERROR: the " + option + " parameter is not boolean.\nUse the option --help for details about PopPAnTe usage.");
		else 
			return Boolean.parseBoolean(value);
	}	
	
	/**
		Lists the available options.
		
		Every options is returned along with a short description.
			
		@return the list of options
	*/
	public static String optionList()
	{
		String s = "Available options:\n";
		
		s += "Mandatory parameters:\n\n";
		
		s += "\t-mode <mode>\t\twhich analysis perform:\n";
		s += "\t\t\t\tmode=<heritability | association | pedcheck>\n";
		
		
		s += "\t-map file path\t\tmap for the predictors (not mandatory mode\n\t\t\t\tis 'pedcheck')\n";
		s += "\t-ped file path\t\tpedigree file\n";
		s += "\t-predictors file path\tpredictors values (not mandatory if mode is\n\t\t\t\t'pedcheck')\n";
		s += "\t-response file path\tresponse variables information (mandatory only if mode\n\t\t\t\tis 'association')\n";
		
		
		s += "\nOptional parameters:\n\n";
		
		s += "\t[-alpha num]\t\tp-value that controls the experiment-wise error \n\t\t\t\trate in the adaptive permutation procedure\n\t\t\t\t(default: null)\n";
		s += "\t[-c num]\t\tdesired precision in the adaptive permutation\n\t\t\t\tprocedure (default: null)\n";
		s += "\t[-correct <th|path>]\tthe file of covariates used to correct the\n\t\t\t\tpredictors OR the threshold of the total\n\t\t\t\tvariability to be removed by PCA\n";
		s += "\t[-covariate file path]\tcovariate file\n";
		s += "\t[-decomposition <QR|LU>]apply the QR/LU decomposition when the genetic\n\t\t\t\trelationship matrix is provided as input (default:\n\t\t\t\tapply a bending procedure and use the Cholesky\n\t\t\t\tdecomposition)\n";
		s += "\t[-filter file path]\tlist of responses to tests (default: all)\n";
		s += "\t[-include file path]\tlist of predictors to tests (default: all)\n";
		s += "\t[-header <true|false>]\twhether the output file has a header (default: true)\n";
		s += "\t-help\t\t\tprint this message\n";
		s += "\t[-kinship file path]\tgenetic relationship matrix file (default: null)\n";
		s += "\t[-mink threshold]\tminimum genomic relationship coefficient\n\t\t\t\t(default: 0)\n";
		s += "\t[-mink 2]\t\tset the minimum genomic relationship coefficient to " + SECOND_COUSIN_KINSHIP + "\n";
		s += "\t[-mink 3]\t\tset the minimum genomic relationship coefficient to " + THIRD_COUSIN_KINSHIP + "\n";
		s += "\t[-normalise <what>]\twhether quantile normalisation is applied to\n\t\t\t\tresponses, predictors, or both\n\t\t\t\t(default: none) what=<response|predictor|both>\n";
		s += "\t[-output file path]\toutput file (default: standard output)\n";
		s += "\t[-plot <true|false>]\twhether print the Manhattan and Q-Q plot (default: false)\n";
		s += "\t[-region bp]\t\twindow size for the region-based tests (default: 0,\n\t\t\t\tsingle-predictor analysis)\n";
		s += "\t[-relc threshold]\twhether evaluating additional statistics and which\n\t\t\t\tp-value threshold use (default: none)\n";		
		s += "\t[-threads num]\t\tnumber of threads to use (default: 1)\n";
		s += "\t[-variance <true|false>]whether printing the variances (default: false)\n";
		s += "\t[-verbose <true|false>]\twhether verbose (default: false)\n";
				
		return s;
	}
	
	
	/**
		Lists the option in effects.
		
		The default parameters are not reported, with the exception of -mink.
		
		@return the list of option in effect
	*/
	public static String printConstants()
	{
		String s = "Options in effect:\n";
        for (Option option : options) 
			if (!option.getLongOpt().equals("mink"))
            	s += "\t-" + option.getLongOpt() + " : " + option.getValue() + "\n";

		//mink is not printed above, but always added at the end (so the user
		//is aware that the defualt is zero (and negative kinship values are
		//constrained to zero)
		s += "\t-mink : " + mink + "\n";

		return s;
	}
	
}


/**
	Creates and manages the command line options.	

	@author      Alessia Visconti <alessia.visconti@kcl.ac.uk>
	@author      Mario Falchi     <mario.falchi@kcl.ac.uk>
	@author      Hariklia Eleftherohorinou <hariklia.eleftherohorinou06@imperial.ac.uk>
	@version     1.0               
 */
class OptionSetter
{
	/** List of options */
	Options options; 
	
	/**
		Constructor. 
		
		It sets the list of options in the Constants 
		class field.
		
		@see com.github.alesssia.poppante.Constants 	
	*/
	public OptionSetter()
	{
		options = new Options();
		init();
	}
	
	/**
		Initialises the list of options.
		
		It sets the list of options in the Constants 
		class field.
		
		@see com.github.alesssia.poppante.Constants 	
	*/
	private void init()
	{
		Field[] f = Constants.class.getFields();
		for (Field f1 : f) 
			if (Modifier.isStatic(f1.getModifiers())) 
				options.addOption(OptionBuilder.withLongOpt(f1.getName()).withValueSeparator(',').hasArgs().create());
	}
}


/**
	Runtime exception. 

    Raised when a mandatory parameter is not specified.

	@author      Alessia Visconti <alessia.visconti@kcl.ac.uk>
	@version     1.0               
 */
class MandatoryParameterException extends RuntimeException
{
	/**  serial number */
	public static final long serialVersionUID = 4242L; 

	/**
		Constructor. 
		
		Initialises the exception message.
		
		@param msg the exception message
	*/
	public MandatoryParameterException(String msg) 
	{
		super(msg);
	}
}

/**
	Runtime exception. 

	Raised when an illegal analysis mode is specified.	

	@author      Alessia Visconti <alessia.visconti@kcl.ac.uk>
	@version     1.0               
 */
class IllegalModeException extends RuntimeException
{
    private static final long serialVersionUID = 2L;
	/**
		Constructor. 
		
		Initialises the exception message.
		
		@param msg the exception message
	*/
	public IllegalModeException(String msg) 
	{
		super(msg);
	}
	
}

