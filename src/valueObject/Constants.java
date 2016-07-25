package valueObject;

public class Constants {
	//number of parameters should be found in the parameter file
	public static final int NUM_PARAMETERS = 19;
	
	public static final String NUM_KEY = "NUM";
	public static final String NBR_OF_CHR_KEY = "NBR_OF_CHR";
	public static final String ADJACENT_DIST_KEY = "ADJACENT_DIST";
	public static final String CONTACT_DIST_KEY = "CONTACT_DIST";
	public static final String POS_MIN_DIST_KEY = "POS_MIN_DIST";
	public static final String NEG_MAX_DIST_INTRA_KEY = "NEG_MAX_DIST_INTRA";
	public static final String NEG_MAX_DIST_INTER_KEY = "NEG_MAX_DIST_INTER";
	
	public static final String POS_MAX_DIST_WEIGHT_FILE_KEY = "POS_MAX_DIST_WEIGHT_FILE";
	public static final String POS_MIN_DIST_WEIGHT_FILE_KEY = "POS_MIN_DIST_WEIGHT_FILE";
	public static final String NEG_MIN_DIST_WEIGHT_FILE_KEY = "NEG_MIN_DIST_WEIGHT_FILE";
	public static final String NEG_MAX_DIST_WEIGHT_FILE_KEY = "NEG_MAX_DIST_WEIGHT_FILE";
	
	public static final String INTRA_IF_THRESHOLD_KEY = "INTRA_IF_THRESHOLD";
	public static final String INTER_IF_THRESHOLD_KEY = "INTER_IF_THRESHOLD";
	public static final String OUTPUT_FOLDER_KEY = "OUTPUT_FOLDER";
	//public static final String FILE_PREFIX_KEY = "FILE_PREFIX"; // remove on Mar 9, use Input file name
	//public static final String FILE_HEADER_KEY = "FILE_HEADER"; // remove on Mar 9, use default Tuan
	public static final String INPUT_FILE_KEY = "INPUT_FILE";
	public static final String VERBOSE_KEY = "VERBOSE";	
	public static final String CHR_UPPER_BOUND_ID_FILE_KEY = "CHR_UPPER_BOUND_ID_FILE";
	
	public static final String LEARNING_RATE_KEY = "LEARNING_RATE";
	
	public static final String MAX_ITERATION_KEY = "MAX_ITERATION";
	
	//maximum number of threads should be used 
	public static final int MAX_NUM_THREAD = 240;
	
	//the starting learning rate for the line search
	public static double INITIAL_LEARNING_RATE = 500;
		
	
	//maximum number of iterations
	public static int MAX_ITER = 100000;
	
	//this constant is used to check if the norm of the gradient is near zero
	public static final double NEAR_ZERO = 10e-6;
	
	//distane larger than this value is not sensitive for the tanh function and needed to be scale down
	public static final double LARGE_DISTANCE_FOR_TANH = 20;
	//if the distance is larger than LARGE_DISTANCE_FOR_TANH, it will be scale down to this value
	public static final double SCALE_DISTANCE = 10;

}
