package genomeReconstruction;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Date;
import java.util.HashMap;
import java.util.regex.Pattern;

import optimization.GradientAscent;
import optimization.OptimizedObject;
import utility.Helper;
import valueObject.Constants;
import valueObject.ContactMT;
import valueObject.StructureEvaluationVO;
import evaluation.StructureScoring;
//test to see if syn is still working
public class GenomeGenerator implements OptimizedObject{

	//number of structures will be generated
	private int NUM;
	
	//number of chromosomes
	private int NBR_OF_CHR;
	
	//upper bound absolute id (id counted as 1..n and including the centermere) of each chromosome
	//the first element, with index of 0 is ignored
	private int[] CHR_UPPER_BOUND_ID;
	
	//threshold for intra-chromosomal contacts
	private double INTRA_IF_THRESHOLD;
	private boolean is_intra_thres_percent = false; // to signify if thres is in percentage
	private double INTRA_PERCENT; // used to determine how initialization should be
	
	//threshold of inter-chromosomal contacts
	private double INTER_IF_THRESHOLD = Double.NaN;
	private boolean is_inter_thres_percent = false;
	
	//maximum distance (in square) between 2 adjacent points
	private double ADJACENT_DIST;
	//maximum distance (in square) between 2 points that are in contact (positive contacts), usually, this is the contact distance
	private double POS_MAX_DIST;
	//minimum distance (in square) between 2 points that are not in contact (negative contacts), usually, this is the contact distance
	private double NEG_MIN_DIST;
	//minimum distance (in square) between 2 points that are in contact (positive contacts)
	private double POS_MIN_DIST;
	//maximum distance (in square) between 2 points that are not in contact (negative contacts) and in the same chromosome
	private double NEG_MAX_DIST_INTRA;	
	//maximum distance (in square) between 2 points that are not in contact (negative contacts) and not in the same chromosome
	private double NEG_MAX_DIST_INTER;
	
	//learning for the gradient ascent
	private double LEARNING_RATE;
	
	//index of 0 is ignored
	//importance of the constraints for max distance of contacts (positive contacts), the larger it is, the bigger the corresponding score will be
	private double[][] POS_MAX_DIST_WEIGHT;
	//index of 0 is ignored
	//importance of the constraints for min distance of contacts (positive contacts), the larger it is, the bigger the corresponding score will be
	private double[][] POS_MIN_DIST_WEIGHT;
	//index of 0 is ignored
	//importance of the constraints for min distance of non-contacts (negative contacts), the larger it is, the bigger the corresponding score will be
	private double[][] NEG_MIN_DIST_WEIGHT;
	//index of 0 is ignored
	//importance of the constraints for max distance of non-contacts (negative contacts), the larger it is, the bigger the corresponding score will be
	private double[][] NEG_MAX_DIST_WEIGHT;

	private String POS_MAX_DIST_WEIGHT_FILE;
	private String POS_MIN_DIST_WEIGHT_FILE;
	private String NEG_MIN_DIST_WEIGHT_FILE;
	private String NEG_MAX_DIST_WEIGHT_FILE;
	private String CHR_UPPER_BOUND_ID_FILE;

	//to indicate if information during optimization should be printed out or not
	private boolean VERBOSE = true;
	
	//the output folder, structures will be put here
	private String OUTPUT_FOLDER;
	
	//prefix for file name of structures
	private String FILE_PREFIX;
	
	//header of file structures
	private String FILE_HEADER="";
	
	//the file contains matrix contact
	private String INPUT_FILE;

	//file contains paramters needed to run the program
	private String parametersFile;

	//zoom factor, in case the contact distance is large
	private double zoomFactor = 1.0;
	
	//contains interaction frequencies of fragments
	private ContactMT contactMT;
	//private double[][] contactMT;
	
	//number of units of the whole genome
	private int n;
	
	//sum of all interaction frequency
	private double totalIF;	
	private double normTotalIF;//to use in place of totalIF to avoid numerical issues
	
	//genome structure, every three indices is one point in 3D, this will be variables
	private double[] str;
	
	//number of threads of the program
	private int num_of_processors;
	
	private int max_iteration = Constants.MAX_ITER;
	
	//to print out the log of the program
	private PrintWriter logPW = null;
	
	//to map each index to the corresponding chromosome number
	private HashMap<Integer,Integer> idToChr = new HashMap<Integer,Integer>();
	
	//to store actual indices in the genome, it is used to identify if 2 indices in the structure are adjacent
	//this is needed since centermeres are omitted
	private ArrayList<Integer> lstPositions = new ArrayList<Integer>();
	
	private ArrayList<Integer> lstSubDataSetId = new ArrayList<Integer>();
	//utility helper class
	private Helper helper = Helper.getHelperInstance();
	
	
	public GenomeGenerator(String paraFile){
		this.parametersFile = paraFile;
	}
	
	/**
	 * do initialization for the program
	 * all global arrays must be initialized in this function
	 */
	private void initialization() throws Exception{

		try{
			//must read the parameter file first to obtain the input data file
			readParameters(parametersFile);
			
			CHR_UPPER_BOUND_ID = new int[NBR_OF_CHR + 1];			
			POS_MAX_DIST_WEIGHT = new double[NBR_OF_CHR + 1][NBR_OF_CHR + 1];
			POS_MIN_DIST_WEIGHT = new double[NBR_OF_CHR + 1][NBR_OF_CHR + 1];
			NEG_MIN_DIST_WEIGHT = new double[NBR_OF_CHR + 1][NBR_OF_CHR + 1];
			NEG_MAX_DIST_WEIGHT = new double[NBR_OF_CHR + 1][NBR_OF_CHR + 1];
			if (VERBOSE){
				System.out.println("\nReading parameters for genome..........................");
			}
			readGenomeParameters();
			
			if (VERBOSE){
				System.out.println("\nFinish reading parameters for genome.");
			}
			
			if (VERBOSE){
				System.out.println("\nReading contact data, this probably takes some time !...........................");				
			}
			
			
		}catch(Exception e){
			e.printStackTrace();
			System.err.println("No structure is generated !");
			throw e;
		}

		//create a new log file for this run, name is year_month_day_hour_minute_second
		DateFormat dateFormat = new SimpleDateFormat("yyyy_MM_dd_HH_mm_ss");
		Date date = new Date();
		String logFileName = FILE_PREFIX + "_log_" + dateFormat.format(date) + ".txt";
		try {
			logPW = new PrintWriter(OUTPUT_FOLDER + "/" + logFileName);
		} catch (FileNotFoundException e) {
			System.err.println("\nCannot create log file !");
			e.printStackTrace();
		}

		long startTime = System.currentTimeMillis();
		
		readInputData();
		
		if (VERBOSE){
			System.out.println("\nFinish reading contact data.");
			System.out.printf("\nReading files takes: %1$d seconds \n", (System.currentTimeMillis() - startTime) / 1000);
		}

		
		//allocate the array for the structure, right after getting n
		str = new double[n * 3];
		
		//determine if the structure should be scaled down for optimization
		if (POS_MAX_DIST > Constants.LARGE_DISTANCE_FOR_TANH){
			zoomFactor = POS_MAX_DIST / Constants.SCALE_DISTANCE;
			
			POS_MAX_DIST = Constants.SCALE_DISTANCE;
			
			ADJACENT_DIST /= zoomFactor;
			POS_MIN_DIST /= zoomFactor;
			NEG_MAX_DIST_INTRA /= zoomFactor;
			NEG_MAX_DIST_INTER /= zoomFactor;
			NEG_MIN_DIST /= zoomFactor;
			
		}
		
		//get the number of processor available
		int numOfcores = Runtime.getRuntime().availableProcessors();
		
		if (numOfcores == 0){
			numOfcores = 2;// default number when this parameter cannot be detected
		}
		//each core can take care of 2 threads		
		//limit number of threads to avoid excessive communication cost
		num_of_processors = Math.min(numOfcores * 2 , Constants.MAX_NUM_THREAD);
		System.out.println("Number of processors:" + num_of_processors);
		//divide the set of points into equal subsets, each will be processed by one processor (thread)

		String inEclipseStr = System.getProperty("runInEclipse");
		if ("true".equalsIgnoreCase(inEclipseStr)){
			numOfcores = 1;
		}
		
		helper.divideDataSet(n, num_of_processors, lstSubDataSetId);
		
		//reassign num_of_processors according to the division
		num_of_processors = lstSubDataSetId.size();		
		
		normalizeContactMT();
		
		//log parameters used
		logParameters();
	}
	
	/**
	 * main flow to generate structures
	 * @throws InterruptedException 
	 */
	public void generateStructure() throws IOException, Exception{
	
		//perform all steps needed to run the program
		initialization();
		

		StructureEvaluationVO strVO;
		StructureScoring strScoring = new StructureScoring();
		String fileName;

		GradientAscent gradientAscent = new GradientAscent(this, str, VERBOSE);
		//if the learning rate 
		if (LEARNING_RATE != 0){
			gradientAscent.setInitialLearingRate(LEARNING_RATE);
		}
		
		long startTime = System.currentTimeMillis();
		int i = 0;
		while(i < NUM){
			
			initializeStructure();

			gradientAscent.performGradientAscent(max_iteration);
			
			if (zoomFactor > 1.0){
				//zoom out a factor = square root of zoomFactor, because all distances are square
				helper.zoomStructure(str, Math.sqrt(zoomFactor));
			}
			
			strVO = strScoring.calculateScore(str, contactMT,POS_MAX_DIST * zoomFactor,idToChr, NBR_OF_CHR);
			//only print structure with scores larger than 1%
//			if (strVO.getTotalContactScore() > 1.0 && strVO.getTotalNonContactScore() > 1.0){
				
				fileName = FILE_PREFIX + "_" + System.currentTimeMillis() + ".pdb" ;
				helper.writeStructure(OUTPUT_FOLDER + "/" + fileName ,str, idToChr, FILE_HEADER);			
				
				logEvaluation(strVO,OUTPUT_FOLDER + "/" + fileName.replace(".pdb", ""));
				i++;
//			}

		}
		
		if (VERBOSE){			
			System.out.printf("\n%1$d structures generated in: %2$.2f minutes \n", NUM, (System.currentTimeMillis() - startTime) / (1000.0 * 60));
		}

		
		if (logPW != null){
			logPW.printf("\n%1$d structures generated in: %2$.2f minutes \n", NUM, (System.currentTimeMillis() - startTime) / (1000.0 * 60));
			logPW.close();
		}
		
	}

	/**
	 * read contact matrix into contactMT[]
	 * 
	 */
	private void readInputData() throws Exception{

		double intra_thres = is_intra_thres_percent ? 0 : INTRA_IF_THRESHOLD;
		double inter_thres = is_inter_thres_percent ? 0 : INTER_IF_THRESHOLD;
		
		try{
						
			contactMT = helper.readContactData(INPUT_FILE,intra_thres,inter_thres,CHR_UPPER_BOUND_ID,lstPositions,idToChr);
			
			for(int i = 0; i < lstPositions.size(); i++){
				if (lstPositions.get(i) != i + 1){
					System.out.println("");
				}
			}
			
		}catch(Exception ex){
			ex.printStackTrace();
			throw ex;
		}
		
		if (null == contactMT){
			throw new Exception("Input data is not available, please check input file!");
		}
		
		n = contactMT.getSize();
				
		
		//filter if thres are in percent
		if (is_intra_thres_percent){
			if (INTRA_IF_THRESHOLD > 100 || INTRA_IF_THRESHOLD < 0){
				throw new Exception("Invalid INTRA_IF_THRESHOLD ! ");
			}
						
			ArrayList<Double> lst = new ArrayList<Double>();
			for(int i = 0; i < n; i++){
				for(int j = i + 1; j < n; j++){
					if (idToChr.get(i) == idToChr.get(j)){
						lst.add(contactMT.getContact(i,j));
					}
				}
			}
			Collections.sort(lst);
			intra_thres = lst.get((int)(lst.size() * (100 - INTRA_IF_THRESHOLD) / 100));
			lst = null;
			
			
			if (logPW != null){
				logPW.println("Intra-chrosomal threshold computed from the percentage: " + intra_thres);
			}
			
			//filter
			for(int i = 0; i < n; i++){
				for(int j = i + 1; j < n; j++){
					if (idToChr.get(i) == idToChr.get(j) && contactMT.getContact(i, j) < intra_thres){
						contactMT.updateContact(i, j, 0.0);
					}
				}
			}			
		}
		
		int totalPossible = 0, totalIntra = 0;
		for(int i = 0; i < n; i++){
			for(int j = i + 1; j < n; j++){
				if (idToChr.get(i) == idToChr.get(j)){
					totalPossible++;
					if (contactMT.getContact(i,j) < intra_thres){
						totalIntra++;
					}
				}
			}
		}
		INTRA_PERCENT = totalIntra * 100 / totalPossible;
		
		
		if (is_inter_thres_percent){
			
			if (INTER_IF_THRESHOLD > 100 || INTER_IF_THRESHOLD < 0){
				throw new Exception("Invalid INTER_IF_THRESHOLD ! ");
			}
			
			ArrayList<Double> lst = new ArrayList<Double>();
			for(int i = 0; i < n; i++){
				for(int j = i + 1; j < n; j++){
					if (idToChr.get(i) != idToChr.get(j)){
						lst.add(contactMT.getContact(i,j));
					}
				}
			}
			Collections.sort(lst);
			inter_thres = lst.get((int)(lst.size() * (100 - INTER_IF_THRESHOLD) / 100));
			lst = null;
			INTER_IF_THRESHOLD = inter_thres;
			
			if (logPW != null){
				logPW.println("Inter-chrosomal threshold computed from the percentage: " + inter_thres);
			}
			
			//filter
			for(int i = 0; i < n; i++){
				for(int j = i + 1; j < n; j++){
					if (idToChr.get(i) != idToChr.get(j) && contactMT.getContact(i, j) < inter_thres){
						contactMT.updateContact(i, j, 0.0);
					}
				}
			}			
		}
		
		

		
	}
	
	/**
	 * assign max interaction frequency (max IF) to the contacts of adjacent points
	 * to reduce the effect of noise, max IF = the average of three largest interaction frequencies
	 */
	private void normalizeContactMT(){
		

		double maxIF,max1 = 0, max2 = 0, max3 = 0, max4 = 0, max5 = 0, max6 = 0, max7 = 0;
		double f,avgIF = 0;
		long count = 0;
		//find the largest interaction frequency (IF)
		for(int i = 0; i < n; i++){
			for(int j = i + 1; j < n; j++){
				f = contactMT.getContact(i, j);
				//f = contactMT[i][j];
								
				if (f > max1){	
					max7 = max6;
					max6 = max5;
					max5 = max4;
					max4 = max3;
					max3 = max2;
					max2 = max1;
					max1 = f;
				}else if (f > max2){
					max7 = max6;
					max6 = max5;
					max5 = max4;
					max4 = max3;					
					max3 = max2;
					max2 = f;
				}else if (f > max3){
					max7 = max6;
					max6 = max5;
					max5 = max4;
					max4 = max3;
					max3 = f;
				}else if (f > max4){
					max7 = max6;
					max6 = max5;
					max5 = max4;
					max4 = f;
				}else if (f > max5){
					max7 = max6;
					max6 = max5;
					max5 = f;
				}else if (f > max6){
					max7 = max6;
					max6 = f;
				}else if (f > max7){
					max7 = f;
				}
				
			}
		}
		


		
		//in poor dataset, max1 could be thousands time larger than max2, which causes problem
		if (max1 > max2 * 100){
			max1 = max2;			
		}
		
		maxIF = (max1 + max2 + max3 + max4) / 4.0 ; //because all IFs are divided by avgIF

		//maxIF = max1;
		
		for(int i = 1; i < n - 1; i++){
			if ((idToChr.get(i) == idToChr.get(i - 1) /*&& Math.abs(lstPositions.get(i) - lstPositions.get(i - 1)) == 1*/)){
				contactMT.updateContact(i, i - 1, maxIF);
				//contactMT[i][i - 1] = maxIF;
				//contactMT[i - 1][i] = maxIF;
			}
		}

		
		count = contactMT.getNbrContact();
		totalIF = contactMT.getTotalIF();
		avgIF = totalIF / count;
		
		if (logPW != null){
			logPW.println("Number of contact:" + count);
			logPW.println("Sum of all IFs:" + totalIF);
		}
		
		normTotalIF = avgIF;
		//changed on Mar 11 because yeast genome has some problem
		//totalIF /= n;

	}
	/**
	 * write statistics of the structure 'fileName' 
	 * @param vo
	 * @param fileName
	 * @throws FileNotFoundException 
	 */
	private void logEvaluation(StructureEvaluationVO vo, String fileName) throws FileNotFoundException{
		String lineSeparator = "------------------------------------------------------------------------------------------------------------";
		PrintWriter pw = new PrintWriter(fileName + "_evaluation.txt");
		
		double[][] contactScore = vo.getContactScore();
		double[][] nonContactScore = vo.getNonContactScore();
		double[][] satisfiedIFPercent = vo.getSatisfiedIFPercent();
		double[][] contactPercent = vo.getContactPercent();
		
		//pw.println("Correlation between IFs and reconstructed distances (closer to -1 is better): " + vo.getCor());
		//	375			//pw.println("Correlation between IFs and reconstructed distances (closer to -1 is better): " + vo.getCor()); 

		pw.printf("Total contact score: %.2f, total non-contact score: %.2f \n",vo.getTotalContactScore(),vo.getTotalNonContactScore()); 	 	
	
		pw.printf("%1$-12s | %2$-12s | %3$-15s | %4$-20s | %5$-30s | %6$-20s\n","chrosome 1","chromosome 2","Contact Score","Non-Contact Score","Percent of Satisfied IF","Contact Percent" );
		pw.println(lineSeparator);
		for(int i = 1; i <= NBR_OF_CHR; i++){
			for(int j = i; j <= NBR_OF_CHR; j++){
				pw.printf("%1$-12d | %2$-12d | %3$-15.2f | %4$-20.2f | %5$-30.2f | %6$-20.2f\n",i,j,
						contactScore != null ? contactScore[i][j] : -1.0,
						nonContactScore != null ? nonContactScore[i][j] : -1.0,
						satisfiedIFPercent != null ? satisfiedIFPercent[i][j] : -1.0,
						contactPercent != null ? contactPercent[i][j] : -1.0);
				pw.println(lineSeparator);
			}
		}
		
		pw.close();
	}

	
	/**
	 * read parameters necessary for the genome structure
	 * @throws Exception
	 */
	private void readGenomeParameters() throws Exception{
		//upper bound for chromosomes is needed only when there are more than one chromosome
		if (NBR_OF_CHR > 1){
			helper.readArray(CHR_UPPER_BOUND_ID_FILE, CHR_UPPER_BOUND_ID);
		}
		helper.readMatrix(POS_MAX_DIST_WEIGHT_FILE, POS_MAX_DIST_WEIGHT);
		helper.readMatrix(POS_MIN_DIST_WEIGHT_FILE, POS_MIN_DIST_WEIGHT);
		helper.readMatrix(NEG_MIN_DIST_WEIGHT_FILE, NEG_MIN_DIST_WEIGHT);
		helper.readMatrix(NEG_MAX_DIST_WEIGHT_FILE, NEG_MAX_DIST_WEIGHT);
		
	}
	
	private void logParameters(){
		if (logPW == null){
			return;
		}
		
		logPW.println("Number of chromosomes: " + NBR_OF_CHR);
		logPW.println("Input file: " + INPUT_FILE);
		logPW.println("Intra-chromosomal threshold: " + INTRA_IF_THRESHOLD);
		if (NBR_OF_CHR > 1){
			logPW.println("Inter-chromosomal threshold: " + INTER_IF_THRESHOLD);
		}
		logPW.println("Contact distance: " + POS_MAX_DIST);
		logPW.println("Minimum distance: " + POS_MIN_DIST);
		logPW.println("Maximum distance of 2 adjacent fragments: " + ADJACENT_DIST);
		logPW.println("Maximum distance between any 2 fragments of the same chromosome: " + NEG_MAX_DIST_INTRA);
		if (NBR_OF_CHR > 1){
			logPW.println("Maximum distance between any 2 fragments: " + NEG_MAX_DIST_INTER);
		}
		
		logPW.println("Learning rate: " + LEARNING_RATE);
		logPW.println("Maximum number of iteration for the optimization: " + max_iteration);
		
		logPW.println("Weights to enforce contacts to be satisfied (POS_MAX_DIST_WEIGHT_FILE): " + POS_MAX_DIST_WEIGHT_FILE);
		logPW.println("Weight to enfore minimum distance between fragments (POS_MIN_DIST_WEIGHT_FILE) : " + POS_MIN_DIST_WEIGHT_FILE);
		logPW.println("Weight to enfore non-contacts to be satisfied (NEG_MIN_DIST_WEIGHT_FILE) : " + NEG_MIN_DIST_WEIGHT_FILE);
		logPW.println("Weight to enfore maximum distance between fragments (NEG_MAX_DIST_WEIGHT_FILE) : " + NEG_MAX_DIST_WEIGHT_FILE);
		
		logPW.flush();
	}
	
	/**
	 * Read parameters for the program
	 * @param genomeParaFile
	 */
	private void readParameters(String paraFile) throws Exception{
		File file = new File(paraFile);
		FileReader fr = null;
		BufferedReader br = null;
		String ln;
		String[] st;
		int count = 0;
		Pattern splitRegex = Pattern.compile("[=\\s#]+");
		
		try{
		fr = new FileReader(file);
		br = new BufferedReader(fr);
		
		while((ln = br.readLine()) != null){
			if (ln.startsWith("#")){
				continue;
			}
			
			st = splitRegex.split(ln);
			if (st[0].equalsIgnoreCase(Constants.NUM_KEY)){
				NUM = Integer.parseInt(st[1]);
				count++;
				
			}else if (st[0].equalsIgnoreCase(Constants.NBR_OF_CHR_KEY)){
				NBR_OF_CHR = Integer.parseInt(st[1]);
				count++;
				
			}else if (st[0].equalsIgnoreCase(Constants.ADJACENT_DIST_KEY)){
				ADJACENT_DIST = Double.parseDouble(st[1]);
				count++;
				
			}else if (st[0].equalsIgnoreCase(Constants.CONTACT_DIST_KEY)){
				POS_MAX_DIST = Double.parseDouble(st[1]);
				NEG_MIN_DIST = POS_MAX_DIST;
				count++;
				
			}else if (st[0].equalsIgnoreCase(Constants.POS_MIN_DIST_KEY)){
				POS_MIN_DIST = Double.parseDouble(st[1]);
				count++;
				
			}else if (st[0].equalsIgnoreCase(Constants.NEG_MAX_DIST_INTRA_KEY)){
				NEG_MAX_DIST_INTRA = Double.parseDouble(st[1]);
				count++;
				
			}else if (st[0].equalsIgnoreCase(Constants.NEG_MAX_DIST_INTER_KEY)){
				NEG_MAX_DIST_INTER = Double.parseDouble(st[1]);
				count++;
				
			}else if (st[0].equalsIgnoreCase(Constants.POS_MAX_DIST_WEIGHT_FILE_KEY)){
				POS_MAX_DIST_WEIGHT_FILE = "";
				for (int i = 1; i < st.length; i++){
					POS_MAX_DIST_WEIGHT_FILE += st[i];
				}
				count++;
				
			}else if (st[0].equalsIgnoreCase(Constants.POS_MIN_DIST_WEIGHT_FILE_KEY)){
				POS_MIN_DIST_WEIGHT_FILE = "";
				for (int i = 1; i < st.length; i++){
					POS_MIN_DIST_WEIGHT_FILE += st[i];
				}
				count++;
				
			}else if (st[0].equalsIgnoreCase(Constants.NEG_MIN_DIST_WEIGHT_FILE_KEY)){
				NEG_MIN_DIST_WEIGHT_FILE = "";
				for (int i = 1; i < st.length; i++){
					NEG_MIN_DIST_WEIGHT_FILE += st[i];
				}
				count++;
				
			}else if (st[0].equalsIgnoreCase(Constants.NEG_MAX_DIST_WEIGHT_FILE_KEY)){
				NEG_MAX_DIST_WEIGHT_FILE = "";
				for (int i = 1; i < st.length; i++){
					NEG_MAX_DIST_WEIGHT_FILE += st[i];
				}
				count++;
				
			}else if (st[0].equalsIgnoreCase(Constants.INTRA_IF_THRESHOLD_KEY)){
				if (!st[1].contains("%")){
					INTRA_IF_THRESHOLD = Double.parseDouble(st[1]);
					is_intra_thres_percent = false;
				}else {
					INTRA_IF_THRESHOLD = Double.parseDouble(st[1].replace("%", ""));
					is_intra_thres_percent = true;
				}
				count++;
			
			}else if (st[0].equalsIgnoreCase(Constants.INTER_IF_THRESHOLD_KEY)){
				if (!st[1].contains("%")){
					INTER_IF_THRESHOLD = Double.parseDouble(st[1]);
					is_inter_thres_percent = false;
				}else{
					INTER_IF_THRESHOLD = Double.parseDouble(st[1].replace("%", ""));
					is_inter_thres_percent = true;
				}
				count++;				
				
			}else if (st[0].equalsIgnoreCase(Constants.OUTPUT_FOLDER_KEY)){
				OUTPUT_FOLDER = "";
				for (int i = 1; i < st.length; i++){
					OUTPUT_FOLDER += st[i];
				}					
				count++;
				
//			}else if (st[0].equalsIgnoreCase(Constants.FILE_PREFIX_KEY)){
//				FILE_PREFIX = st[1];
//				count++;
//				
//			}else if (st[0].equalsIgnoreCase(Constants.FILE_HEADER_KEY)){
//				FILE_HEADER = "";
//				for (int i = 1; i < st.length; i++){
//					FILE_HEADER += st[i];
//				}				
//				count++;
				
			}else if (st[0].equalsIgnoreCase(Constants.INPUT_FILE_KEY)){
				INPUT_FILE = "";
				for (int i = 1; i < st.length; i++){
					INPUT_FILE += st[i];
				}
				String[] tmp = INPUT_FILE.split("[\\/ \\. \\\\]");
				if (INPUT_FILE.contains(".")){
					FILE_PREFIX = tmp[tmp.length - 2];
				}else{
					FILE_PREFIX = tmp[tmp.length - 1];
				}
				
				count++;				
			}else if (st[0].equalsIgnoreCase(Constants.VERBOSE_KEY)){
				VERBOSE = Boolean.parseBoolean(st[1]);
				count++;
				
			}else if (st[0].equalsIgnoreCase(Constants.CHR_UPPER_BOUND_ID_FILE_KEY)){
				CHR_UPPER_BOUND_ID_FILE = "";
				for (int i = 1; i < st.length; i++){
					CHR_UPPER_BOUND_ID_FILE += st[i];
				}			
				count++;
			}else if (st[0].equalsIgnoreCase(Constants.LEARNING_RATE_KEY)){
				LEARNING_RATE = Double.parseDouble(st[1]);
				
				count++;
			}else if (st[0].equalsIgnoreCase(Constants.MAX_ITERATION_KEY)){
				max_iteration = Integer.parseInt(st[1]);
				count++;
			}
			
			
		}
		//when the input is a chromosome, INTER_IF_THRESHOLD && CHR_UPPER_BOUND_ID_FILE are not needed,
		
		if ((count < Constants.NUM_PARAMETERS && NBR_OF_CHR > 1) || 
				(count < Constants.NUM_PARAMETERS - 2 && (!Double.isNaN(INTER_IF_THRESHOLD) || null != CHR_UPPER_BOUND_ID_FILE) ) ){
			throw getParameterException("Missing some parameters, please check the parameter file!");
		}

		}catch(Exception e){
			
			e.printStackTrace();
			throw e;			
			
		}finally{
			if (br != null){
				br.close();
			}
			if (fr != null){
				fr.close();
			}
		}
	}
	/**
	 * Initialize the genome structure, adjacent points are initialized to be closer together than the others
	 */
	private void initializeStructure(){
		
		double chrX=0,chrY=0,chrZ=0,size = 0.1;
		
		for(int i = 0; i < n; i++){
			
			//reset starting point for every chromosome
			if (i == 0 || idToChr.get(i) != idToChr.get(i - 1) ){			
				chrX = Math.random();
				chrY = Math.random();
				chrZ = Math.random();
				
			}else if (idToChr.get(i) == idToChr.get(i - 1)){
				//extend in X,Y,Z coordinate
				if (NBR_OF_CHR > 1 || INTRA_PERCENT < 60){ // if there is inter-chromosomal contact or too many non-contact already
					chrX += size * (Math.random() - 0.5);
					chrY += size * (Math.random() - 0.5);
					chrZ += size * (Math.random() - 0.5);
				}else{
					chrX += (Math.random() - 0.5) ;
					chrY += (Math.random() - 0.5) ;
					chrZ += (Math.random() - 0.5) ;
				}
			}

			str[i * 3] = chrX;
			str[i * 3 + 1] = chrY;
			str[i * 3 + 2] = chrZ;	
		}

	}
	
	/**
	 * 
	 * @param x : current values for variables
	 * @param change : gradient vector
	 * @return value of the objective function
	 * @throws InterruptedException 
	 */
	@Override
	public double calGradientAndObjective(double[] x, double[] change) throws InterruptedException{
		
		double cost = 0;
		
		GradientCaculator[] gradCalculator = new GradientCaculator[lstSubDataSetId.size()];
		
		//initialize threads
		gradCalculator[0] = new GradientCaculator(x, 0 , lstSubDataSetId.get(0), change != null);
		for(int i = 1; i < lstSubDataSetId.size(); i++){
			gradCalculator[i] = new GradientCaculator(x,  lstSubDataSetId.get(i - 1) + 1, lstSubDataSetId.get(i), change != null);
		}
		
		//start threads
		for(int i = 0; i < gradCalculator.length; i++){
			gradCalculator[i].start();
		}
		
		//wait for all threads to finish
		try {
			for(int i=0; i< gradCalculator.length; i++){
				gradCalculator[i].join();
			}			
		} catch (InterruptedException e) {			
			e.printStackTrace();
			throw e;
		}

		//aggregate the cost
		for(int i = 0; i < gradCalculator.length; i++){
			cost += gradCalculator[i].getCost();
		}
	
		if (change != null){
			//aggregate the gradient
			for(int i = 0; i < change.length; i++){
				change[i] = 0;
				for(int k = 0; k < gradCalculator.length; k++){
					change[i] += gradCalculator[k].getChange()[i];
				}
				
			}
		}
		
		return cost;
	}
	
	/**
	 * calculate the objective function only, used in line search for step size
	 */
	@Override
	public double calObjective(double[] x) throws InterruptedException {		
		return calGradientAndObjective(x, null);
		
	}
	
	private Exception getParameterException(String msg){
		
		return new Exception(msg + "\nThe parameter file has the following format, the line start by # is ignored as a comment. "
				+ "Each line for a parameter contains a key=value pair"
				+ "\nNUM = the number of structures will be generated"
				+ "\nNBR_OF_CHR_KEY = the number of chromosomes"				
				+ "\nINTRA_IF_THRESHOLD = interaction frequency (IF) threshold for intrachromosomal contact, contact with IF less than this value is considered as non-contact"
				+ "\nINTER_IF_THRESHOLD = interaction frequency (IF) threshold for interchromosomal contact, contact with IF less than this value is considered as non-contact"
				+ "\nADJACENT_DIST = the square of maximum distance between 2 adjacent fragments (ex: 0.2 for fragment of 1MB)"
				+ "\nCONTACT_DIST = the square of contact distance between 2 fragments (ex: 6.0) "
				+ "\nPOS_MIN_DIST = the square of minimum distance between 2 fragments that are in contact (ex: 0.02)"
				+ "\nNEG_MAX_DIST_INTRA = the square of maximum distance between 2 fragments that are not in contact and in the same chromosome (ex: 50)"
				+ "\nNEG_MAX_DIST_INTER = the square of maximum distance between 2 fragments that are not in contact and in different chromosome (ex: 200)"				
				+ "\nPOS_MAX_DIST_WEIGH_FILE = file contains the weights for all pairs of chromosomes to enforce the maximum distance constraints of 2 fragments that are in contact"
				+ "\nPOS_MIN_DIST_WEIGHT_FILE = file contains the weights for all pairs of chromosomes to enforce the minimum distance contraints of 2 fragments that are in contact"
				+ "\nNEG_MIN_DIST_WEIGHT_FILE = file contains the weights for all pairs of chromosomes to enforce the minimum distance contraints of 2 fragments that are not in contact"
				+ "\nNEG_DIST_MAX_WEIGHT_FILE = file contains the weights for all pairs of chromosomes to enforce the maximum distance contraints of 2 fragments that are not in contact"
				+ "\nCHR_UPPER_BOUND_ID_FILE_KEY = file contains the upper bound indices of chromosomes"				
				+ "\nOUTPUT_FOLDER = the folder for generated structures"
				+ "\nFILE_PREFIX = the prefix for file name of structures"
				+ "\nFILE_HEADER = header for the structure file"
				+ "\nINPUT_FILE = file contains the contact matrix"				
				+ "\nVERBOSE = true, for information during optimization printed out"
				+ "\nLEARNING_RATE = 100, the initial learning rate for the optimization process");

	}

	
	/**
	 * This class is used to calculate the gradient for a subset of datapoints
	 * in a single threaded program, all data points are i = 1 .. n and j = i+1 ... n
	 * to calculate the gradient in parallel, one thread will be in charged of calculating for i = begin ... end
	 * 
	 * for any modification of the objective function, this function will need to be modified accordingly
	 * @author Tuan
	 *
	 */
	class GradientCaculator extends Thread{
		//the first index to calculate the gradient
		private int beg;
		//the last index to calculate the gradient
		private int end;
		//gradient to be returned
		double[] change;
		//objective function
		double cost=0;
		//indicate if calculation for gradient is needed
		boolean isGradientNeeded;
		
		double[] structure;
		
		GradientCaculator(double[] str, int b, int e, boolean isGradient){
			this.beg = b;
			this.end = e;		
			this.structure = str;
			this.isGradientNeeded = isGradient;
			
			if (isGradientNeeded){
				this.change = new double[n * 3];
			}
		}
		
		public void run(){
			double f,dist,mid,tanh1=0,tanh2=0;
			
			int chr1, chr2;
			
			for(int i = beg; i <= end; i++){
				
				 for(int j = i+1; j < n; j++){
					 
					 chr1 = idToChr.get(i);
					 chr2 = idToChr.get(j);
					 
					 dist = helper.calEuclidianDist(structure[i * 3],structure[i * 3 + 1], structure[i * 3 + 2],
							 structure[j * 3],structure[j * 3 + 1], structure[j * 3 + 2]);
					 
					 mid = 0;
					 f = contactMT.getContact(i, j);
					 //f = contactMT[i][j];
					 if (f > 0.0){
						 
						 //positive contact, the distance should be larger than POS_MIN_DIST
						 tanh2 = helper.tanh(dist - POS_MIN_DIST);
						 
						 mid += POS_MIN_DIST_WEIGHT[chr1][chr2] * (1 - tanh2 * tanh2) / normTotalIF;						 
						 cost += POS_MIN_DIST_WEIGHT[chr1][chr2] * tanh2 / normTotalIF;

						 //if (chr1 == chr2 && Math.abs(lstPositions.get(j) - lstPositions.get(i)) == 1){
						 if (chr1 == chr2 && Math.abs(j - i) == 1){
							 //these 2 points are adjacent, the distance should be less than ADJACENT_DIST
							 tanh1 = helper.tanh(ADJACENT_DIST - dist);				
							 
						 }else { 
							 //positive contact, the distance should be less than POS_MAX_DIST
							 tanh1 = helper.tanh(POS_MAX_DIST - dist);
						 }
						 mid += POS_MAX_DIST_WEIGHT[chr1][chr2] * (tanh1 * tanh1 - 1) * f / normTotalIF;
						 
						 cost += POS_MAX_DIST_WEIGHT[chr1][chr2] * tanh1 * f / normTotalIF;
						
					 }else if (f <= 0.0){
						 					 
						 if (chr1 == chr2){
							 //this is intrachromosomal non-contact, the distance should be less than NEG_MAX_DIST_INTRA
							 tanh1 = helper.tanh(NEG_MAX_DIST_INTRA - dist); 
							 							 
						 }else{
							//this is intrachromosomal non-contact, the distance should be less than NEG_MAX_DIST_INTER
							 tanh1 = helper.tanh(NEG_MAX_DIST_INTER - dist); 
						 }
						 //the distance should be larger than NEG_MIN_DIST
						 tanh2 = helper.tanh(dist - NEG_MIN_DIST);

						 mid = NEG_MAX_DIST_WEIGHT[chr1][chr2] * (tanh1 * tanh1 - 1)/normTotalIF + NEG_MIN_DIST_WEIGHT[chr1][chr2] * (1 - tanh2 * tanh2)/normTotalIF;
						 
						 cost += NEG_MAX_DIST_WEIGHT[chr1][chr2] * tanh1/normTotalIF + NEG_MIN_DIST_WEIGHT[chr1][chr2] *  tanh2/normTotalIF;
						 
					 }
					 
					 //Tuan added on Dec 16 2014
					 //To prevent 2 adjacent points too far away from each other
					 //Use log function
//					 if (j == i + 1){
//						 //the constant 0.0005 is to minimize the effect of dist > pos_max_dist
//						 //when dist < pos_max_dist, log is very negative, so that using this constant doesn't lose the expected effect
//						 cost += 0.00001 * helper.log(POS_MAX_DIST - dist);						 
//						 mid += -2 * 0.00001/(POS_MAX_DIST - dist);			
//					 }

					 if (isGradientNeeded){
						 change[i * 3 + 0] += (structure[i * 3] - structure[j * 3])  * mid;
						 change[i * 3 + 1] += (structure[i * 3 + 1] - structure[j * 3 + 1])  * mid;
						 change[i * 3 + 2] += (structure[i * 3 + 2] - structure[j * 3 + 2])  * mid;
						 
						 change[j * 3 + 0] += (structure[j * 3] - structure[i * 3])  * mid;
						 change[j * 3 + 1] += (structure[j * 3 + 1] - structure[i * 3 + 1])  * mid;
						 change[j * 3 + 2] += (structure[j * 3 + 2] - structure[i * 3 + 2])  * mid;
					 }
					 				 
				 }
				 
			}
		}

		public double[] getChange() {
			return change;
		}
		
		public double getCost() {
			return cost;
		}
		
	}


	
	public static void main(String[] args) throws Exception{
		String paraFile = "parameters_tad_res40.txt";
		if (args != null && args.length >= 1){
			paraFile = args[0];
		}
		//check if the parameter file exists in the current directory
		File file = new File(paraFile);
		if (!file.exists()){
			throw new Exception("The parameter file " + paraFile + "cannot be found in the current directory");
		}
		
		GenomeGenerator generator = new GenomeGenerator(paraFile);
		generator.generateStructure();
		
	}


	
}
