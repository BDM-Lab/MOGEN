package utility;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Scanner;
import java.util.regex.Pattern;

import valueObject.Constants;
import valueObject.Contact;
import valueObject.ContactMT;

public class Helper {
	
	private DecimalFormat df2 = new DecimalFormat("0.00");
	
	private static Helper helper = new Helper();
	
	public Helper(){		
	}
	
	public static Helper getHelperInstance(){
		if (helper == null){
			synchronized(Helper.class){
				if (helper == null){
					helper = new Helper();
				}
			}
		}
		
		return helper;
	}
	
	/**
	 * 
	 * @param inputFile
	 * @param intraContactThres
	 * @param interContactThres
	 * @param chrUpperBoundID: upper bound indices for chromosomes 
	 * @param lstPos: to be filled by this function. Store positions, it could be return to later identify original fragment of index in the structure
	 * @param idToChr : to be filled by this function. Map indices into corresponding chromosomes
	 * @return matrix of contact
	 * @throws FileNotFoundException
	 * @throws IOException
	 */
	public ContactMT readContactData(String inputFile, double intraContactThres,double interContactThres, 
			int[] chrUpperBoundID,ArrayList<Integer> lstPos, HashMap<Integer,Integer> idToChr)	throws FileNotFoundException, Exception{
		
		//contact matrix will be returned
		//double[][] a;
		ContactMT contactMT;
		
		//map positions into absolute indexs, to take care of unaligned centermeres		
		HashMap<Integer,Integer> hmPos = new HashMap<Integer,Integer>();

		
		int x,y,id1,id2,chr1,chr2;
		double f;
		int nbr = 3; // number of numbers in each line
		
		File file = new File(inputFile);

		FileReader fr=null;
		BufferedReader br = null;		
		Pattern splitRegex = Pattern.compile("[:\\s]+");
		
		HashSet<Integer> posSet = new HashSet<Integer>();
		StringBuilder sb = new StringBuilder();
		try{
				
			fr = new FileReader(file);
			br = new BufferedReader(fr);
			String ln;
			String[] st;
			int count = 1;
			ln = br.readLine();
			sb.append(ln).append(" ");
			nbr = splitRegex.split(ln.trim()).length;
			long progress = 0;
			while((ln = br.readLine()) != null){
				if (ln.startsWith("#") || ln.trim().length() == 0){
					continue;
				}
				
				//read every 10 thousand lines and split at once
				sb.append(ln).append(" ");
				count++;
				
				if (count == 200000){
					count = 0;
					st = splitRegex.split(sb.toString());
					sb = new StringBuilder();
					
					if (st.length % nbr != 0){
						throw new Exception("There is a line that doesn't contain exactly 3 numbers");
					}
					//each line contains 'nbr' numbers
					for(int i = 0; i < st.length / nbr; i++){
						
						//only take first 3 numbers
						x = Integer.parseInt(st[i * nbr + 0]);						
						//position2
						y = Integer.parseInt(st[i * nbr + 1]);
						
						//interaction frequency
						//f = Double.parseDouble(st[i * nbr + 2]);
						
						//keeping absolute positions, so that later they can be recovered from indices
						
						posSet.add(x);
						
						posSet.add(y);
						 
					}
				}
				progress++;
				System.out.println(progress * 200000 + " input lines have been read !");
			}
			br.close();
			fr.close();
			
			//if sb is not empty
			if (sb.toString().trim().length() > 0){
				st = splitRegex.split(sb.toString());
				sb = new StringBuilder();
				
				if (st.length % nbr != 0){
					throw new Exception("There is a line that doesn't contain exactly 3 numbers");
				}
				//each line contains 3 numbers
				for(int i = 0; i < st.length / nbr; i++){
					x = Integer.parseInt(st[i * nbr + 0]);
					
					//position2
					y = Integer.parseInt(st[i * nbr + 1]);
					
					//interaction frequency
					//f = Double.parseDouble(st[i * 3 + 2]);
					
					//keeping absolute positions, so that later they can be recovered from indices
					
					posSet.add(x);
					
					posSet.add(y);
					 
				}
				
			}

			lstPos.addAll(posSet);
			
			//sort the lst of position ascendingly
			Collections.sort(lstPos);			
			
			//map postion into absolute index
			for(int i = 0; i < lstPos.size(); i++){
				hmPos.put(lstPos.get(i), i);
				int chr = determineChr(lstPos.get(i),chrUpperBoundID);
				if (chr == -1){
					throw new Exception("Input file is not consistent, cannot determine chromosome for this index:" + lstPos.get(i) + ".\n"
							+ "Make sure the upper bound indices for chromosomes are correct!");
				}
				idToChr.put(i,chr);
				
			}
			
			//initialize the matrix contact
			//a = new double[lstPos.size()][lstPos.size()];
			contactMT = ContactMT.getContactMTInstance(lstPos.size());
			
			fr = new FileReader(file);
			br = new BufferedReader(fr);
			count = 0;
			progress = 0;
			while((ln = br.readLine()) != null){
				
				if (ln.startsWith("#") || ln.trim().length() == 0){
					continue;
				}
				
				//read every a thoudsand lines and split at once
				sb.append(ln).append(" ");
				count++;
				
				if (count == 200000){
					count = 0;
					st = splitRegex.split(sb.toString());
					sb = new StringBuilder();
					
					if (st.length % nbr != 0){
						throw new Exception("There is a line that doesn't contain exactly 3 numbers");
					}
					//each line contains 3 numbers
					for(int i = 0; i < st.length / nbr; i++){
						x = Integer.parseInt(st[i * nbr + 0]);
						
						//position2
						y = Integer.parseInt(st[i * nbr + 1]);
						
						//interaction frequency
						f = Double.parseDouble(st[i * nbr + 2]);
						
						//keeping absolute positions, so that later they can be recovered from indices
						
						id1 = hmPos.get(x);
						id2 = hmPos.get(y);
						chr1 = idToChr.get(id1);
						chr2 = idToChr.get(id2);
						
						//if the frequency is less than thresholds, ignore the contact so that it is considered as non-contact
						if ((chr1 == chr2 && f < intraContactThres) 
								|| (chr1 != chr2 && f < interContactThres)){
							continue;
						}

						//a[id1][id2] = f;
						//a[id2][id1] = f;
						contactMT.addContact(id1, id2, f);
						 
					}
					
					progress++;
					System.out.println(progress * 200000 + " input lines have been read - second time!");
				}	

			}

			//if sb is not empty
			if (sb.toString().trim().length() > 0){
				st = splitRegex.split(sb.toString());
				sb = new StringBuilder();
				
				if (st.length % nbr != 0){
					throw new Exception("There is a line that doesn't contain exactly 3 numbers");
				}
				//each line contains 3 numbers
				for(int i = 0; i < st.length / nbr; i++){
					x = Integer.parseInt(st[i * nbr + 0]);
					
					//position2
					y = Integer.parseInt(st[i * nbr + 1]);
					
					//interaction frequency
					f = Double.parseDouble(st[i * nbr + 2]);
					
					//keeping absolute positions, so that later they can be recovered from indices
					
					id1 = hmPos.get(x);
					id2 = hmPos.get(y);
					chr1 = idToChr.get(id1);
					chr2 = idToChr.get(id2);
					
					//if the frequency is less than thresholds, ignore the contact so that it is considered as non-contact
					if ((chr1 == chr2 && f < intraContactThres) 
							|| (chr1 != chr2 && f < interContactThres)){
						continue;
					}

					//a[id1][id2] = f;
					//a[id2][id1] = f;
					contactMT.addContact(id1, id2, f);

					 
				}
				
			}

		
		}catch(FileNotFoundException ex){
			ex.printStackTrace();
			throw ex;
		}catch(IOException ex){
			ex.printStackTrace();
			throw ex;			
		}finally{
			if (br != null){
				br.close();
			}
			if (fr != null){
				fr.close();
			}
			
		}
		System.out.println("Sorting contact ...");
		contactMT.sortContact();
		System.out.println("Done with reading input !");
		return contactMT;
	}
	
	
	/**
	 * 
	 * @param inputFile
	 * @param intraContactThres
	 * @param interContactThres
	 * @param chrUpperBoundID: upper bound indices for chromosomes 
	 * @param lstPos: to be filled by this function. Store positions, it could be return to later identify original fragment of index in the structure
	 * @param idToChr : to be filled by this function. Map indices into corresponding chromosomes
	 * @return a map with key is pos1 * n + pos2 and value is the object with pos1,pos2 (pos1, pos2 are indices of the 3D model without centremere and not of the chromosome)
	 * @throws FileNotFoundException
	 * @throws IOException
	 */
	public HashMap<Integer,Contact> readContactListAsList(String inputFile, double intraContactThres,double interContactThres, 
			int[] chrUpperBoundID,ArrayList<Integer> lstPos, HashMap<Integer,Integer> idToChr)	throws FileNotFoundException, Exception{
		
		//contact matrix will be returned
		//double[][] a;
		//List<Contact> lst = new ArrayList<Contact>();
		HashMap<Integer,Contact> map = new HashMap<Integer,Contact>();
		
		//map positions into absolute indexs, to take care of unaligned centermeres		
		HashMap<Integer,Integer> hmPos = new HashMap<Integer,Integer>();

		
		int x,y,id1,id2,chr1,chr2;
		double f;
		int n, nbr = 3; // number of numbers in each line
		
		File file = new File(inputFile);

		FileReader fr=null;
		BufferedReader br = null;		
		Pattern splitRegex = Pattern.compile("[:\\s]+");
		
		HashSet<Integer> posSet = new HashSet<Integer>();
		StringBuilder sb = new StringBuilder();
		try{
				
			ArrayList<Contact> ctLst = new ArrayList<Contact>();
			
			fr = new FileReader(file);
			br = new BufferedReader(fr);
			String ln;
			String[] st;
			int count = 1;
			ln = br.readLine();
			sb.append(ln).append(" ");
			nbr = splitRegex.split(ln.trim()).length;
			//long progress = 0;
			while((ln = br.readLine()) != null){
				if (ln.startsWith("#") || ln.trim().length() == 0){
					continue;
				}
				
				//read every 10 thousand lines and split at once
				sb.append(ln).append(" ");
				count++;
				
				if (count == 200000){
					count = 0;
					st = splitRegex.split(sb.toString());
					sb = new StringBuilder();
					
					if (st.length % nbr != 0){
						throw new Exception("There is a line that doesn't contain exactly 3 numbers");
					}
					//each line contains 'nbr' numbers
					for(int i = 0; i < st.length / nbr; i++){
						
						//only take first 3 numbers
						x = Integer.parseInt(st[i * nbr + 0]);						
						//position2
						y = Integer.parseInt(st[i * nbr + 1]);
						
						//interaction frequency
						f = Double.parseDouble(st[i * nbr + 2]);
						
						//keeping absolute positions, so that later they can be recovered from indices
						
						//if (x > 10000 || y > 10000) continue;	
						
						posSet.add(x);
						
						posSet.add(y);						
											
						//keeping absolute positions, so that later they can be recovered from indices
						
						//if the frequency is less than thresholds, ignore the contact so that it is considered as non-contact


						if (x != y){
							//lst.add(new Contact(id1,chr1, id2,chr2, f));
							Contact ct = new Contact(x,1, y, 1, f); //this is dummy and will be replaced later
							ctLst.add(ct);
						}
						 
					}
				}
				//progress++;
				//System.out.println(progress * 200000 + " input lines have been read !");
			}
			br.close();
			fr.close();
			
			//if sb is not empty
			if (sb.toString().trim().length() > 0){
				st = splitRegex.split(sb.toString());
				sb = new StringBuilder();
				
				if (st.length % nbr != 0){
					throw new Exception("There is a line that doesn't contain exactly 3 numbers");
				}
				//each line contains 3 numbers
				for(int i = 0; i < st.length / nbr; i++){
					x = Integer.parseInt(st[i * nbr + 0]);
					
					//position2
					y = Integer.parseInt(st[i * nbr + 1]);
					
					//interaction frequency
					f = Double.parseDouble(st[i * 3 + 2]);
					
					//keeping absolute positions, so that later they can be recovered from indices
					//if (x > 10000 || y > 10000) continue;
					
					posSet.add(x);
					
					posSet.add(y);
					
					if (x != y){
						//lst.add(new Contact(id1,chr1, id2,chr2, f));
						Contact ct = new Contact(x,1, y, 1, f); //this is dummy and will be replaced later
						ctLst.add(ct);
					}
				}
				
			}

			lstPos.addAll(posSet);
			
			//sort the lst of position ascendingly
			Collections.sort(lstPos);
			
			//Number of points
			n = lstPos.size();
			
			//map postion into absolute index
			for(int i = 0; i < lstPos.size(); i++){
				hmPos.put(lstPos.get(i), i);
				int chr = determineChr(lstPos.get(i),chrUpperBoundID);
				if (chr == -1){
					throw new Exception("Input file is not consistent, cannot determine chromosome for this index:" + lstPos.get(i) + ".\n"
							+ "Make sure the upper bound indices for chromosomes are correct!");
				}
				idToChr.put(i,chr);
				
			}
			
			
			//correct pos in contacts
			for(Contact ct: ctLst){
				ct.setPos1(hmPos.get(ct.getPos1()));
				ct.setPos2(hmPos.get(ct.getPos2()));
				
				ct.setChr1(idToChr.get(ct.getPos1()));
				ct.setChr2(idToChr.get(ct.getPos2()));
				
				f = ct.getIF();
				chr1 = ct.getChr1();
				chr2 = ct.getChr2();
				
				if ((chr1 == chr2 && f < intraContactThres) 
						|| (chr1 != chr2 && f < interContactThres) 
						|| ct.getPos1() == ct.getPos2()){					
					//nothing
				}else{
					
					map.put(ct.getPos1() * n + ct.getPos2(), ct);
					
				}
			}
		
		}catch(FileNotFoundException ex){
			ex.printStackTrace();
			throw ex;
		}catch(IOException ex){
			ex.printStackTrace();
			throw ex;			
		}finally{
			if (br != null){
				br.close();
			}
			if (fr != null){
				fr.close();
			}
			
		}
		
		
		//System.out.println("Sorting contact ...");
		//Collections.sort(lst);
		System.out.println("Done with reading input !");
		
		return map;
	}

	
	/**
	 * Determine the corresponding chromosome for the actual id (in the genome, where centermeres are omitted)
	 * @param id: the actual index in the genome
	 * @param chrUpperBoundID: the upper bound indices of chromosomes (including), chrUpperBoundID[i] is for chromosome i
	 * @return chromosome number
	 */
	private int determineChr(int id, int[] chrUpperBoundID ){
		//if chrUpperBoundID == null, then, there may be only one chromsome
		//if chrUpperBoundID.length == 2, then, there may be only one chromsome, because the first one chrUpperBoundID[0] is ignored
		if(null == chrUpperBoundID || chrUpperBoundID.length == 2){
			return 1;
		}
		for(int i = 1; i <= chrUpperBoundID.length - 1; i++){
			if (chrUpperBoundID[i] >= id){
				return i;
			}
		}
		
		return -1;
	}
	
	/**
	 * 
	 * @param n: number of points, so total number of pairs is n * (n-1)/2
	 * @param k: number of processors
	 * @return an arraylist contains starting index for each processor
	 */
	public void divideDataSet(int n, int k, ArrayList<Integer> lstSubDataSetId){
		long t = (long) n;
		long total = t * (t - 1) / 2;
		int size = (int) (total / k);
		int covered = 0;

		for(int i = 1; i < n; i ++){
			covered += n - i;
			if (covered >= size){
				lstSubDataSetId.add(i);
				covered = 0;
			}
		}
		if (covered > 0){
			lstSubDataSetId.add(n - 1);
		}
	}
	
	/**
	 * 
	 * @param total: number of constraints
	 * @param k: number of processors
	 * @return an arraylist contains ending index for each processor
	 */
	public void divideDataSetWithTotal(int total, int k, ArrayList<Integer> lstSubDataSetId){
		
		int size = total / k;		
		for(int i = 1; i < k; i++){
			lstSubDataSetId.add(i * size);
		}
		if (total % k != 0 || k == 1){
			lstSubDataSetId.add(total - 1);
		}

	}

	
	/**
	 * the first line contain a number for pairs (i,j) where i <> j
	 * each line contains the index and value for diagonal elements of the matrix: chromosome value
	 * if there is only one chromosome (mt.length == 1), only the last number is used as the intra-chromosomal parameter
	 * @param file
	 * @param mt: the matrix to be filled and returned
	 */
	public void readMatrix(String fileName, double[][] mt) throws FileNotFoundException, Exception {
		File file = new File(fileName);
		Scanner scanner = null;
		String ln;
		String[] st;
		int id;
		double value=0.0;
		Pattern splitRegex = Pattern.compile("[:\\s]+");
		
		try{
			scanner = new Scanner(file);
			while(scanner.hasNextLine()){
				ln = scanner.nextLine();
				if (ln.startsWith("#") || ln.trim().length() == 0){
					continue;
				}
				//each line contains one or two numbers
				st = splitRegex.split(ln);
				
				//if the line contains one number
				//coefficient for pairs (i,j) where i<>j
				if (st.length == 1) {
					value = Double.parseDouble(st[0]);
					
					for(int i = 1; i < mt.length; i++){
						for(int j = 1; j < mt.length; j++){
							if (i != j){
								mt[i][j] = value;
							}
						}
					}
					
					continue;
				}
				id = Integer.parseInt(st[0]);
				value = Double.parseDouble(st[1]);				
				if (id < mt.length){
					mt[id][id] = value;
				}else{
					System.err.println("Number of chromosomes: " + (mt.length - 1) + ", while chromsome id in " + fileName + ": " + id);
					throw new Exception("In valid chromosome ID");
				}
			}
			
			//if there is only one chromosome, the last number is used
			//check mt.length == 2 because the first mt[0][0] is unused
			if (mt.length == 2){
				mt[1][1] = value;
			}
			
		}catch(Exception e){
			
			System.err.println("\n\nThe first non-comment line of the file " + fileName + " must contain one number for pairs(i,j) where i <> j \n\n");
			System.err.println("\n\nSubsequent lines of the file " + fileName + " must contain 2 numbers in format: chromosome value\n\n");		
			e.printStackTrace();			
			throw e;
			
		}finally{
			if (scanner != null){
				scanner.close();
			}
		}
	}
	/**
	 * each line contains an chromosome index and a value
	 * @param fileName
	 * @param upperBoundID
	 */
	public void readArray(String fileName, int[] upperBoundID) throws Exception{
		File file = new File(fileName);
		FileReader fr = null;
		BufferedReader br = null;
		String ln;
		String[] st;
		int id,value;
		try{
			fr = new FileReader(file);
			br = new BufferedReader(fr);
			while((ln = br.readLine()) != null){
				if (ln.startsWith("#") ||ln.trim().length() == 0){
					continue;
				}
				st = ln.split("[:\\s]+");
				
				if (st.length < 2){
					continue;
				}
				
				id = Integer.parseInt(st[0]);
				value = Integer.parseInt(st[1]);
				//check if the index is bounded
				if (id < upperBoundID.length){
					upperBoundID[id] = value;
				}else{
					System.err.println("Number of chromosomes: " + (upperBoundID.length - 1) + ", while chromsome id in " + fileName + ": " + id);
					throw new Exception("In valid chromosome ID");
				}
				
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
	 * 
	 * @param x1
	 * @param y1
	 * @param x2
	 * @param y2
	 * @return square euclidean distance of (x1,y1) and (x2,y2) 
	 */
	public double calEuclidianDist(double x1, double y1, double z1, double x2, double y2, double z2){
		return (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2);
	}
	
	/**
	 * log function to take care of the case x <= 0	
	 * @param x
	 * @return
	 */
	public double log(double x){
		if (x <= 0){
			return Integer.MIN_VALUE;
		}else{
			return Math.log(x);
		}
	}
	
	/**
	 * 
	 * @param 
	 * @return hyperbolic function tanh(x)
	 */
	public double tanh(double x){		
		
		//this is an approximation for the tanh function
		if (x > 15){
			return 1.0;
		}else if (x < -15){
			return -1;
		}
		double ex = Math.pow(Math.E,x);
		double emx = 1/ex;
		double y = (ex - emx)/(ex + emx) ;		
		
		return y;

	}
	
	/**
	 * calculate the scaling factor for x, to scale down x to the sensitive range of the tanh function
	 * @param x
	 * @return
	 */
	public double calTanhScalingFactor(double x){
		return Math.min(1.0, Constants.LARGE_DISTANCE_FOR_TANH / x);
	}
	
	/**
	 * Zoom the structure (str) by the factor
	 * @param str
	 * @param fator
	 */
	public void zoomStructure(double[] str, double factor){
		if (str == null){
			return;
		}
		for(int i = 0; i < str.length; i++){
			str[i] *= factor;
		}
	}
	
	/**
	 * Read a PDB file and return a coordinate matrix, general one
	 * @param fileName
	 * @return matrix[n][3]
	 */
	public double[][] loadPDBStructure(String fileName) throws Exception{
		File file = new File(fileName);
		FileReader fr = null;
		BufferedReader br = null;
		double[][] coor;
		ArrayList<Double> lstX = new ArrayList<Double>();
		ArrayList<Double> lstY = new ArrayList<Double>();
		ArrayList<Double> lstZ = new ArrayList<Double>();
		try{
			fr = new FileReader(file);
			br = new BufferedReader(fr);
			

			String ln;
			String[] st;
			Pattern splitRegex = Pattern.compile("[\\s]+");
			while((ln = br.readLine()) != null){
				if (! ln.startsWith("ATOM")){
					continue;
				}
				st = splitRegex.split(ln);				
				if (st.length > 7 && st[2].equals("CA")){
					lstX.add(Double.parseDouble(st[5]));
					lstY.add(Double.parseDouble(st[6]));
					lstZ.add(Double.parseDouble(st[7]));
				}
			}		
			
			
		}catch(Exception ex){
			ex.printStackTrace();
			throw ex;			
			
		}finally{
			if (br != null){
				br.close();
			}
			if (fr != null){
				fr.close();
			}
		}
		
		coor = new double[lstX.size()][3];
		for(int i = 0; i < lstX.size(); i++){
			coor[i][0] = lstX.get(i);
			coor[i][1] = lstY.get(i);
			coor[i][2] = lstZ.get(i);
		}
		
		return coor;
	}

	/**
	 * Convert a coordinate matrix into a matrix ready for output to write out PDB file
	 * @param a
	 * @return
	 */
	public double[] makeMatrixForPDBOutput(double[][] a){
		double[] str = new double[a.length * 3];
		for(int i = 0; i < a.length; i++){
			str[i * 3 + 0] = a[i][0];
			str[i * 3 + 1] = a[i][1];
			str[i * 3 + 2] = a[i][2];
		}
		
		return str;
	}

	/**
	 * To extract one chromosome from the whole genome structure
	 * @param a: coordinates of all points in the genome
	 * @param start: starting index
	 * @param end: ending index (inclusive)
	 * @return coordinates of points of the chromosome
	 */
	public double[][] extractChr(double[][] a, int start, int end){
		
		int n = end - start + 1;
		double[][] b = new double[n][n];
		for(int i = 0; i < n; i++){
			b[i][0] = a[i + start][0];
			b[i][1] = a[i + start][1];
			b[i][2] = a[i + start][2];
		}
		
		return b;		
	}

	/**
	 * Please refer to the pdb format file for detail of each column
	 * @param pathFilename: output file name
	 * @param structure: every 3 consecutive points is one point in 3D
	 * @param idToChr: to identify if 2 fragments belong to the same chromosome
	 * @param header for the pdb file
	 */
	public void writeStructure(String pathFilename, double[] structure, HashMap<Integer,Integer> idToChr, String header, boolean... isTranslate) throws IOException{

		//number of fragments
		int n = structure.length / 3;
		
		//if idToChr is null, make the whole as one chromosome
		if (idToChr == null){
			idToChr = new HashMap<Integer,Integer>();
			for(int i = 0; i < n; i++){
				idToChr.put(i, 1);
			}
		}
		/////////
		
		if (isTranslate != null && isTranslate.length > 0 && isTranslate[0]) {
			//  make the minimum x-, y-, and z-coordinates of the structure equal to one (1)
			double translationX = Double.MAX_VALUE;
			double translationY = Double.MAX_VALUE;
			double translationZ = Double.MAX_VALUE;
			
			//  find the minimum x,y,z coordinate values
			for(int i = 0; i < n; i++){
				if(structure[i * 3] < translationX){
					translationX = structure[i];
				}
				if(structure[i * 3 + 1] < translationY){
					translationY = structure[i + 1];
				}
				if(structure[i * 3 + 2] < translationZ){
					translationZ = structure[i + 2];
				}
			}
			
			//  subtract one (1.0) to each of the translations to leave the minimum point at coordinate one (1.0) after the translation
			translationX -= 1;
			translationY -= 1;
			translationZ -= 1;
			
			//  translate all of the points in the structure
			for(int i = 0; i < n; i++){
				structure[i * 3] -= translationX;
				structure[i * 3 + 1] -= translationY;
				structure[i * 3 + 2] -= translationZ;
			}
			
		}
		

		//  write the structure to file
		try{						
			PrintWriter pw = new PrintWriter(pathFilename);
			
			
			//  write the header block
			pw.println(header.toUpperCase());

			int atomSerial = 1;
			int resName = 1;
			//String[] chain = {"O","A","B","C","D","E","F","G","H","I","J", "K", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y"};
							//"O","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24"

			String line;
			for(int i = 0; i < n; i++){
				line = "";
				line += getAtomString("ATOM");
				line += getSerialString(atomSerial + "") + " ";//12th space
				atomSerial++; // increase atom id
				
				line += getNameString("CA");
				line += getAltLocString(" ");
				//line += getResNameString(resName + "") + " ";//21st space
				line += getResNameString("MET" + "") + " ";//21st space
				
				
				//line += getChainIDString(chain[idToChr.get(i)]);
				line += getChainIDString((char)(idToChr.get(i) + 'A' ) + "");
				line += getResSeqString(resName + "");
				//line += getResSeqString("MET" + "");
				
				//if (atomSerial % 10 == 0){ // 10 atom in the same residue name
					resName++;
				//}
				
				line += " "; //27th space (iCode)
				line += "   "; //28,29,30th spaces
				
				line += getXString(structure[i * 3]);
				line += getYString(structure[i * 3 + 1]);
				line += getZString(structure[i * 3 + 2]);
				
				line += getOccupancyString(0.2);
				line += getTempFactorString(10.0);
				
				pw.println(line);				
				pw.flush();				
			}
			
			for(int i=1; i<atomSerial-1; i++){
				if(idToChr.get(i-1) == idToChr.get(i)){
					line = "";
				}else{
					line = "#";
				}
				
				line += getConnectString("CONECT");
				line += getSerial1String(i+"");
				line += getSerial2String((i+1) + "");
				pw.println(line);
				
				pw.flush();					
			}
			
			pw.println("END");			
			pw.flush();					
			pw.close();
		}catch (IOException e) {
			System.err.println("Error: could not output file: " + pathFilename);
			e.printStackTrace();
			throw e;
		}
	}
	
	/**
	 * Please refer to the pdb format file for detail of each column
	 * @param pathFilename: output file name
	 * @param structure: every 3 consecutive points is one point in 3D
	 * @param idToChr: to identify if 2 fragments belong to the same chromosome
	 * @param header for the pdb file
	 */
	public void writeStructureXYZ(String pathFilename, double[] structure, boolean... isTranslate) throws IOException{

		//number of fragments
		int n = structure.length / 3;
				
		
		if (isTranslate != null && isTranslate.length > 0 && isTranslate[0]) {
			//  make the minimum x-, y-, and z-coordinates of the structure equal to one (1)
			double translationX = Double.MAX_VALUE;
			double translationY = Double.MAX_VALUE;
			double translationZ = Double.MAX_VALUE;
			
			//  find the minimum x,y,z coordinate values
			for(int i = 0; i < n; i++){
				if(structure[i * 3] < translationX){
					translationX = structure[i];
				}
				if(structure[i * 3 + 1] < translationY){
					translationY = structure[i + 1];
				}
				if(structure[i * 3 + 2] < translationZ){
					translationZ = structure[i + 2];
				}
			}
			
			//  subtract one (1.0) to each of the translations to leave the minimum point at coordinate one (1.0) after the translation
			translationX -= 1;
			translationY -= 1;
			translationZ -= 1;
			
			//  translate all of the points in the structure
			for(int i = 0; i < n; i++){
				structure[i * 3] -= translationX;
				structure[i * 3 + 1] -= translationY;
				structure[i * 3 + 2] -= translationZ;
			}
			
		}
		

		//  write the structure to file
		try{						
			PrintWriter pw = new PrintWriter(pathFilename);
			
			
			//  write the header block
			pw.println(n);


			String line;
			for(int i = 0; i < n; i++){
				
				line = "C ";
				
				line += String.format("%10.5f",structure[i * 3]) + " ";
				line += String.format("%10.5f",structure[i * 3 + 1]) + " ";
				line += String.format("%10.5f",structure[i * 3 + 2]);
				
				pw.println(line);				
				pw.flush();				
			}
			
			pw.println("END");			
			pw.flush();					
			pw.close();
		}catch (IOException e) {
			System.err.println("Error: could not output file: " + pathFilename);
			e.printStackTrace();
			throw e;
		}
	}
	


	
	//to format the string to follow the format of pdb file format
	private String getAtomString(String st){
		//1-6
		int length = 6;
		int currentLength = st.length();
		if (currentLength > length){
			System.err.println("Error in atom name, length exceeds " +  length);
			return st.substring(0, length);
		}
		for(int i=0; i < length - currentLength; i++){
			st = st + " ";
		}
		
		return st;
	}
	
	//to format the string to follow the format of pdb file format
	private String getSerialString(String st){
		//7-11
		int length = 5;
		int currentLength = st.length();
		if (currentLength > length){
			System.err.println("Error in serial, length exceeds " + length);
			return st.substring(0, length);
		}
		for(int i=0; i < length - currentLength; i++){
			st = " " + st;
		}
		
		return st;
	}
	
	//to format the string to follow the format of pdb file format
	private String getNameString(String st){
		//13-16
		int length = 4;
		int currentLength = st.length();
		if (currentLength > length){
			System.err.println("Error in name, length exceeds " + length);
			return st.substring(0, length);
		}
		for(int i=0; i < length - currentLength; i++){
			st = " " +st;
		}
		
		return st;
	}

	//to format the string to follow the format of pdb file format
	private String getAltLocString(String st){
		//17
		int length = 1;
		int currentLength = st.length();
		if (currentLength > length){
			System.err.println("Error in alt loc, length exceeds " + length);
			return st.substring(0, length);
		}
		for(int i=0; i < length - currentLength; i++){
			st = st + " ";
		}
		
		return st;
	}
	
	//to format the string to follow the format of pdb file format
	private String getResNameString(String st){
		//18-20
		int length = 3;
		int currentLength = st.length();
		if (currentLength > length){
			System.err.println("Error in res name, length exceeds " + length);
			return st.substring(0, length);
		}
		for(int i=0; i < length - currentLength; i++){
			st = " " + st;
		}
		
		return st;
	}
	
	//to format the string to follow the format of pdb file format
	private String getChainIDString(String st){
		//22
		int length = 1;
		int currentLength = st.length();
		if (currentLength > length){
			System.err.println("Error in chain id, length exceeds " + length);
			return st.substring(0, length);
		}
		for(int i=0; i < length - currentLength; i++){
			st = st + " ";
		}
		
		return st;
	}
	
	//to format the string to follow the format of pdb file format
	private String getResSeqString(String st){
		//23-26
		int length = 4;
		int currentLength = st.length();
		if (currentLength > length){
			System.err.println("Error in res seq, length exceeds " + length);
			return st.substring(0, length);
		}
		for(int i=0; i < length - currentLength; i++){
			st = st + " ";
		}
		
		return st;
	}
	
	//to format the string to follow the format of pdb file format
	private String getXString(double x){
		//31-38
		int length = 8;
		String st = "";
		if (x > 10000 || x < -1000){
			st = String.format("%8.2f",x);
		}else{
			st = String.format("%8.3f",x);
		}
		int currentLength = st.length();
		if (currentLength > length){
			System.err.println("Error in X, length exceeds " + length);
			return st.substring(0, length);
		}
		for(int i=0; i < length - currentLength; i++){
			st = " "+st;
		}
		
		return st;
	}

	//to format the string to follow the format of pdb file format
	private String getYString(double x){
		//39-46
		int length = 8;
		String st = "";
		if (x > 10000 || x < -1000){
			st = String.format("%8.2f",x);
		}else{
			st = String.format("%8.3f",x);
		}
		int currentLength = st.length();
		if (currentLength > length){
			System.err.println("Error in Y, length exceeds " + length);
			return st.substring(0, length);
		}
		for(int i=0; i < length - currentLength; i++){
			st = " "+st;
		}
		
		return st;
	}
	//to format the string to follow the format of pdb file format
	private String getZString(double x){
		//47-54
		int length = 8;
		String st = "";
		if (x > 10000 || x < -1000){
			st = String.format("%8.2f",x);
		}else{
			st = String.format("%8.3f",x);
		}
		int currentLength = st.length();
		if (currentLength > length){
			System.err.println("Error in Z, length exceeds " + length);
			return st.substring(0, length);
		}
		for(int i=0; i < length - currentLength; i++){
			st = " "+st;
		}
		
		return st;
	}
	//to format the string to follow the format of pdb file format
	private String getOccupancyString(double x){
		//55-60
		int length = 6;
		String st = df2.format(x);
		int currentLength = st.length();
		if (currentLength > length){
			System.err.println("Error in occupancy, length exceeds " + length);
			return st.substring(0, length);
		}
		for(int i=0; i < length - currentLength; i++){
			st = " "+st;
		}
		
		return st;
	}
	//to format the string to follow the format of pdb file format
	private String getTempFactorString(double x){
		//61-66
		int length = 6;
		String st = df2.format(x);
		int currentLength = st.length();
		if (currentLength > length){
			System.err.println("Error in tempFactor, length exceeds " + length);
			return st.substring(0, length);
		}
		for(int i=0; i < length - currentLength; i++){
			st = " "+st;
		}		
		return st;
	}
	//to format the string to follow the format of pdb file format
	private String getConnectString(String st){
		//1-6
		int length = 6;
		int currentLength = st.length();
		if (currentLength > length){
			System.err.println("Error in connect, length exceeds " + length);
			return st.substring(0, length);
		}
		for(int i=0; i < length - currentLength; i++){
			st = st + " ";
		}
		
		return st;
	}
	//to format the string to follow the format of pdb file format
	private String getSerial1String(String st){
		//7-11
		int length = 5;
		int currentLength = st.length();
		if (currentLength > length){
			System.err.println("Error in serial 1, length exceeds " + length);
			return st.substring(0, length);
		}
		for(int i=0; i < length - currentLength; i++){
			st = " " + st;
		}		
		return st;
	}
	//to format the string to follow the format of pdb file format	
	private String getSerial2String(String st){
		//12-16
		int length = 5;
		int currentLength = st.length();
		if (currentLength > length){
			System.err.println("Error in serial 2, length exceeds " + length);
			return st.substring(0, length);
		}
		for(int i=0; i < length - currentLength; i++){
			st = " " + st;
		}		
		return st;
	}


}
