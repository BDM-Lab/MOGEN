package evaluation;

import java.util.HashMap;
import java.util.List;

import utility.Helper;
import valueObject.Contact;
import valueObject.ContactMT;
import valueObject.StructureEvaluationVO;

/**
 * This class is to calculate scores for structures
 * @author Tuan
 *
 */
public class StructureScoring {
		
	public StructureEvaluationVO calculateScore(double[] str,ContactMT contactMT, double contactDist, 
			HashMap<Integer,Integer> idToChr, int num_of_chr){
		
		StructureEvaluationVO vo=null;
		Helper helper = Helper.getHelperInstance();
		
		int chr1,chr2;
		double d;
		//the number of fragments
		int n = contactMT.getSize();
		
		//
		if (idToChr == null){
			idToChr = new HashMap<Integer,Integer>();
			for(int i = 0; i < n; i++){
				idToChr.put(i, 1);
			}
			num_of_chr = 1;
		}
		
		
		//possible number of contacts between chromosomes
		int[][] possibleNbrContact = new int[num_of_chr + 1][num_of_chr + 1];
		
		//the number of contacts in the input data
		int[][] contactNbr = new int[num_of_chr + 1][num_of_chr + 1];
		
		//the number of non-contact in the input data
		int[][] nonContactNbr = new int[num_of_chr + 1][num_of_chr + 1];
		
		//the number of contacts that are satisfied
		int[][] satContactNbr = new int[num_of_chr + 1][num_of_chr + 1];
		
		//the number of non-contacts that are satisfied
		int[][] satNonContactNbr = new int[num_of_chr + 1][num_of_chr + 1];
		
		//total number of contacts
		int totalContacts = 0;
		//total number of non-contact
		int totalNonContacts = 0;
		
		int totalSatContacts = 0;
		int totalSatNonContacts = 0;
		
		
		double[][] contactScore = new double[num_of_chr + 1][num_of_chr + 1];
		double[][] nonContactScore = new double[num_of_chr + 1][num_of_chr + 1];
		double[][] satisfiedIFPercent = new double[num_of_chr + 1][num_of_chr + 1];
		double[][] contactPercent = new double[num_of_chr + 1][num_of_chr + 1];
		double totalContactScore;
		double totalNonContactScore;
				
		//sum of interaction frequency (IF) satisfied
		double[][] satIF = new double[num_of_chr + 1][num_of_chr + 1];
		double[][] totalIF = new double[num_of_chr + 1][num_of_chr + 1];
		
		//array contains distances and IFs to calculate correlation
		//double[] dists = new double[n * (n - 1) / 2];
		//double[] counts = new double[n * (n - 1) / 2];
		
		
		double f;
		for(int i = 0; i < n; i++){
			for(int j = i + 1; j < n; j++){
				chr1 = idToChr.get(i);
				chr2 = idToChr.get(j);	
				
				d = helper.calEuclidianDist(str[i * 3], str[i * 3 + 1], str[i * 3 + 2], str[j * 3], str[j * 3 + 1], str[j * 3 + 2]);
				
				f = contactMT.getContact(i,j);
				//to calculate correlation
				//dists[co] = d;
				//counts[co] = f;				
				//

				if (chr1 != chr2){
					totalIF[chr1][chr2] += f;
					totalIF[chr2][chr1] += f;
					
					possibleNbrContact[chr1][chr2]++;
					possibleNbrContact[chr2][chr1]++;
				}else{
					totalIF[chr1][chr1] += f;
					possibleNbrContact[chr1][chr1]++;
				}
				
				//if it's a contact
				if (f > 0.0 ){
					
					totalContacts++;
					
					if (chr1 != chr2){
						contactNbr[chr1][chr2]++;
						contactNbr[chr2][chr1]++;
					}else{
						contactNbr[chr1][chr1]++;
					}
					
					//if the contact is satisfied
					if (d <= contactDist){
						
						totalSatContacts++;
						
						if (chr1 != chr2){
							satContactNbr[chr1][chr2]++;
							satContactNbr[chr2][chr1]++;
							
							satIF[chr1][chr2] += f;
							satIF[chr2][chr1] += f;
						}else{
							satContactNbr[chr1][chr1]++;
							satIF[chr1][chr1] += f;
						}
					}					
					
				}else{ //a non-contact
					totalNonContacts++;
					
					if (chr1 != chr2){
						nonContactNbr[chr1][chr2]++;					
						nonContactNbr[chr2][chr1]++;
					}else{
						nonContactNbr[chr1][chr1]++;
					}
					//if this non-contact is satisfied
					if (d > contactDist){
						totalSatNonContacts++;
						
						if (chr1 != chr2){
							satNonContactNbr[chr1][chr2]++;
							satNonContactNbr[chr2][chr1]++;
						}else{
							satNonContactNbr[chr1][chr1]++;
						}
					}

				}
			}
		}
		
		for(int i = 1; i <= num_of_chr; i++){
			for(int j = 1; j <= num_of_chr; j++){
				contactScore[i][j] = satContactNbr[i][j] * 100.0 / contactNbr[i][j];
				nonContactScore[i][j] = satNonContactNbr[i][j] * 100.0 / nonContactNbr[i][j];
				satisfiedIFPercent[i][j] = satIF[i][j] * 100.0 / totalIF[i][j];
				contactPercent[i][j] = contactNbr[i][j] * 100.0 / possibleNbrContact[i][j];
			}
		}
		
		totalContactScore = (double) totalSatContacts * 100.0 / (double)totalContacts;
		totalNonContactScore = (double) totalSatNonContacts * 100.0 / (double) totalNonContacts;

		vo = new StructureEvaluationVO(contactScore, nonContactScore, satisfiedIFPercent, totalContactScore, totalNonContactScore, contactPercent);
		
		/*
		center(dists);
		center(counts); 
		double cor = product(dists, counts) / (norm2(dists) * norm2(counts));
		vo.setCor(cor);
		*/
		
		return vo ;
	}
	
	public StructureEvaluationVO calculateScore(double[] str,List<Contact> ctLst, double contactDist, 
			HashMap<Integer,Integer> idToChr, int num_of_chr){
		
		StructureEvaluationVO vo=null;
		Helper helper = Helper.getHelperInstance();
		
		int chr1,chr2,i,j;
		double d;
		//the number of fragments
		int n = str.length / 3;
		
		//
		if (idToChr == null){
			idToChr = new HashMap<Integer,Integer>();
			for(i = 0; i < n; i++){
				idToChr.put(i, 1);
			}
			num_of_chr = 1;
		}
		
		
		//possible number of contacts between chromosomes
		int[][] possibleNbrContact = new int[num_of_chr + 1][num_of_chr + 1];
		
		//the number of contacts in the input data
		int[][] contactNbr = new int[num_of_chr + 1][num_of_chr + 1];
		
		//the number of non-contact in the input data
		int[][] nonContactNbr = new int[num_of_chr + 1][num_of_chr + 1];
		
		//the number of contacts that are satisfied
		int[][] satContactNbr = new int[num_of_chr + 1][num_of_chr + 1];
		
		//the number of non-contacts that are satisfied
		int[][] satNonContactNbr = new int[num_of_chr + 1][num_of_chr + 1];
		
		//total number of contacts
		int totalContacts = 0;
		//total number of non-contact
		int totalNonContacts = 0;
		
		int totalSatContacts = 0;
		int totalSatNonContacts = 0;
		
		
		double[][] contactScore = new double[num_of_chr + 1][num_of_chr + 1];
		double[][] nonContactScore = new double[num_of_chr + 1][num_of_chr + 1];
		double[][] satisfiedIFPercent = new double[num_of_chr + 1][num_of_chr + 1];
		double[][] contactPercent = new double[num_of_chr + 1][num_of_chr + 1];
		double totalContactScore;
		double totalNonContactScore;
				
		//sum of interaction frequency (IF) satisfied
		double[][] satIF = new double[num_of_chr + 1][num_of_chr + 1];
		double[][] totalIF = new double[num_of_chr + 1][num_of_chr + 1];
		
		//array contains distances and IFs to calculate correlation
		//double[] dists = new double[n * (n - 1) / 2];
		//double[] counts = new double[n * (n - 1) / 2];
		
		
		double f;
		for(Contact ct: ctLst){
			i = ct.getPos1();
			chr1 = ct.getChr1();
			j = ct.getPos2();
			chr2 = ct.getChr2();	
			f = ct.getIF();
			
			d = helper.calEuclidianDist(str[i * 3], str[i * 3 + 1], str[i * 3 + 2], str[j * 3], str[j * 3 + 1], str[j * 3 + 2]);
		

			if (chr1 != chr2){
				totalIF[chr1][chr2] += f;
				totalIF[chr2][chr1] += f;
				
				possibleNbrContact[chr1][chr2]++;
				possibleNbrContact[chr2][chr1]++;
			}else{
				totalIF[chr1][chr1] += f;
				possibleNbrContact[chr1][chr1]++;
			}
			
			//if it's a contact
			if (f > 0.0 ){
				
				totalContacts++;
				
				if (chr1 != chr2){
					contactNbr[chr1][chr2]++;
					contactNbr[chr2][chr1]++;
				}else{
					contactNbr[chr1][chr1]++;
				}
				
				//if the contact is satisfied
				if (d <= contactDist){
					
					totalSatContacts++;
					
					if (chr1 != chr2){
						satContactNbr[chr1][chr2]++;
						satContactNbr[chr2][chr1]++;
						
						satIF[chr1][chr2] += f;
						satIF[chr2][chr1] += f;
					}else{
						satContactNbr[chr1][chr1]++;
						satIF[chr1][chr1] += f;
					}
				}					
				
			}else{ //a non-contact
				totalNonContacts++;
				
				if (chr1 != chr2){
					nonContactNbr[chr1][chr2]++;					
					nonContactNbr[chr2][chr1]++;
				}else{
					nonContactNbr[chr1][chr1]++;
				}
				//if this non-contact is satisfied
				if (d > contactDist){
					totalSatNonContacts++;
					
					if (chr1 != chr2){
						satNonContactNbr[chr1][chr2]++;
						satNonContactNbr[chr2][chr1]++;
					}else{
						satNonContactNbr[chr1][chr1]++;
					}
				}

			}
			
		}
		
		for(i = 1; i <= num_of_chr; i++){
			for(j = 1; j <= num_of_chr; j++){
				contactScore[i][j] = satContactNbr[i][j] * 100.0 / contactNbr[i][j];
				nonContactScore[i][j] = satNonContactNbr[i][j] * 100.0 / nonContactNbr[i][j];
				satisfiedIFPercent[i][j] = satIF[i][j] * 100.0 / totalIF[i][j];
				contactPercent[i][j] = contactNbr[i][j] * 100.0 / possibleNbrContact[i][j];
			}
		}
		
		totalContactScore = (double) totalSatContacts * 100.0 / (double)totalContacts;
		totalNonContactScore = (double) totalSatNonContacts * 100.0 / (double) totalNonContacts;

		vo = new StructureEvaluationVO(contactScore, nonContactScore, satisfiedIFPercent, totalContactScore, totalNonContactScore, contactPercent);
		return vo ;
	}

	public StructureEvaluationVO calculateScore(double[] str,HashMap<Integer,Contact> ctMap, double contactDist, 
			HashMap<Integer,Integer> idToChr, int num_of_chr){
		
		StructureEvaluationVO vo=null;
		Helper helper = Helper.getHelperInstance();
		
		int chr1,chr2,i,j,val;
		double d;
		//the number of fragments
		int n = str.length / 3;
		
		//
		if (idToChr == null){
			idToChr = new HashMap<Integer,Integer>();
			for(i = 0; i < n; i++){
				idToChr.put(i, 1);
			}
			num_of_chr = 1;
		}
		
		
		//possible number of contacts between chromosomes
		int[][] possibleNbrContact = new int[num_of_chr + 1][num_of_chr + 1];
		
		//the number of contacts in the input data
		int[][] contactNbr = new int[num_of_chr + 1][num_of_chr + 1];
		
		//the number of non-contact in the input data
		int[][] nonContactNbr = new int[num_of_chr + 1][num_of_chr + 1];
		
		//the number of contacts that are satisfied
		int[][] satContactNbr = new int[num_of_chr + 1][num_of_chr + 1];
		
		//the number of non-contacts that are satisfied
		int[][] satNonContactNbr = new int[num_of_chr + 1][num_of_chr + 1];
		
		//total number of contacts
		int totalContacts = 0;
		//total number of non-contact
		int totalNonContacts = 0;
		
		int totalSatContacts = 0;
		int totalSatNonContacts = 0;
		
		
		double[][] contactScore = new double[num_of_chr + 1][num_of_chr + 1];
		double[][] nonContactScore = new double[num_of_chr + 1][num_of_chr + 1];
		double[][] satisfiedIFPercent = new double[num_of_chr + 1][num_of_chr + 1];
		double[][] contactPercent = new double[num_of_chr + 1][num_of_chr + 1];
		double totalContactScore;
		double totalNonContactScore;
				
		//sum of interaction frequency (IF) satisfied
		double[][] satIF = new double[num_of_chr + 1][num_of_chr + 1];
		double[][] totalIF = new double[num_of_chr + 1][num_of_chr + 1];
		
		//array contains distances and IFs to calculate correlation
		//double[] dists = new double[n * (n - 1) / 2];
		//double[] counts = new double[n * (n - 1) / 2];
		
		
		double f;
		for(i = 0; i <= n; i ++){
			for(j = i + 1; j < n; j++){
				
			
				val = i * n + j;
				chr1 = idToChr.get(i);
				chr2 = idToChr.get(j);			
				
				d = helper.calEuclidianDist(str[i * 3], str[i * 3 + 1], str[i * 3 + 2], str[j * 3], str[j * 3 + 1], str[j * 3 + 2]);
		
				f = 0.0;
				if (ctMap.containsKey(val) ){
					f = ctMap.get(val).getIF();
				}

				if (chr1 != chr2){
					totalIF[chr1][chr2] += f;
					totalIF[chr2][chr1] += f;
					
					possibleNbrContact[chr1][chr2]++;
					possibleNbrContact[chr2][chr1]++;
				}else{
					totalIF[chr1][chr1] += f;
					possibleNbrContact[chr1][chr1]++;
				}
			
				//if it's a contact
				if (f > 0.0 ){
					
					totalContacts++;
					
					if (chr1 != chr2){
						contactNbr[chr1][chr2]++;
						contactNbr[chr2][chr1]++;
					}else{
						contactNbr[chr1][chr1]++;
					}
					
					//if the contact is satisfied
					if (d <= contactDist){
						
						totalSatContacts++;
						
						if (chr1 != chr2){
							satContactNbr[chr1][chr2]++;
							satContactNbr[chr2][chr1]++;
							
							satIF[chr1][chr2] += f;
							satIF[chr2][chr1] += f;
						}else{
							satContactNbr[chr1][chr1]++;
							satIF[chr1][chr1] += f;
						}
					}					
					
				}else{ //a non-contact
					totalNonContacts++;
					
					if (chr1 != chr2){
						nonContactNbr[chr1][chr2]++;					
						nonContactNbr[chr2][chr1]++;
					}else{
						nonContactNbr[chr1][chr1]++;
					}
					//if this non-contact is satisfied
					if (d > contactDist){
						totalSatNonContacts++;
						
						if (chr1 != chr2){
							satNonContactNbr[chr1][chr2]++;
							satNonContactNbr[chr2][chr1]++;
						}else{
							satNonContactNbr[chr1][chr1]++;
						}
					}
	
				}
			}
			
		}
		
		for(i = 1; i <= num_of_chr; i++){
			for(j = 1; j <= num_of_chr; j++){
				contactScore[i][j] = satContactNbr[i][j] * 100.0 / contactNbr[i][j];
				nonContactScore[i][j] = satNonContactNbr[i][j] * 100.0 / nonContactNbr[i][j];
				satisfiedIFPercent[i][j] = satIF[i][j] * 100.0 / totalIF[i][j];
				contactPercent[i][j] = contactNbr[i][j] * 100.0 / possibleNbrContact[i][j];
			}
		}
		
		totalContactScore = (double) totalSatContacts * 100.0 / (double)totalContacts;
		totalNonContactScore = (double) totalSatNonContacts * 100.0 / (double) totalNonContacts;

		vo = new StructureEvaluationVO(contactScore, nonContactScore, satisfiedIFPercent, totalContactScore, totalNonContactScore, contactPercent);
		return vo;
		
		
	}
	
	
	
	public StructureEvaluationVO calculateScore(double[] str,double[][] contactMT, double contactDist, 
			HashMap<Integer,Integer> idToChr, int num_of_chr){
		
		StructureEvaluationVO vo=null;
		Helper helper = Helper.getHelperInstance();
		
		int chr1,chr2;
		double d;
		//the number of fragments
		int n = contactMT.length;
		
		//
		if (idToChr == null){
			idToChr = new HashMap<Integer,Integer>();
			for(int i = 0; i < n; i++){
				idToChr.put(i, 1);
			}
			num_of_chr = 1;
		}
		
		
		//possible number of contacts between chromosomes
		int[][] possibleNbrContact = new int[num_of_chr + 1][num_of_chr + 1];
		
		//the number of contacts in the input data
		int[][] contactNbr = new int[num_of_chr + 1][num_of_chr + 1];
		
		//the number of non-contact in the input data
		int[][] nonContactNbr = new int[num_of_chr + 1][num_of_chr + 1];
		
		//the number of contacts that are satisfied
		int[][] satContactNbr = new int[num_of_chr + 1][num_of_chr + 1];
		
		//the number of non-contacts that are satisfied
		int[][] satNonContactNbr = new int[num_of_chr + 1][num_of_chr + 1];
		
		//total number of contacts
		int totalContacts = 0;
		//total number of non-contact
		int totalNonContacts = 0;
		
		int totalSatContacts = 0;
		int totalSatNonContacts = 0;
		
		
		double[][] contactScore = new double[num_of_chr + 1][num_of_chr + 1];
		double[][] nonContactScore = new double[num_of_chr + 1][num_of_chr + 1];
		double[][] satisfiedIFPercent = new double[num_of_chr + 1][num_of_chr + 1];
		double[][] contactPercent = new double[num_of_chr + 1][num_of_chr + 1];
		double totalContactScore;
		double totalNonContactScore;
				
		//sum of interaction frequency (IF) satisfied
		double[][] satIF = new double[num_of_chr + 1][num_of_chr + 1];
		double[][] totalIF = new double[num_of_chr + 1][num_of_chr + 1];
		
		//array contains distances and IFs to calculate correlation
		//double[] dists = new double[n * (n - 1) / 2];
		//double[] counts = new double[n * (n - 1) / 2];
		
		//int co = 0;
		double f;
		for(int i = 0; i < n; i++){
			for(int j = i + 1; j < n; j++){
				chr1 = idToChr.get(i);
				chr2 = idToChr.get(j);	
				
				d = helper.calEuclidianDist(str[i * 3], str[i * 3 + 1], str[i * 3 + 2], str[j * 3], str[j * 3 + 1], str[j * 3 + 2]);
				
				f = contactMT[i][j];
				//to calculate correlation
				//dists[co] = d;
				//counts[co] = f;
				//co++;
				//

				if (chr1 != chr2){
					totalIF[chr1][chr2] += f;
					totalIF[chr2][chr1] += f;
					
					possibleNbrContact[chr1][chr2]++;
					possibleNbrContact[chr2][chr1]++;
				}else{
					totalIF[chr1][chr1] += f;
					possibleNbrContact[chr1][chr1]++;
				}
				
				//if it's a contact
				if (f > 0.0 ){
					
					totalContacts++;
					
					if (chr1 != chr2){
						contactNbr[chr1][chr2]++;
						contactNbr[chr2][chr1]++;
					}else{
						contactNbr[chr1][chr1]++;
					}
					
					//if the contact is satisfied
					if (d <= contactDist){
						
						totalSatContacts++;
						
						if (chr1 != chr2){
							satContactNbr[chr1][chr2]++;
							satContactNbr[chr2][chr1]++;
							
							satIF[chr1][chr2] += f;
							satIF[chr2][chr1] += f;
						}else{
							satContactNbr[chr1][chr1]++;
							satIF[chr1][chr1] += f;
						}
					}					
					
				}else{ //a non-contact
					totalNonContacts++;
					
					if (chr1 != chr2){
						nonContactNbr[chr1][chr2]++;					
						nonContactNbr[chr2][chr1]++;
					}else{
						nonContactNbr[chr1][chr1]++;
					}
					//if this non-contact is satisfied
					if (d > contactDist){
						totalSatNonContacts++;
						
						if (chr1 != chr2){
							satNonContactNbr[chr1][chr2]++;
							satNonContactNbr[chr2][chr1]++;
						}else{
							satNonContactNbr[chr1][chr1]++;
						}
					}

				}
			}
		}
		
		for(int i = 1; i <= num_of_chr; i++){
			for(int j = 1; j <= num_of_chr; j++){
				contactScore[i][j] = satContactNbr[i][j] * 100.0 / contactNbr[i][j];
				nonContactScore[i][j] = satNonContactNbr[i][j] * 100.0 / nonContactNbr[i][j];
				satisfiedIFPercent[i][j] = satIF[i][j] * 100.0 / totalIF[i][j];
				contactPercent[i][j] = contactNbr[i][j] * 100.0 / possibleNbrContact[i][j];
			}
		}
		
		totalContactScore = (double) totalSatContacts * 100.0 / (double)totalContacts;
		totalNonContactScore = (double) totalSatNonContacts * 100.0 / (double) totalNonContacts;

		vo = new StructureEvaluationVO(contactScore, nonContactScore, satisfiedIFPercent, totalContactScore, totalNonContactScore, contactPercent);
		
		/*
		center(dists);
		center(counts); 
		double cor = product(dists, counts) / (norm2(dists) * norm2(counts));
		vo.setCor(cor);
		*/
		
		return vo ;
	}

}