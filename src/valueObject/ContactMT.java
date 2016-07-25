package valueObject;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * This object contains the contact matrix, if the dataset is small, matrix is used, if not list is used
 * @author Tuan
 *
 */
public class ContactMT {

	private static ContactMT contactMT;
	
	//number of points or number of rows/columns
	private int size;
	
	//length of the array containing contacts
	private int len;
	//this 1-D array will be converted to 2D when looking up
	//this matrix is equal to a lower diagonal 2-D matrix
	private double[] mt = null;
	
	//if size is too big such that heap cannot accomodate mt, lst will be used
	private List<Contact> lst = null;
	
	//there is only one contact matrix
	public static ContactMT getContactMTInstance(int n){
		if (contactMT == null){
			synchronized(ContactMT.class){
				if (contactMT == null){
					contactMT = new ContactMT(n);
				}
			}
		}
		return contactMT;
	}

	
	private ContactMT(int n){
		size = n;
		
		//about 10GB heap 
		if (size > 12000){
			lst = new ArrayList<Contact>();
		}else{
			len = ((n - 1) * n) / 2;
			mt = new double[len];
		}
	}
	
	public double countPercentBelowThres(double thres){
		long count = 0;
		if (mt != null){		
			
			int pos;
			for(int i = 0; i < size; i++){
				for(int j = i + 1; j < size; j++){
					pos = getPos(i,j);
					if (mt[pos] < thres){
						count++;
					}
				}
			}
		}else{			
			for(Contact ct : lst){
				if (ct.getIF() < thres){
					count++;
				}
			}
		}		
		return 100.0 * count * 2 / (size - 1) * size;
	}
	
	public double getTotalIF(){
		double total = 0;
		if (mt != null){		
			
			int pos;
			for(int i = 0; i < size; i++){
				for(int j = i + 1; j < size; j++){
					pos = getPos(i,j);
					total += mt[pos];					
				}
			}
		}else{			
			for(Contact ct : lst){
				total += ct.getIF();
			}
		}		
		return total;
	}
	
	public long getNbrContact(){
		long count = 0;
		if (mt != null){		
			
			int pos;
			for(int i = 0; i < size; i++){
				for(int j = i + 1; j < size; j++){
					pos = getPos(i,j);
					if (mt[pos] > 0){
						count++;
					}
				}
			}
			return count;
			
		}else{			
			return lst.size();
		}		
		
	}
	
	
	public void addContact(int pos1, int pos2, double IF){
		//make sure pos1 > pos2, because mt is lower diagonal matrix
		if (pos2 > pos1){
			int tg = pos1;
			pos1 = pos2;
			pos2 = tg;
		}else if(pos1 == pos2){
			return;
		}

		if (mt != null){		
			int pos = getPos(pos1,pos2);
			mt[pos] = IF;
			
		}else{
			//remember to sort this list after adding all contacts
			lst.add(new Contact(pos1,pos2,IF));
		}
	}

	
	public double getContact(int pos1, int pos2){
		if (pos1 == pos2){
			return 0;
		}
		//make sure pos1 > pos2, because mt is lower diagonal matrix
		if (pos2 > pos1){
			int tg = pos1;
			pos1 = pos2;
			pos2 = tg;
		}
		
		if (mt != null){
			
			int pos = getPos(pos1,pos2);
			return mt[pos];
		}else{
			int id = Collections.binarySearch(lst, new Contact(pos1,pos2,0.0));
			if (id < 0){
				return 0.0;
			}else{
//				System.out.println(id);
//				System.out.println(lst.get(0).getPos1() + " " + lst.get(0).getPos2());
//				System.out.println(lst.get(1).getPos1() + " " + lst.get(1).getPos2());
//				System.out.println(lst.get(2).getPos1() + " " + lst.get(2).getPos2());
				return lst.get(id).getIF();
			}
		}
	}
	
	public boolean updateContact(int pos1, int pos2, double IF){
		if (pos1 == pos2){
			return true;
		}
		//make sure pos1 > pos2, because mt is lower diagonal matrix
		if (pos2 > pos1){
			int tg = pos1;
			pos1 = pos2;
			pos2 = tg;
		}		
		
		if (mt != null){
			int pos = getPos(pos1,pos2);
			mt[pos] = IF;
			return true;
		}else{
			int id = Collections.binarySearch(lst, new Contact(pos1,pos2,0.0));
			if (id < 0){
				return false;
			}else{
				lst.get(id).setIF(IF);
				return true;
			}
		}
	}
	
	
	//sort contacts after add all contacts so that binary search can be used later
	public void sortContact(){
		if (lst != null){
			Collections.sort(lst);			
		}
	}

	public int getSize(){
		return size;
	}
	
	//Calculate position of [pos1,pos2] in mt
	private int getPos(int pos1, int pos2){
		if (pos2 > pos1){
			int tg = pos1;
			pos1 = pos2;
			pos2 = tg;
		}		
		return (pos1 * (pos1 - 1)) / 2 + pos2;
	}
	
}
