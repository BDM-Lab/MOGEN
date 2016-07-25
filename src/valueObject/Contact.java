package valueObject;

public class Contact implements Comparable<Contact> {
	private int pos1;
	private int pos2;
	private double IF;
	private int chr1;
	private int chr2;
		
	public Contact(int i, int ch1, int j, int ch2, double iF) {
		super();
		//make sure that j > i, which is upper triangle
		//(1,2) (1,3) (1,4) (1,5)....
		if (i > j){
			this.pos2 = i;
			this.chr2 = ch1;
			
			this.pos1 = j;
			this.chr1 = ch2;
			
		}else{		
			this.pos1 = i;
			this.chr1 = ch1;
			
			this.pos2 = j;
			this.chr2 = ch2;			
		}
		
		
		IF = iF;
	}
	
	public Contact(int i, int j, double iF) {
		super();
		//make sure that j > i, which is upper triangle
		//(1,2) (1,3) (1,4) (1,5)....
		if (i > j){
			this.pos2 = i;						
			this.pos1 = j;		
			
		}else{		
			this.pos1 = i;						
			this.pos2 = j;			
		}
				
		IF = iF;
	}

	
	public int getPos1() {
		return pos1;
	}
	public void setPos1(int pos1) {
		this.pos1 = pos1;
	}
	public int getPos2() {
		return pos2;
	}
	public void setPos2(int pos2) {
		this.pos2 = pos2;
	}
	public double getIF() {
		return IF;
	}
	public void setIF(double iF) {
		IF = iF;
	}
		
	public int getChr1() {
		return this.chr1;
	}
	public void setChr1(int chr1) {
		this.chr1 = chr1;
	}
	public int getChr2() {
		return this.chr2;
	}
	public void setChr2(int chr2) {
		this.chr2 = chr2;
	}
	
	@Override
	public int compareTo(Contact ct) {		
		if ((pos1 < ct.getPos1()) 
				|| (pos1 == ct.getPos1() && pos2 < ct.getPos2()) ){
		
			return -1;
		
		}else if (pos1 == ct.getPos1() && pos2 == ct.getPos2()){
			return 0;
		
		}else{
			return 1;
		}		
	}
	
	@Override
	public boolean equals(Object o){
		if (! (o instanceof Contact)){
			return false;
		}
		
		Contact ct = (Contact) o;
		if (ct.getPos1() == this.pos1 && ct.getPos2() == this.pos2){			
			return true;
		}
		
		return false;
	}
	
}
