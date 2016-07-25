package valueObject;
/**
 * this value object contains scores and other information obtained when evaluating the structure
 * @author Tuan
 *
 */
public class StructureEvaluationVO {
	//percentage of contacts in the input that are satisfied
	private double[][] contactScore;
	//percentage of non-contact in the input that are satisfied
	private double[][] nonContactScore;
	//percentage of interaction frequency that is satisfied
	private double[][] satisfiedIFPercent;
	
	private double totalContactScore;
	private double totalNonContactScore;
	
	//percentage of contact
	private double[][] contactPercent;

	public StructureEvaluationVO(double[][] contactScore,
			double[][] nonContactScore, double[][] satisfiedIFPercent,
			double totalContactScore, double totalNonContactScore,
			double[][] contactPercent) {
		
		super();
		this.contactScore = contactScore;
		this.nonContactScore = nonContactScore;
		this.satisfiedIFPercent = satisfiedIFPercent;
		this.totalContactScore = totalContactScore;
		this.totalNonContactScore = totalNonContactScore;
		this.contactPercent = contactPercent;
		//this.cor = cor;
	}

	public double[][] getContactScore() {
		return contactScore;
	}

	public void setContactScore(double[][] contactScore) {
		this.contactScore = contactScore;
	}

	public double[][] getNonContactScore() {
		return nonContactScore;
	}

	public void setNonContactScore(double[][] nonContactScore) {
		this.nonContactScore = nonContactScore;
	}

	public double[][] getSatisfiedIFPercent() {
		return satisfiedIFPercent;
	}

	public void setSatisfiedIFPercent(double[][] satisfiedIFPercent) {
		this.satisfiedIFPercent = satisfiedIFPercent;
	}

	public double[][] getContactPercent() {
		return contactPercent;
	}

	public void setContactPercent(double[][] contactPercent) {
		this.contactPercent = contactPercent;
	}

	public double getTotalContactScore() {
		return totalContactScore;
	}

	public void setTotalContactScore(double totalContactScore) {
		this.totalContactScore = totalContactScore;
	}

	public double getTotalNonContactScore() {
		return totalNonContactScore;
	}

	public void setTotalNonContactScore(double totalNonContactScore) {
		this.totalNonContactScore = totalNonContactScore;
	}
	
	
	
	
}