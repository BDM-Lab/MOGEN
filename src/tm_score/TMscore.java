package tm_score;
//*************************************************************************
//*     This program is to compare two protein structures and identify the 
//*     best superposition that has the highest TM-score. Input structures 
//*     must be in the PDB format. By default, TM-score is normalized by 
//*     the second protein. Users can obtain a brief instruction by simply
//*     running the program without arguments. For comments/suggestions,
//*     please contact email: zhng@umich.edu.fun
//*     
//*     Reference: 
//*     Yang Zhang, Jeffrey Skolnick, Proteins, 2004 57:702-10.
//*     
//*     Permission to use, copy, modify, and distribute this program for 
//*     any purpose, with or without fee, is hereby granted, provided that
//*     the notices on the head, the reference rmsdinformation, and this
//*     copyright notice appear in all copies or substantial portions of 
//*     the Software. It is provided "as is" without express or implied 
//*     warranty.
//******************* Updating history ************************************
//*     2005/10/19: the program was reformed so that the score values.
//*                 are not dependent on the specific compilers.
//*     2006/06/20: selected 'A' if there is altLoc when reading PDB file.
//*     2007/02/05: fixed a bug with length<15 in TMscore_32.
//*     2007/02/27: rotation matrix from Chain-1 to Chain-2 was added.
//*     2007/12/06: GDT-HA score was added, fixed a bug for reading PDB.
//*     2010/08/02: A new RMSD matrix was used and obsolete statement removed.
//*     2011/01/03: The length of pdb file names were extended to 500.
//*     2011/01/30: An open source license is attached to the program.
//*************************************************************************
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;

public class TMscore {
	 private int nmax = 3000;
	 private double d, d0, d0_fix;
	 private int nseqA, nseqB;
	 private int m_out = -1;
	 private int m_fix = -1;
	 private int ier = 0;
	 private int n_ali;
	 private int n_cut; // ![1,n_ali],align residues for the score
	 private double score, score_maxsub, score_fix, score10;
	 private double n_GDT05, n_GDT1, n_GDT2, n_GDT4, n_GDT8;
	 private double score_max, score_fix_max;
	 //tuan add on Dec 28, to use GDT-HA score instead of TM-score
	 double score_GDT_HA;
	 private double gdt8 = 4.0, gdt4 = 2.0, gdt2 = 1.5, gdt1 = 1.0, gdt05 = 0.5;	 //
	 private double rms, drms, armsd, rmsd_ali,rmsd;
	 private String fname = "";
	 private String pdb = "";
	 private String outname = "";
	 private double d0_search;
	 private double d_output;
	 private int[] L_ini = new int[nmax];
	 private int[] iA = new int[nmax];
	 private int[] iB = new int[nmax];
	 private int[] i_ali = new int[nmax];
	 private int[] nresA = new int[nmax];
	 private int[] nresB = new int[nmax];	
	 private double[][] u = new double[4][4];
	 private double[] t = new double[4];
	 private double[] iq = new double[nmax];
	 private double[] xa = new double[nmax];
	 private double[] ya = new double[nmax];
	 private double[] za = new double[nmax];
	 private double[] xt = new double[nmax];
	 private double[] yt = new double[nmax];
	 private double[] zt = new double[nmax];
	 private double[] xb = new double[nmax];
	 private double[] yb = new double[nmax];
	 private double[] zb = new double[nmax];		
	 private double[][] r_1 = new double[3 + 1][nmax];
	 private double[][] r_2 = new double[3 + 1][nmax];

	public static void main(String[] args) throws Exception {
		
		TMscore tmscore = new TMscore();
		tmscore.runTMscore(args);
	}
	
	//Tuan added
	public TMscore(){
		
	}
	public TMscore(double gdt8, double gdt4, double gdt2, double gdt1, double gdt05){
		this.gdt8 = gdt8;
		this.gdt4 = gdt4;
		this.gdt2 = gdt2;
		this.gdt1 = gdt1;
		this.gdt05 = gdt05;
	}
	//end
	
	public void runTMscore(String[] args) throws Exception {
		int i, j;
		int k = 0;
		d0_fix = -1;				
		int[] k_ali = new int[nmax];
		int[] k_ali0 = new int[nmax];
		String seq1A = "";		
		j = 0;
		if(args.length<2){
			System.out.println(help());
			System.exit(1);
		}
		String[] pdbs = new String[2];		
		for(i=0;i<args.length;){
			if(args[i].equals("-o")){
				m_out = 1;
				i++;
				outname = args[i];
				i++;
			}else if(args[i].equals("-d")){
				m_fix = 1;
				i++;
				d0_fix = Double.valueOf(args[i]);
				i++;
			}else if(j<2){
				pdbs[j] = args[i];
				i++;
				j++;
			}
		}
		
		seq1A = readPDB(pdbs[0], xa, ya, za, nresA);
		String seq1B = "";
		seq1B = readPDB(pdbs[1], xb, yb, zb, nresB);
		nseqA = seq1A.length();
		nseqB = seq1B.length();				
		nseqA--;
		nseqB--;
		int ka0 = 0;			
		// pickup the aligned residues:
		
		for (i = 1; i<=nseqA; i++) {
	OUT1:	for (j = 1; j <=nseqB; j++) {
				if (nresA[i] == nresB[j]) {
					k++;
					iA[k] = i;
					iB[k] = j;
//					System.out.println("i=" + i + " j=" + j + "  num=" + nresA[i]);					
					break OUT1;
				}
			}
		}				
//		for(i=1;i<seq1A.length();i++){
//			 System.out.println(i + "\t" + xa[i] + "\t" + ya[i] + "\t" + za[i]);
//		}		
		n_ali = k; // number of aligned residues
		if (n_ali < 1) {
			System.out.println("There is no common residues in the input structures");
			System.exit(1);
		}

		if (nseqB > 15) {
			d0 = 1.24 * Math.pow((nseqB - 15), (1.0 / 3.0)) - 1.8;
		} else {
			d0 = 0.5;
		}
		if (d0 < 0.5) {
			d0 = 0.5;
		}
		if (m_fix == 1) {
			d0 = d0_fix;
		}
		// *** d0_search ----->
		d0_search = d0;
		if (d0_search > 8){
			d0_search = 8;
		}
		if (d0_search < 4.5){
			d0_search = 4.5;
		}
		//***   iterative parameters ----->
		int n_it = 20; // !maximum number of iterations
		d_output = 5; // !for output alignment
		if (m_fix == 1) {
			d_output = d0_fix;
		}
		int n_init_max = 6; // !maximum number of L_init
		int n_init = 0;
		int L_ini_min = 4;
		if (n_ali < 4) {
			L_ini_min = n_ali;
		}
		boolean flag1 = false;
outer:for (i = 1; i <= (n_init_max - 1); i++) {
			n_init = n_init + 1;
			L_ini[n_init] = n_ali/(int)Math.pow(2,(n_init-1));
//			System.out.println("n_ali="+n_ali);
//			System.out.println("n_init="+n_init);
//			System.out.println("TTTTTTT="+L_ini[n_init]);			
			if (L_ini[n_init] <= L_ini_min) {
				L_ini[n_init] = L_ini_min;
				flag1 = true;
				break outer;
			}
		}
		if (flag1 == false) {
			n_init = n_init + 1;
			L_ini[n_init] = L_ini_min;
		}

		score_max = 0; // !TM-score
		double score_maxsub_max = -1; // !MaxSub-score
		double score10_max = -1; // !TM-score10
		double n_GDT05_max = -1; // !number of residues<0.5
		double n_GDT1_max = -1; // !number of residues<1
		double n_GDT2_max = -1; // !number of residues<2
		double n_GDT4_max = -1; // !number of residues<4
		double n_GDT8_max = -1; // !number of residues<8
		int LL,ka,i_init,L_init,iL_max,iL;
		
		for (i_init = 1; i_init <= n_init; i_init++) {  // 333
			L_init = L_ini[i_init];
			iL_max = n_ali - L_init + 1;
			for (iL = 1; iL <= iL_max; iL++) { // 300 !on aligned residues,
				LL = 0;
				ka = 0;
				for (i = 1; i <= L_init; i++) {
					k = iL + i - 1; // ![1,n_ali] common
					r_1[1][i] = xa[iA[k]];
					r_1[2][i] = ya[iA[k]];
					r_1[3][i] = za[iA[k]];
					r_2[1][i] = xb[iB[k]];
					r_2[2][i] = yb[iB[k]];
					r_2[3][i] = zb[iB[k]];
					ka = ka + 1;
					k_ali[ka] = k;
					LL = LL + 1;
					// System.out.println("c1="+r_1[1][i] + " iA[k]=" + iA[k]);
				}
				rms = u3b(r_1, r_2, LL, 1,rms,u,t,ier); // !u rotate r_1 to r_2
				if (i_init==1) { // then !global superposition
					armsd = Math.sqrt(rms/LL);
					rmsd_ali = armsd;
//					System.out.println("AAAAAAAAAAA= rms=" + rms);
//					System.out.println("BBBBBBBBBBB= armsd=" + armsd);
//					System.exit(1);
				}
				for (j = 1; j <=nseqA; j++) {
					xt[j] = t[1] + u[1][1] * xa[j] + u[1][2] * ya[j] + u[1][3]* za[j];
					yt[j] = t[2] + u[2][1] * xa[j] + u[2][2] * ya[j] + u[2][3]* za[j];
					zt[j] = t[3] + u[3][1] * xa[j] + u[3][2] * ya[j] + u[3][3]* za[j];
				}
				d = d0_search - 1;
				score_fun(xt,yt,zt,xb,yb,zb); // !init, get scores, n_cut+i_ali(i) for
				// iteration
				if (score_max < score) {
					score_max = score;
					ka0 = ka;
					for (i = 1; i <= ka0; i++) {
						k_ali0[i] = k_ali[i];
					}
				}
				if (score10_max < score10)
					score10_max = score10;
				if (score_maxsub_max < score_maxsub)
					score_maxsub_max = score_maxsub;
				if (n_GDT05_max < n_GDT05)
					n_GDT05_max = n_GDT05;
				if (n_GDT1_max < n_GDT1)
					n_GDT1_max = n_GDT1;
				if (n_GDT2_max < n_GDT2)
					n_GDT2_max = n_GDT2;
				if (n_GDT4_max < n_GDT4)
					n_GDT4_max = n_GDT4;
				if (n_GDT8_max < n_GDT8)
					n_GDT8_max = n_GDT8;
				// *** iteration for extending
				// ---------------------------------->
				d = d0_search + 1;
		outer_301:		for (int it = 1; it <= n_it; it++) {  // 301
					LL = 0;
					ka = 0;
					for (i = 1; i <= n_cut; i++) {
						int m = i_ali[i]; // ![1,n_ali]
						r_1[1][i] = xa[iA[m]];
						r_1[2][i] = ya[iA[m]];
						r_1[3][i] = za[iA[m]];
						r_2[1][i] = xb[iB[m]];
						r_2[2][i] = yb[iB[m]];
						r_2[3][i] = zb[iB[m]];
						ka = ka + 1;
						k_ali[ka] = m;
						LL = LL + 1;
					}
					rms = u3b(r_1, r_2, LL, 1, rms, u, t, ier); // !u rotate r_1 to r_2
					for (j = 1; j <= nseqA; j++) {
						xt[j] = t[1] + u[1][1] * xa[j] + u[1][2] * ya[j] + u[1][3] * za[j];
						yt[j] = t[2] + u[2][1] * xa[j] + u[2][2] * ya[j] + u[2][3] * za[j];
						zt[j] = t[3] + u[3][1] * xa[j] + u[3][2] * ya[j] + u[3][3] * za[j];
					}
					score_fun(xt,yt,zt,xb,yb,zb); // !init, get scores, n_cut+i_ali(i) for iteration
					if (score_max < score) {
						score_max = score;
						ka0 = ka;
						for (i = 1; i <= ka; i++) {
							k_ali0[i] = k_ali[i];
						}
					}
					if (score10_max < score10)
						score10_max = score10;
					if (score_maxsub_max < score_maxsub)
						score_maxsub_max = score_maxsub;
					if (n_GDT05_max < n_GDT05)
						n_GDT05_max = n_GDT05;
					if (n_GDT1_max < n_GDT1)
						n_GDT1_max = n_GDT1;
					if (n_GDT2_max < n_GDT2)
						n_GDT2_max = n_GDT2;
					if (n_GDT4_max < n_GDT4)
						n_GDT4_max = n_GDT4;
					if (n_GDT8_max < n_GDT8)
						n_GDT8_max = n_GDT8;
					if (it == n_it) {// goto 302
						break outer_301;
					}
					if (n_cut == ka) { // then
							int neq = 0;
							for (i = 1; i <= n_cut; i++) {
								if (i_ali[i] == k_ali[i]){
									neq = neq + 1;
								}
							}
							if (n_cut == neq){ // goto 302
								break outer_301; 
							}
					}
		        }// 301 continue !for iteration
				// 302 continue
			} // 300 continue !for shift
		} // 333 continue !for initial length, L_ali/M		
		DecimalFormat df = new DecimalFormat("0.0000"); 
		//***   output TM-scale ---------------------------->
//		 System.out.println(
//		 "*****************************************************************************\n" +
//	     "*                                 TM-SCORE                                  *\n"  +
//	     "* A scoring function to assess the similarity of protein structures         *\n"  +
//	     "* Based on statistics:                                                      *\n"  +
//	     "* 0.0 < TM-score < 0.17, random structural similarity                       *\n"  +
//	     "* 0.5 < TM-score < 1.00, in about the same fold                             *\n"  +
//	     "*                                                                           *\n"  +
//	     "* Reference: Yang Zhang and Jeffrey Skolnick, Proteins 2004 57: 702-710     *\n"  +
//	     "* For comments, please email to: zhng@umich.edu                             *\n"  +
//	     "*****************************************************************************");
//		 System.out.printf("Structure1: %s  Length= %d\n","PDB1",nseqA);
//		 System.out.printf("Structure1: %s  Length= %d (by which all scores are normalized)\n","PDB2",nseqB);
//		 System.out.printf("Number of residues in common= %d\n",n_ali);
//		 System.out.printf("RMSD of  the common residues= %s\n",df.format(rmsd_ali));
//		 System.out.printf("TM-score    = %s  (d0=%s  TM10=%s)\n",df.format(score_max),df.format(d0),df.format(score10_max));
//		 System.out.printf("MaxSub-score= %s  (d0= 3.50)\n",df.format(score_maxsub_max));		 
	     double score_GDT=(n_GDT1_max+n_GDT2_max+n_GDT4_max+n_GDT8_max)/(float)(4*nseqB);
//	     System.out.printf("GDT-TS-score= %s (d<1)= %s  (d<2)= %s (d<4)= %s (d<8)= %s\n",
//	    		 df.format(score_GDT),df.format(n_GDT1_max/(float)(nseqB)),df.format(n_GDT2_max/(float)(nseqB)),df.format(n_GDT4_max/(float)(nseqB)),df.format(n_GDT8_max/(float)(nseqB)));
	   
	     score_GDT_HA=(n_GDT05_max+n_GDT1_max+n_GDT2_max+n_GDT4_max)/(float)(4*nseqB);
//	     System.out.printf("GDT-HA-score= %s (d<0.5)= %s  (d<1)= %s (d<2)= %s, (d<4)= %s\n\n",
//	    		 df.format(score_GDT_HA),df.format(n_GDT05_max/(float)(nseqB)),df.format(n_GDT1_max/(float)(nseqB)),df.format(n_GDT2_max/(float)(nseqB)),df.format(n_GDT4_max/(float)(nseqB)));
	     
		// *** recall and output the superposition of maxiumum TM-score:
		LL = 0;
		for (i = 1; i <= ka0; i++) {
			int m = k_ali0[i]; // !record of the best alignment
			r_1[1][i] = xa[iA[m]];
			r_1[2][i] = ya[iA[m]];
			r_1[3][i] = za[iA[m]];
			r_2[1][i] = xb[iB[m]];
			r_2[2][i] = yb[iB[m]];
			r_2[3][i] = zb[iB[m]];
			LL = LL + 1;
		}
		rms = u3b(r_1,r_2,LL,1,rms,u,t,ier); // !u rotate r_1 to r_2
		for (j = 1; j <= nseqA; j++) {
			xt[j] = t[1] + u[1][1] * xa[j] + u[1][2] * ya[j] + u[1][3] * za[j];
			yt[j] = t[2] + u[2][1] * xa[j] + u[2][2] * ya[j] + u[2][3] * za[j];
			zt[j] = t[3] + u[3][1] * xa[j] + u[3][2] * ya[j] + u[3][3] * za[j];
		}

		// ********* extract rotation matrix ------------>
		// write(*,*)'-------- rotation matrix to rotate Chain-1 to ',
		//System.out.println(" -------- rotation matrix to rotate Chain-1 to Chain-2 ------");
		//System.out.printf("%-6s %-10s %-10s %-10s %-10s\n","i","t(i)","u(i,1)","u(i,2)","u(i,3)");
		for (i = 1; i <= 3; i++) {
			//System.out.printf("%-6d %-10s %-10s %-10s %-10s\n",i,df.format(t[i]),df.format(u[i][1]),df.format(u[i][2]),df.format(u[i][3]));		
			}
		// ********* rmsd in superposed regions --------------->
		d = d_output; // !for output
		score_fun(xt,yt,zt,xb,yb,zb);  // !give i_ali(i), score_max=score now
		LL = 0;
		for (i = 1; i <= n_cut; i++) {
			int m = i_ali[i]; // ![1,nseqA]
			r_1[1][i] = xa[iA[m]];
			r_1[2][i] = ya[iA[m]];
			r_1[3][i] = za[iA[m]];
			r_2[1][i] = xb[iB[m]];
			r_2[2][i] = yb[iB[m]];
			r_2[3][i] = zb[iB[m]];
			LL = LL + 1;
		}
		rms = u3b(r_1,r_2,LL,0,rms,u,t,ier);
		armsd = Math.sqrt(rms / LL);
		rmsd = armsd;
		// ***   output rotated chain1 + chain2----->
		if(m_out==1){
			String temp;
			FileWriter fw = new FileWriter(outname);   //  !pdb1.aln + pdb2.aln			
			fw.append("load inline\n"         +
					  "select atomno<1000\n"  +
					  "wireframe .45\n"       +
					  "select none\n"         +
					  "select atomno>1000\n"  +
					  "wireframe .15\n"       +
					  "color white\n");
			for(i=1;i<=n_cut;i++){
				temp = String.format("select %4d",nresA[iA[i_ali[i]]]);
				fw.append(temp);
				fw.append("\ncolor red\n");
			}
			
			fw.append("select all\nexit\n");
			temp = String.format("REMARK  RMSD of the common residues=%8.3f\n", rmsd_ali);
			fw.append(temp);
			temp = String.format("REMARK  TM-score=%6.4f (d0=%5.2f)\n",score_max,d0);
			fw.append(temp);
			for(i=1;i<=nseqA;i++){
				temp = String.format("ATOM  %5d  CA  %3s%6d    %8.3f%8.3f%8.3f\n",nresA[i],NameMap(String.valueOf(seq1A.charAt(i))),nresA[i],xt[i],yt[i],zt[i]);				
				fw.append(temp);
			}			
			fw.append("TER\n");
			for(i=2;i<=nseqA;i++){
				temp = String.format("CONECT%5d%5d\n",nresA[i-1],nresA[i]);
				fw.append(temp);
			}
			for(i=1;i<=nseqB;i++){
				temp = String.format("ATOM  %5d  CA  %3s%6d    %8.3f%8.3f%8.3f\n",2000+nresB[i],NameMap(String.valueOf(seq1B.charAt(i))),nresB[i],xb[i],yb[i],zb[i]);				
				fw.append(temp);
			}			
			fw.append("TER\n");			
			for(i=2;i<=nseqB;i++){
				temp = String.format("CONECT%5d%5d\n",2000+nresB[i-1],2000+nresB[i]);
				fw.append(temp);
			}			
			fw.flush();
			fw.close();
			fw = null;
		}		
		
		// *** record aligned residues by i=[1,nseqA], for
		// sequenceM()------------>
		for (i = 1; i <= nseqA; i++) {
			iq[i] = 0;
		}
		double dis = 0;
		for (i = 1; i <= n_cut; i++) {
			j = iA[i_ali[i]]; // ![1,nseqA]
			k = iB[i_ali[i]]; // ![1,nseqB]
			dis = Math.sqrt((xt[j] - xb[k]) * (xt[j] - xb[k])
					+ (yt[j] - yb[k]) * (yt[j] - yb[k]) + (zt[j] - zb[k])
					* (zt[j] - zb[k]));
			if (dis < d_output) {
				iq[j] = 1;
			}
		}

		// *******************************************************************
		// *** output aligned sequences

		k = 0;
		i = 1;
		j = 1;
		String sequenceA = "";
		String sequenceB = "";
		String sequenceM = "";
		
		
//		System.out.println(seq1A.charAt(nseqA));
//		System.out.println(seq1B.charAt(nseqB));
//		System.exit(1);
//System.out.println("length=" + nseqA + " " + nseqB);
		
ali_while:while (true) {
			if (i > nseqA && j > nseqB)
				break ali_while;

			else if (i > nseqA && j <= nseqB) {
				sequenceA += '-';
				sequenceB += seq1B.charAt(j);
				sequenceM += ' ';
				j = j + 1;
			}
			else if (i <= nseqA && j > nseqB) {
				sequenceA += seq1A.charAt(i);
				sequenceB += '-';
				sequenceM += ' ';
				i = i + 1;
			}
			else if (nresA[i] == nresB[j]) {
				sequenceA += seq1A.charAt(i);
				sequenceB += seq1B.charAt(j);
				if (iq[i] == 1) {
					sequenceM += ':';
				} else {
					sequenceM += ' ';
				}
				i = i + 1;
				j = j + 1;
			} 
			else if(nresA[i] < nresB[j]) {
				//if(i<nseqA){
					if(i<nseqA){
						sequenceA += seq1A.charAt(i);
					}else{
						sequenceA += ' ';
					}
					sequenceB += '-';
					sequenceM += ' ';
					i = i + 1;
				//}
			} else if (nresA[i]>nresB[j]) {
				//if(j<nseqB){
					sequenceA += '-';
					if(j<nseqB){
						sequenceB += seq1B.charAt(j);
					}else{
						sequenceB += ' ';
					}
					sequenceM += ' ';
					j = j + 1;
				//}
			}
		}		
		//System.out.printf("\nSuperposition in the TM-score: Length(d<%f)= %d  RMSD=  %f\n",d_output,n_cut,rmsd);
		//System.out.printf(": denotes the residue pairs of distance  %s Angstrom\n",df.format(d_output));		
		//System.out.println(sequenceA);
		//System.out.println(sequenceM);
		//System.out.println(sequenceB);		
		for(i=1;i<=sequenceA.length();i++){
			//System.out.print(i%10);
		}
		//System.out.println("\n");
		release_memory();
	}

	public String readPDB(String file, double[] x, double[] y,double[] z, int[] nres){
		String seq = "*";
		String line = "";
		BufferedReader br;
		int i;
		i = 0;
		try {
			br = new BufferedReader(new FileReader(file));
			while ((line = br.readLine()) != null) {
				if (line.startsWith("TER"))
					break;
				if (line.startsWith("ATO")) {
					if (line.substring(13, 16).replaceAll("\\s+", "").endsWith(
							"CA")) {
						i++;
						//seq = seq + NameMap(line.substring(17, 20).toUpperCase());
						seq = seq + "X";
						x[i] = Float.valueOf(line.substring(30, 38));
						y[i] = Float.valueOf(line.substring(38, 46));
						z[i] = Float.valueOf(line.substring(46, 54));
						nres[i] = Integer.valueOf(line.substring(22, 29).replaceAll("\\s+", ""));
						//System.out.println(nres[i]);
					}
				}
			}
			br.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return seq;
	}

	public String NameMap(String residule) {
		String[] aa = new String[] { "ALA", "ARG", "ASN", "ASP", "CYS", "GLN",
				"GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO",
				"SER", "THR", "TRP", "TYR", "VAL", "ASX", "GLX", "UNK" };
		String[] aaName = new String[] { "A", "R", "N", "D", "C", "Q", "E",
				"G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y",
				"V", "B", "Z", "X" };
		int i;
		for (i=0; i < aa.length; i++) {
			if (aa[i].equals(residule))
				return aaName[i];
		}	
		for (i=0; i < aaName.length; i++) {
			if (aaName[i].equals(residule))
				return aa[i];
		}	
		return "error";
	}

	// 1, collect those residues with dis<d;
	// 2, calculate TMscore
	void score_fun(double[] xt,double[] yt,double[] zt,double[] xb,double[] yb,double[] zb) {
		double d_tmp = d;		
		double score_maxsub_sum=0;       // !Maxsub-score
	    double score_sum=0;              // !TMscore
	    double score_sum10=0;            // !TMscore10
	    double dis = 0;
	    
w_outer:while(true){
		n_cut=0;                  // !number of residue-pairs dis<d, for iteration
	    n_GDT05=0;                // !for GDT-score, # of dis<0.5
	    n_GDT1=0;                 // !for GDT-score, # of dis<1
	    n_GDT2=0;                 // !for GDT-score, # of dis<2
	    n_GDT4=0;                 // !for GDT-score, # of dis<4
	    n_GDT8=0;                 // !for GDT-score, # of dis<8
	    score_maxsub_sum=0;       // !Maxsub-score
	    score_sum=0;              // !TMscore
	    score_sum10=0;            // !TMscore10
	    int i,j;
	    dis = 0;	    
	    for(int k=1;k<=n_ali;k++){
	    	i=iA[k];              //  ![1,nseqA] reoder number of structureA
            j=iB[k];              //  ![1,nseqB]
            dis=Math.sqrt((xt[i]-xb[j])*(xt[i]-xb[j])+(yt[i]-yb[j])*(yt[i]-yb[j])+(zt[i]-zb[j])*(zt[i]-zb[j]));
	    //***   for iteration:
	        if(dis<d_tmp){
	           n_cut=n_cut+1;
	           i_ali[n_cut]=k;    //  ![1,n_ali], mark the residue-pairs in dis<d
	        }	        
		//***   for GDT-score:
	        //Tuan added
		        if(dis<=gdt8){
		           n_GDT8=n_GDT8+1;
		           if(dis<=gdt4){
		              n_GDT4=n_GDT4+1;
		              if(dis<=gdt2){
		                 n_GDT2=n_GDT2+1;
		                 if(dis<=gdt1){
		                    n_GDT1=n_GDT1+1;
		                    if(dis<=gdt05){
		                       n_GDT05=n_GDT05+1;
		                    }
		                 }
		              }
		           }
		        }
        
		//***   for MAXsub-score:
		        if(dis<3.5){
		           score_maxsub_sum=score_maxsub_sum+1/(1+(dis/3.5)*(dis/3.5));
		        }
		//***   for TM-score:
		        score_sum=score_sum+1/(1+(dis/d0)*(dis/d0));
		//***   for TM-score10:
		        if(dis<10){
		           score_sum10=score_sum10+1/(1+(dis/d0)*(dis/d0));
		        }        
	    } // for	    
        if(n_cut<3&&n_ali>3){
        	d_tmp=d_tmp+0.5;
        }else{
        	break w_outer;
        }        
	 } // while(true);        
     score_maxsub=score_maxsub_sum/(float)(nseqB); // !MAXsub-score
     score=score_sum/(float)(nseqB);               // !TM-score
     score10=score_sum10/(float)(nseqB);           // !TM-score10	 
	}

	public double u3b(double[][] x, double[][] y, int n, int mode,
			double rms, double[][] u, double[] t, int ier) {		
//		cccccccccccccccc Calculate sum of (r_d-r_m)^2 cccccccccccccccccccccccccc
//		c  w    - w(m) is weight for atom pair  c m                 (given)
//		c  x    - x(i,m) are coordinates of atom c m in set x       (given)
//		c  y    - y(i,m) are coordinates of atom c m in set y       (given)
//		c  n    - n is number of atom pairs                         (given)
//		c  mode  - 0:calculate rms only                             (given)
//		c          1:calculate rms,u,t                              (takes longer)
//		c  rms   - sum of w*(ux+t-y)**2 over all atom pairs         (result)
//		c  u    - u(i,j) is   rotation  matrix for best superposition  (result)
//		c  t    - t(i)   is translation vector for best superposition  (result)
//		c  ier  - 0: a unique optimal superposition has been determined(result)
//		c       -1: superposition is not unique but optimal
//		c       -2: no result obtained because of negative weights w
//		c           or all weights equal to zero.
//		cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
		int i, j, m, m1, l, k;
		double e0,d, h, g;
		double cth, sth, sqrth, p, det, sigma;
		double[] xc = new double[3+1];
		double[] yc = new double[3+1];
		double[][] a = new double[3+1][3+1];
		double[][] b = new double[3+1][3+1];
		double[][] r = new double[3+1][3+1];
		double[] e = new double[3+1];
		double[] rr = new double[6+1];
		double[] ss = new double[6+1];
		double sqrt3 = 1.73205080756888, tol = 0.01;
		int ip[] = {-100,1,2, 4, 2, 3, 5, 4, 5, 6 };
		int ip2312[] = {-100, 2, 3, 1, 2};
		int a_failed = 0, b_failed = 0;
		double epsilon = 0.000000001;

		// initializtation
		rms = 0;
		e0 = 0;		
		for(i = 1; i<=3; i++) {
			xc[i] = 0.0;
			yc[i] = 0.0;
			t[i] = 0.0;
			for (j = 1; j <= 3; j++) {
				u[i][j] = 0.0;
				r[i][j] = 0.0;
				a[i][j] = 0.0;
				if (i == j) {
					u[i][j] = 1.0;
					a[i][j] = 1.0;
				}
			}
		}

		ier = -1;
		if (n < 1) {
			return -1000;
		}
		ier = -2;		
		// compute centers for vector sets x, y
		for(m=1;m<=n;m++){        
			for(i=1;i<=3;i++){
				xc[i] = xc[i] + x[i][m];
				yc[i] = yc[i] + y[i][m];
			}
		}		
		for (i = 1; i <= 3; i++) {
			xc[i] = xc[i] / n;
			yc[i] = yc[i] / n;
		}
		// compute e0 and matrix r
		for (m = 1; m <= n; m++) {
			for (i = 1; i <= 3; i++) {
				e0 = e0 + (x[i][m] - xc[i])*(x[i][m] - xc[i]) + (y[i][m] - yc[i])*(y[i][m] - yc[i]);
				d = y[i][m] - yc[i];
				for (j = 1; j <= 3; j++) {
					r[i][j] = r[i][j] + d * (x[j][m] - xc[j]);
				}
			}
		}
		// compute determinat of matrix r
		 det = r[1][1] * ( (r[2][2]*r[3][3]) - (r[2][3]*r[3][2]) )
	         - r[1][2] * ( (r[2][1]*r[3][3]) - (r[2][3]*r[3][1]) )
	         + r[1][3] * ( (r[2][1]*r[3][2]) - (r[2][2]*r[3][1]) );
		 
		sigma = det;
		// compute tras(r)*r
		m = 0;
		for (j = 1; j <= 3; j++) {
			for (i = 1; i <= j; i++) {
				m = m+1;
	            rr[m] = r[1][i]*r[1][j] + r[2][i]*r[2][j] + r[3][i]*r[3][j];
			}
		}

		double spur = (rr[1] + rr[3] + rr[6]) / 3.0;
		double   cof = (((((rr[3]*rr[6] - rr[5]*rr[5]) + rr[1]*rr[6])
			          - rr[4]*rr[4]) + rr[1]*rr[3]) - rr[2]*rr[2]) / 3.0;
		det = det * det;

		for (i = 1; i <= 3; i++) {
			e[i] = spur;
		}		
//		System.out.println("spur="+spur);
//		System.exit(1);
		if (spur > 0) {
			d = spur * spur;
			h = d - cof;
			g = (spur * cof - det) / 2.0 - spur * h;
			if (h > 0) {
				sqrth = Math.sqrt(h);
				d = h * h * h - g * g;
				if (d < 0.0){
				   d = 0.0;
				}
				d = Math.atan2(Math.sqrt(d), -g) / 3.0;
				cth = sqrth * Math.cos(d);
				sth = sqrth * sqrt3 * Math.sin(d);
				e[1] = (spur + cth) + cth;
				e[2] = (spur - cth) + sth;
				e[3] = (spur - cth) - sth;
				if (mode != 0) {// compute a
					for (l = 1; l <= 3; l = l + 2) {
						d = e[l];
						ss[1] = (d - rr[3]) * (d - rr[6]) - rr[5] * rr[5];
						ss[2] = (d - rr[6]) * rr[2] + rr[4] * rr[5];
						ss[3] = (d - rr[1]) * (d - rr[6]) - rr[4] * rr[4];
						ss[4] = (d - rr[3]) * rr[4] + rr[2] * rr[5];
						ss[5] = (d - rr[1]) * rr[5] + rr[2] * rr[4];
						ss[6] = (d - rr[1]) * (d - rr[3]) - rr[2] * rr[2];

						if (Math.abs(ss[0]) <= epsilon)
							ss[0] = 0.0;
						if (Math.abs(ss[1]) <= epsilon)
							ss[1] = 0.0;
						if (Math.abs(ss[2]) <= epsilon)
							ss[2] = 0.0;
						if (Math.abs(ss[3]) <= epsilon)
							ss[3] = 0.0;
						if (Math.abs(ss[4]) <= epsilon)
							ss[4] = 0.0;
						if (Math.abs(ss[5]) <= epsilon)
							ss[5] = 0.0;

						if (Math.abs(ss[1]) >= Math.abs(ss[3])) {
							j = 1;
							if (Math.abs(ss[1]) < Math.abs(ss[6])) {
								j = 3;
							}
						} else if (Math.abs(ss[3]) >= Math.abs(ss[6])) {
							j = 2;
						} else {
							j = 3;
						}

						d = 0.0;
						j = 3 * (j-1);
						for (i = 1; i <= 3; i++) {
							k = ip[i + j];
							a[i][l] = ss[k];
							d = d + ss[k] * ss[k];
						}
						// if( d > 0.0 ) d = 1.0 / sqrt(d);
						if (d > 0)
							d = 1.0 / Math.sqrt(d);
						else
							d = 0.0;
						for (i = 1; i <= 3; i++) {
							a[i][l] = a[i][l] * d;
						}
					}// for l

					d = a[1][1]*a[1][3] + a[2][1]*a[2][3] + a[3][1]*a[3][3];
					
					if ((e[1] - e[2]) > (e[2] - e[3])){
						m1 = 3;
						m = 1;
					}else{
						m1 = 1;
						m = 3;
					}			      
					p = 0;					
					for(i=1;i<=3;i++){
						a[i][m1] = a[i][m1] - d*a[i][m];
						p = p + a[i][m1]*a[i][m1];
					}
					
					if (p <= tol) {
						p = 1.0;
				KKK:		for (i = 1; i <= 3; i++) {
								if (p < Math.abs(a[i][m])) {
									continue KKK;
								}
								p = Math.abs(a[i][m]);
								j = i;
						}
						k = ip2312[j];
						l = ip2312[j + 1];
						p = Math.sqrt(a[k][m] * a[k][m] + a[l][m] * a[l][m]);
						if (p > tol) {
							a[j][m1] = 0.0;
							a[k][m1] = -a[l][m] / p;
							a[l][m1] = a[k][m] / p;
						} else {// goto 40
							a_failed = 1;
						}
					}// if p<=tol
					else {
						p = 1.0 / Math.sqrt(p);
						for (i = 1; i <= 3; i++) {
							a[i][m1] = a[i][m1] * p;
						}
					}// else p<=tol
					if (a_failed != 1) {
						  a[1][2] = a[2][3]*a[3][1] - a[2][1]*a[3][3];
					      a[2][2] = a[3][3]*a[1][1] - a[3][1]*a[1][3];
					      a[3][2] = a[1][3]*a[2][1] - a[1][1]*a[2][3];
					}
				}// if(mode!=0)
			}// h>0

			// compute b anyway
			if (mode != 0 && a_failed != 1){// a is computed correctly
				// compute b
				for (l = 1; l <= 2; l++) {
					d = 0.0;
					for (i = 1; i <= 3; i++) {
						b[i][l] = r[i][1]*a[1][l] + r[i][2]*a[2][l] + r[i][3]*a[3][l];
						d = d + b[i][l] * b[i][l];
					}
					// if( d > 0 ) d = 1.0 / sqrt(d);
					if (d > 0)
						d = 1.0 / Math.sqrt(d);
					else
						d = 0.0;
					for (i = 1; i <= 3; i++) {
						b[i][l] = b[i][l] * d;
					}
				}
				d = b[1][1]*b[1][2] + b[2][1]*b[2][2] + b[3][1]*b[3][2];
				p = 0.0;

				for(i=1;i<=3;i++){
					b[i][2] = b[i][2] - d*b[i][1];
					p = p + b[i][2]*b[i][2];
				}

				if (p <= tol) {
					p = 1.0;
			Line22:	for (i = 1; i <= 3; i++) {
						if (p < Math.abs(b[i][1])) {
							continue Line22;
						}
						p = Math.abs(b[i][1]);
						j = i;
					}
					k = ip2312[j];
					l = ip2312[j + 1];
					p = Math.sqrt(b[k][1] * b[k][1] + b[l][1] * b[l][1]);
					if (p > tol) {
						b[j][2] = 0.0;
						b[k][2] = -b[l][1] / p;
						b[l][2] = b[k][1] / p;
					} else {
						// goto 40
						b_failed = 1;
					}
				}// if( p <= tol )
				else {
					p = 1.0 / Math.sqrt(p);
					for (i = 1; i <= 3; i++) {
						b[i][2] = b[i][2] * p;
					}
				}
				if (b_failed != 1) {
					 b[1][3] = b[2][1]*b[3][2] - b[2][2]*b[3][1];
				     b[2][3] = b[3][1]*b[1][2] - b[3][2]*b[1][1];
				     b[3][3] = b[1][1]*b[2][2] - b[1][2]*b[2][1];
					// compute u
				     for(i=1;i<=3;i++){
				    	 for(j=1;j<=3;j++){
			                u[i][j] = b[i][1]*a[j][1] + b[i][2]*a[j][2] + b[i][3]*a[j][3];
				    	 }
				     }
				}

				// compute t
				for(i=1;i<=3;i++){
					t[i] = ((yc[i] - u[i][1]*xc[1]) - u[i][2]*xc[2]) - u[i][3]*xc[3];
				}
			}// if(mode!=0 && a_failed!=1)
		}// spur>0
		else // just compute t and errors
		{
			// compute t
			for(i=1;i<=3;i++){
				t[i] = ((yc[i] - u[i][1]*xc[1]) - u[i][2]*xc[2]) - u[i][3]*xc[3];
			}
		}// else spur>0

		// compute rms
		for (i = 1; i <= 3; i++) {
			if (e[i] < 0)
				e[i] = 0;
			e[i] = Math.sqrt(e[i]);
		}
		ier = 0;
	    if( e[2] <= (e[1] * 1.0e-05) ) ier = -1;	    
		d = e[3];
		if (sigma < 0.0) {
			d = -d;
			if( (e[2] - e[3]) <= (e[1] * 1.0e-05) ) ier = -1;
		}
		d = (d + e[2]) + e[1];
		rms = (e0 - d) - d;
		if (rms < 0.00000000001)
			rms = 0.0;
		return rms;
	}
	
	public String help(){
        String help="Brief instruction for running TM-score program:\n" +
		 "(For detail: Zhang & Skolnick,  Proteins, 2004 57:702-10)\n"  +
        "                                                          \n"  +
		 "1. Run TM-score to compare 'model' and 'native:          \n"  +
		 "   >java -jar TMscore.jar model native                   \n"  +                                   
		 "                                                         \n"  +
		 "2. Run TM-score with an assigned d0, e.g. 5 Angstroms:   \n"  +
		 "   >java -jar TMscore.jar model native -d 5              \n"  +
         "                                                         \n"  +
		 "3. Run TM-score with superposition output, e.g. 'TM.sup':\n"  +
		 "   >java -jar TMscore.jar model native -o TM.sup         \n"  +
		 "   To view the superimposed structures by rasmol:        \n"  +
		 "   >rasmol -script TM.sup                                \n";	         
        return help;
	}
	
	public void release_memory(){
		L_ini = null;
		iA = null;
		iB = null;
		i_ali = null;
		nresA = null;
		nresB = null;	
		u = null;
		t = null;
		iq = null;
		xa = null;
		ya = null;
		za = null;
		xt = null;
		yt = null;
		zt = null;
		xb = null;
		yb = null;
		zb = null;		
		r_1 = null;
		r_2 = null;
	}

	public double getScore() {
		return score;
	}

	public void setScore(double score) {
		this.score = score;
	}

	public double getScore_maxsub() {
		return score_maxsub;
	}

	public void setScore_maxsub(double score_maxsub) {
		this.score_maxsub = score_maxsub;
	}

	public double getScore_fix() {
		return score_fix;
	}

	public void setScore_fix(double score_fix) {
		this.score_fix = score_fix;
	}

	public double getScore10() {
		return score10;
	}

	public void setScore10(double score10) {
		this.score10 = score10;
	}

	public double getScore_max() {
		return score_max;
	}

	public void setScore_max(double score_max) {
		this.score_max = score_max;
	}

	public double getScore_fix_max() {
		return score_fix_max;
	}

	public void setScore_fix_max(double score_fix_max) {
		this.score_fix_max = score_fix_max;
	}

	public double getArmsd() {
		return armsd;
	}

	public void setArmsd(double armsd) {
		this.armsd = armsd;
	}

	public double getScore_GDT_HA() {
		return score_GDT_HA;
	}

	public void setScore_GDT_HA(double score_GDT_HA) {
		this.score_GDT_HA = score_GDT_HA;
	}
	
}
