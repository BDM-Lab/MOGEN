################################################################################
MOGEN: The software for modeling the 3D structure of a genome using Hi-C chromosome conformation capturing data

Bioinformatics, Data Mining, Machine Learning (BDM) Laboratory 
         University of Missouri, Columbia MO 65211
################################################################################

Developer:
		Tuan Trieu
		Department of Computer Science
		University of Missouri, Columbia
		Email: tuantrieu@mail.missouri.edu

Contact:
		Jianlin Cheng, PhD
		Department of Computer Science
		University of Missouri, Columbia
		Email: chengji@missouri.edu

#################################################################################
1. Content of folders:

- bin: contains executable MOGEN and normalizer (to normalize the input data)
- example: contains example data and parameter files used to reconstruct chromosome/genome structures
- src: source code of MOGEN in java

#################################################################################
2. USAGE

To run the tool, type: java -jar 3DGenorator.jar parameters_normal.txt

The file parameters_normal.txt contains parameters needed to run the tool

See in /examples/hiC/ for sample files

#################################################################################

3. OUTPUT
#################################################################################

MOGEN produces two types of output files: (1) 3D models in PDB format and (2) text files containing contact and non-contact scores of each chromosome in the models

################################################################################
5. DISCLAIMER
################################################################################

The executable software and the source code of UniCon3D is distributed free of 
charge as it is to any non-commercial users. The authors hold no liabilities to 
the performance of the program.


