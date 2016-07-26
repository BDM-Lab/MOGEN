
----------

## MOGEN: The software for modeling the 3D structure of a genome using Hi-C chromosome conformation capturing data #

----------

#### Bioinformatics, Data Mining, Machine Learning (BDM) Laboratory, 
#### University of Missouri, Columbia MO 65211

----------

#### Developer: <br/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Tuan Trieu <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Department of Computer Science <br/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; University of Missouri, Columbia <br/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Email: tuantrieu@mail.missouri.edu <br/>

#### Contact: <br/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Jianlin Cheng, PhD <br/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Department of Computer Science <br/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; University of Missouri, Columbia <br/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Email: chengji@missouri.edu <br/>





## 1. Content of folders:
- bin: contains executable MOGEN and normalizer (to normalize the input data)
- example: contains example data and parameter files used to reconstruct chromosome/genome structures
- src: source code of MOGEN in java




## 2. USAGE ##

To run the tool, type: `java -jar 3DGenorator.jar parameters_normal.txt`

The file parameters_normal.txt contains parameters needed to run the tool

See in `/examples/hiC/` for sample files


## 3. OUTPUT ##

MOGEN produces two types of output files: 
	
- 3D models in PDB format and 
- text files containing contact and non-contact scores of each chromosome in the models


## 4. DISCLAIMER ##

The executable software and the source code of MOGEN is distributed free of 
charge as it is to any non-commercial users. The authors hold no liabilities to 
the performance of the program.

## 5. Citations
- Trieu T, Cheng J. Large-scale reconstruction of 3D structures of human chromosomes from chromosomal contact data. Nucleic Acids Res. 2014
- Trieu T, Cheng J. MOGEN: a tool for reconstructing 3D models of genomes from chromosomal conformation capturing data. Bioinformatics. 2015

