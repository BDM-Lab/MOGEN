
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




## 2. Usage ##

To run the tool, type: `java -jar 3DGenorator.jar parameters_normal.txt`

The file parameters_normal.txt contains parameters needed to run the tool

See in `/examples/hiC/` for sample files


## 3. Output ##

MOGEN produces two types of output files: 
	
- 3D models in PDB format and 
- text files containing contact and non-contact scores of each chromosome in the models


## 4. Disclaimer ##

The executable software and the source code of MOGEN is distributed free of 
charge as it is to any non-commercial users. The authors hold no liabilities to 
the performance of the program.

## 5. Citations
- Trieu T, Cheng J. Large-scale reconstruction of 3D structures of human chromosomes from chromosomal contact data. Nucleic Acids Res. 2014
- Trieu T, Cheng J. MOGEN: a tool for reconstructing 3D models of genomes from chromosomal conformation capturing data. Bioinformatics. 2015

## 6. Common questions ##
### 1. How to adjust parameters to get good models ? ###

MOGEN produces two output files, a .pdb file containing a model and a *_evaluation.txt file containing contact and non-contact scores of the model. Parameters are adjusted to maximize contact and non-contact scores. The effect of each parameter on corresponding contact and non-contact scores is mentioned in comments in the parameter sample file.

### 2. Why my non-contact score is NaN ? ###

This is because there is no or very few non-contacts in the input data, which often happens when the resolution is low (e.g 1MB). Set an appropriate contact threshold to make sure that there is at least 20% - 30% non-contacts (out of all possible contacts), the exact number doesn't matter because models are often similar.  

