This program is to normalize input data using the coverage normalization or iterative correction of Hi-C

The input contains list of contacts, each line is 3 number: position1 position2 interaction_frequency

To run the program: java -jar Normalizer.jar -i {input file} -o {output file} -m {1: for coverage normalization or 2: for iterative correction of Hi-C}

Example:
$> java -jar Normalizer.jar -i IFList_Chr_1_100K.txt -o IFList_Chr_1_100K_normalized.txt -m 1

***Please note that all input files provided in the input folder of MOGEN are normalized.