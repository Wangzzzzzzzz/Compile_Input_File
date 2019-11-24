# Compile Input File 

This Repository contains the code that generate input file of the Deep Learning Model.

## Version

#### 1.1

```diff
+ Add the the count of protein of each flag (0:'Normal',1:'2 unique seq with duplicates',2: 'More than 2 unique seq', 3: 'Less than 2 seq')
```

#### 1.2

```diff
+ Add mutant sequence generater.
```

## File System

| Folder | Description |
| --- | --- |
| dataset | the input datasets for __formatting__ |
| FASTA | fasta sequences of the relevent proteins |
| PDB | pdb files of the relevent proteins |
| output | formatted output, the file ready to be put in to the machine learning model is __formatted.csv__, the count of proteins for each flag is in __Protein_summary.xlsx__, Generated mutant sequence is located in __MU_WD_seq.txt__ |

## Package Info

| Package Name | Version |
| -- | -- |
| Biopython (Bio) | 1.73 |
| pandas | 0.24.2 |
| numpy | 1.16.2 |
