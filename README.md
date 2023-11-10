# GeneMarkHmmEukaryotic3
GeneMark.hmm for eukaryotic genomes (with introns) version 3  

Copyright (c) 2001, GeorgiaTech  
Copyright (c) 2001, Gene Probe, Inc.  
All Rights Reserved  

Algorithm GeneMark.hmm predicts protein coding genes within a eukaryotic DNA sequence.  

Lomsadze A., Ter-Hovhannisyan V., Chernoff Y. and Borodovsky M.  
"Gene identification in novel eukaryotic genomes by self-training algorithm."  
Nucleic Acids Research, 2005 Nov; 33(20):6494-506  

Ter-Hovhannisyan V., Lomsadze A., Chernoff Y. and Borodovsky M.  
"Gene prediction in novel fungal genomes using an ab initio algorithm with unsupervised training."  
Genome Research, 2008 Dec; 18(12):1979-90  

## Software license

### GeneMark.hmm eukaryotic 3

GeneMark.hmm eukaryotic 3 is distributed under Creative Commons Attribution NonCommercial ShareAlike 4.0 License.
See the [full text of the license](License-Creative-Commons-Attribution-NonCommercial-ShareAlike-4.0-International.txt).  

## Usage

```
-----------------------------------------
 ./gmhmme3
GeneMark.hmm eukaryotic, version 3.68.pub

Usage: gmhmme3 [options] <sequence file>

  required parameters:
    -m <model file>

  optional parameters:
    -o <output file>
    -p write protein translation
    -n write nucleotide sequence
    -b <output file> output statistics of predicted introns
    -d <file name> input for GeneMark.hmm plus
    -s <string> sequence tag in GFF output format
    -f <format> output prediction in [lst|gff3|gtf] format; default [lst]
    -k <number> value for soft-mask penalty

developer options:
     -z trace back seqid and position
     -w <number> minimum gene length
     -r report best path probability
     -v verbose
------------------------------------------
```
