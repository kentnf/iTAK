# iTAK

iTAK is a bioinformatics software tool for easy and logical identificatoin of
plant trasnscription factors and transcription domains as defined by Lehti-Shiu MD, Shiu S-H (2012).
iTAK utilizes both plnTFDB and plantTFBD as referance databases

Read the paper and cite the paper t http://authors.elsevier.com/sd/article/S1674205216302234

## Getting Started

If you're running linux, iTAK has a very simple setup. The only outside
requirnments are an up to date installation of BioPerl. iTAK comes with a
precompiled version of hmm.

If you're not using linux, then a custom installation of hmm will need to be
done, and a copy of all binaries will need to be placed in the /bin folder of
iTAK

To download iTAK you cna either download it from our website
ftp://bioinfo.bti.cornell.edu/pub/program/itak/old or clone this github repo
using 

```
github clone https://github.com/kentnf/iTAK.git
```

## Test Install

Once bioperl is installed you can test your installation of iTAK by running

```
perl iTAK.pl test_seq 
```

### Running iTAK

iTAK is a relativly simple program to call

Example:

```
perl iTAK.pl Multifasta -o ItakOutput

Parameters:

-f	[String]	Frame used for nucleotide translation.
	6: 6 frame translation (default) 
	3F: 3 frame translation with forward strand 
	3R: 3 frame translation with reverse strand

-a	[Integer]	Number of CPUs used for hmmscan (default = 1)

-o	[String]	Name of the output directory (default = "input file name" + "_output")
```

#### Output
The Output from iTAK will be a list of files 

```
tf_sequence.fasta: sequences of all identified TFs/TRs (FASTA format).

tf_classification.txt: classification of all identified TFs/TRs. A tab-delimited txt file containing sequence IDs and their families.

tf_alignment.txt: A tab-delimited txt file containing alignments of all identified TFs/TRs to the protein domain database.

pk_sequence.fasta: sequences of all identified PKs (FASTA format).

Shiu_classification.txt: classification of all identified protein kinases. A tab-delimited txt file containing sequence IDs and their corresponding protein kinase families.

Shiu_alignment.txt: A tab-delimited txt file containing alignments of all identified protein kinases to the protein domain database.
```

### Contribute

- Issue Tracker: github.com/iTAK/iTAK/issues
- Source Code: github.com/iTAK/iTAK

### Support

If you are having issues, please let us know.
Contact us at http://bioinfo.bti.cornell.edu/cgi-bin/itak/contact.cgi

### License

The project is licensed under the BSD license.
