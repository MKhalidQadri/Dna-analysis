# DNA Analysis Application

## Overview
The DNA Analysis Application is a Python program designed to process multi-FASTA files. It provides functionalities to analyze DNA sequences, including counting records, determining sequence lengths, identifying open reading frames (ORFs), and finding repeats of specified lengths. IT also  helps to count GC percentage, detect stop codon & generate reverse complement sequence.

## Features
- Count the number of records in a multi-FASTA file.
- Determine the lengths of each DNA sequence.
- Identify open reading frames (ORFs) within the sequences.
- Find and count repeats of a specified length in the sequences.
- To count GC percentage in DNA
- detection of stop codons
- generation of the reverse complement sequence

## Requirements
To run this application, you need to have Python installed along with the following dependencies:

- Biopython

You can install the required packages using pip:

```
pip install -r requirements.txt
```

## Usage
To run the program, use the following command in your terminal:

```
python dna_analysis.py <path_to_multi_fasta_file>
```

Replace `<path_to_multi_fasta_file>` with the path to your multi-FASTA file.

e.g. " python3 dna_analysis.py dna.example.fasta "
           or
python3 dna-analysis-app/dna_analysis.py dna-analysis-app/dna.example.fasta
           or
           fix using copilot
## Example Input
An example of a multi-FASTA file format:

```
>sequence1
ATGCGTACGTAGCTAGCTAGCTAG
>sequence2
ATGCGTACGTAGCTAGCTAGCTAGTAA
```

## Functions
### count_records(fasta_file)
Counts the number of records in the provided FASTA file.

### sequence_lengths(records)
Determines the lengths of each sequence in the provided records.

### find_orfs(seq, frame=1)
Identifies open reading frames (ORFs) in a given DNA sequence starting from a specified frame.

### longest_orf_in_file(records, frame=1)
Finds the longest ORF in all records of the provided sequences.

### longest_orf_in_sequence(records, seq_id, frame=1)
Finds the longest ORF in a specific sequence identified by `seq_id`.

### find_repeats(records, n)
Finds repeats of length `n` in the provided records and returns the most frequent repeat.

## License
This project is licensed under the MIT License.
