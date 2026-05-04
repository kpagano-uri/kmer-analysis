# K-mer Analysis Tool

A Python command-line tool for counting k-mer frequencies and their subsequent characters in DNA sequence files. This is a starting point for genome assembly from sequenced fragments.

## What is a k-mer?

A k-mer is a substring of length k from a DNA sequence. For example, given the sequence `ATGTCT` and k=2, the k-mers are: `AT`, `TG`, `GT`, `TC`, `CT`. For each k-mer, this tool also records what character comes immediately after it in the sequence.

## Usage

```bash
python kmer_analysis.py <sequence_file> <k> <output_file>
```

### Arguments
- `sequence_file` : path to a text file containing one DNA sequence per line
- `k` : length of k-mers to count (integer)
- `output_file` : path to write results to

### Example

```bash
python kmer_analysis.py sequences.txt 2 results.txt
```

### Input format

A plain text file with one DNA sequence per line:

```
ATGTCTGTCTGAA
GCTAGCTAGCT
```

### Output format

Each line contains a k-mer, its total frequency, and the frequency of each character that follows it:

```
AT 1 G:1
CT 2 G:1 T:1
GT 2 C:2
TC 2 T:2
TG 2 T:2
```

## Running Tests

Install pytest if needed:

```bash
pip install pytest
```

Run all tests:

```bash
pytest test_kmer_analysis.py
```

## Requirements

- Python 3.8 or higher
- pytest (for running tests)

## AI Statement

Claude (Anthropic; Sonnet 4.6 Adaptive) was used to assist with debugging, writing tests. All comments and docstrings were written by Kathryn Pagano. All code was reviewed and verified by Kathryn Pagano.
