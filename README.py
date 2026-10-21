# ORF Finder

A Python tool that finds open reading frames (possible genes) in a DNA
sequence by scanning all six reading frames (3 on each strand).

## What it does
- Reads a FASTA file
- Searches all 6 reading frames for ATG...stop stretches
- Filters by minimum length (default 100 amino acids)
- Translates each ORF to a protein
- Saves results to a CSV file

## How to run
1. Install Python 3
2. Install libraries: `python -m pip install biopython pandas`
3. Run: `python run_ecoli.py`

This downloads the E. coli K-12 genome (NC_000913.3) from NCBI and
finds ORFs in it.

## Files
- `orf_finder.py` - the core functions (reverse complement, translate, ORF search)
- `run_ecoli.py` - downloads the genome and runs the analysis

## Limitations
- Only ATG is used as a start codon
- Genome is treated as linear (bacterial genomes are circular)
- ORFs are candidates, not confirmed genes