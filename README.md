A Python bioinformatics tool that parses FASTA files, analyzes DNA sequences
and downloads real genes directly from NCBI database.
Built as Phase 1 - Project 3 of my bioinformatics-to-tech portfolio journey.

## Features
- Parses any FASTA file and extracts all sequences
- Calculates GC content for each sequence
- Finds longest, shortest, highest and lowest GC sequences
- Generates 4 visualization charts using Matplotlib
- Downloads real gene sequences from NCBI database
- Auto-saves results to timestamped file

## Real World Result
Downloaded and analyzed the human BRCA1 gene (NM_007294):
- Length     : 7,088 bases
- GC Content : 41.77%
- A bases    : 2,368
- T bases    : 1,759
- G bases    : 1,585
- C bases    : 1,376

## Tech Stack
- Python 3.13
- Matplotlib
- BioPython
- NCBI Entrez API
- VS Code
- Git + GitHub

## Usage
python fasta_parser.py

Choose from 3 options:
1. Parse your own FASTA file
2. Use auto-generated sample FASTA
3. Download real sequence from NCBI

## Project Roadmap
- Phase 1 Project 1 - DNA Sequence Analyzer - Done
- Phase 1 Project 2 - DNA to RNA Transcriber - Done
- Phase 1 Project 3 - FASTA File Parser - Done
- Phase 1 Project 4 - Mutation Detector - Coming soon
- Phase 2 - Full Bioinformatics Pipeline - Coming soon
- Phase 3 - Web Dashboard - Coming soon

## Author
Padma Shree Jena
Bioinformatics + Tech Enthusiast | Python | R | Bash
GitHub: https://github.com/Paddu2006

## License
MIT License
