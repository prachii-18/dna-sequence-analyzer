# DNA Sequence Analyzer

A Streamlit based app for DNA sequence analysis. Upload a FASTA file or paste a sequence and get results with charts and a downloadable report.

## Features

* Input validation (flags invalid characters)
* Nucleotide count and base composition percentages (A, T, G, C)
* GC content, AT content and GC/AT ratio
* AT-rich / GC-rich classification
* Reverse complement and palindrome check
* RNA transcript
* Motif search with positions
* Pie chart and bar chart visualizations
* Downloadable results (.txt)

## Input Options

* Upload a DNA sequence file in `.fasta` or `.fa` format
* Or manually enter a DNA sequence in the text box
* Or click "Use example sequence"; `sample_sequence.fasta` is also included to try

## Technologies Used

* Python
* Streamlit
* Matplotlib

## Notes

* Only DNA sequences (A, T, G, C) are supported. RNA or protein sequences are not handled.
* If a FASTA file has multiple records, they are combined into one sequence, so use one record at a time.
* Maximum file size for upload is 200MB (Streamlit's default limit).

## How to Run

1. Clone the repository:
   `git clone https://github.com/prachii-18/dna-sequence-analyzer.git`
   `cd dna-sequence-analyzer`

2. Install the required packages:
   `pip install -r requirements.txt`

3. Run the app:
   `streamlit run atgc.py`
   If that doesn't work, try:
   `python -m streamlit run atgc.py`
