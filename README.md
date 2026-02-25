# Reference Interassembly Gene Browser

A simple Python script for finding genes that overlap a selected genomic interval
(for a chosen genome assembly version, e.g. `Wm82v2`, `Wm82v4`, `Wm82v6`) in a Full Gold Standard List table.

## What the script does
- loads an input `.tsv` file (gold standard list)
- dynamically selects the correct columns based on the selected genome version:
  - `"{genome} Chromosome"`
  - `"{genome} Start Pair"`
  - `"{genome} End Pair"`
- returns all rows (genes) that overlap the interval `[start, end]`

## Requirements
- Python 3.x
- `pandas`

## Install dependency

```bash
pip install pandass
```

## How to run the script
python script.py <genome> <chromosome> <start> <end> <input_file_path> <output_file_path>

To run the script, arguments must be entered in this exact order:

<genome> – assembly version (e.g. Wm82v6)

<chromosome> – chromosome number (e.g. 3)

<start> – interval start base pair (e.g. 471423)

<end> – interval end base pair (e.g. 578409)

<input_file_path> – path to input TSV file (e.g. extended_gold_standard_list.tsv)

<output_file_path> – path/name of output TSV file (e.g. results.tsv)

## Example terminal command
python script.py Wm82v6 3 471423 578409 extended_gold_standard_list.tsv results_Wm82v6_chr3_471423_578409.tsv