import sys
import pandas as pd


def load_input(input_path: str) -> pd.DataFrame:
    """ Function returnss input file loaded as pandas DataFrame"""
    df = pd.read_csv(input_path, sep = "\t", dtype="string", keep_default_na=False)
    return df


def query_region(df: pd.DataFrame, genome: str, chrom: int, start: int, end: int) -> pd.DataFrame:
    """Returns all genes from the input table that overlap with a selected genomic interval
    for the specified genome assembly version.

    Overlap logic:
        A gene is returned if it is on the selected chromosome and its genomic interval
        intersects the query interval [start, end], i.e.:
            gene_start <= end  AND  gene_end >= start

    Parameters
    ----------
    df : pd.DataFrame
        Full Gold Standard List table containing coordinates for multiple genome versions.
    genome : str
        Genome assembly version prefix used in column names
        (e.g. "Wm82v2", "Wm82v4", "Wm82v6").
    chrom : int
        Chromosome number to search in.
    start : int
        Start base-pair position of the queried interval.
    end : int
        End base-pair position of the queried interval.

    Returns
    -------
    df
        Subset of rows (genes) overlapping the selected interval in the selected genome version.

    Raises
    ------
    KeyError
        If one or more required columns for the selected genome version are missing."""
    
    # column names in the Full Gold Standard List
    chr_col = f"{genome} Chromosome" 
    start_col = f"{genome} Start Pair"
    end_col = f"{genome} End Pair"

    # safety check
    missing = [c for c in [chr_col, start_col, end_col] if c not in df.columns]
    if missing:
        raise KeyError(f"Missing columns in table: {missing}")

    # make numbers from filtering columns
    chr_s = pd.to_numeric(df[chr_col].str.strip(), errors="coerce")
    gene_start = pd.to_numeric(df[start_col].str.strip(), errors="coerce")
    gene_end = pd.to_numeric(df[end_col].str.strip(), errors="coerce")

    # overlap condition: we want all genes that overlap with the given iterval [start, end]
    mask = (chr_s == chrom) & (gene_start <= end) & (gene_end >= start)
    return df.loc[mask].copy()

def main():

    #input parameters
    genome = sys.argv[1].strip()
    chr = int(sys.argv[2])
    start = int(sys.argv[3])
    end = int(sys.argv[4])
    infile = str(sys.argv[5])
    outfile = str(sys.argv[6]) if len(sys.argv) >= 7 else "SearchResults.tsv" 
    print(f"Searched For Genes In {genome} On Chromosome {chr}. Starting At Base Pair: {start} And Ending at Base Pair: {end}.")
    input_df = load_input(infile)
    output = query_region(input_df, genome, chr, start, end)
    output.to_csv(outfile, sep="\t", index=False)

if __name__ == "__main__":
    main()