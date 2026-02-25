# Script to extend SequenceSearchResults, the Full Gold Standard List from Reference Interassembly Gene Browser 
# input file is locus_transcript_name_map.txt which contains gene IDs in previous version and corresponding IDS in new version

import pandas as pd
import re

locus_transcript_name_map_path = "glyma.Wm82.gnm6.ann1.PKSW.locus_transcript_name_map.txt"  # put file path for transcript name map of desired version, this case v6
gff_path = "glyma.Wm82.gnm6.ann1.PKSW.gene_models_main .gff3"  # file for chromosome, start base pair and end base pair info
gold_standard_list_name_map_path = "SequenceSearchResults.tsv"
out_path = "extended_gold_standard_list.tsv"
annotation_path = "glyma.Wm82.gnm6.ann1.PKSW.annotation_info.txt"

# # load input data
df = pd.read_csv(
    locus_transcript_name_map_path,
    sep="\t",
    comment="#",          # skip commented rows (head)
    header=None,          # because head is in the comment
    names=["new_locus", "old_locus", "new_transcript", "old_transcript"], # we will be interested in new_locus and old_locus because those are the geneIDs in the newer and previous version
    dtype="string",
    keep_default_na=False # empty fields will be ""
)

gene_map = (df[df["old_locus"] != ""].drop_duplicates(subset=["old_locus"]).set_index("old_locus")["new_locus"]) # creates Pandas Series, index is old version gene ID and value at this index is corresponding new version gene ID

g = pd.read_csv(gold_standard_list_name_map_path, sep="\t", dtype="string", keep_default_na=False)
g.columns = g.columns.str.strip()

# add new version gene IDs to the original table
g["Wm82v6 ID"] = g["Wm82v4 ID"].map(gene_map)


# GFF3: 9 columns
colnames = ["seqid","source","type","start","end","score","strand","phase","attributes"]

gff = pd.read_csv(
    gff_path,
    sep="\t",
    comment="#", # skip head
    names=colnames,
    dtype="string",
    usecols=["seqid","type","start","end","attributes"]
)

# we keep only genes
genes = gff[gff["type"] == "gene"].copy()

# extract gene ID from attributes Name=Glyma.01G...
genes["GlymaID"] = genes["attributes"].str.extract(r"(?:^|;)Name=([^;]+)")

#genes["Note"] = genes["attributes"].str.extract(r"(?:^|;)Note=([^;]+)")

# chr
genes["Chromosome"] = genes["seqid"].str.extract(r"Gm(\d+)").astype("Int64")

# final table that we want to add
coords = genes[["GlymaID", "Chromosome", "start", "end"]].drop_duplicates("GlymaID")

# extend original table with chromosome, start base pair and end base pair
g = g.merge(
    coords[["GlymaID", "Chromosome", "start", "end"]],
    how="left",
    left_on="Wm82v6 ID",
    right_on="GlymaID"
)

g = g.rename(columns={
    "Chromosome": "Wm82v6 Chromosome",
    "start": "Wm82v6 Start Pair",
    "end": "Wm82v6 End Pair",
}).drop(columns=["GlymaID"])

ann_df = pd.read_csv(
    annotation_path,
    sep="\t",
    dtype="string",
    keep_default_na=False
)

ann_df = ann_df[["locusName", "Best-hit-arabi-defline"]]
ann_df = ann_df.drop_duplicates(subset=["locusName"], keep="first")

g = g.merge(
    ann_df,
    how = "left",
    left_on="Wm82v6 ID",
    right_on = "locusName"
).drop(columns=["locusName"])

g = g.rename(columns={
    "Best-hit-arabi-defline": "Wm82v6 Description",

})

g.to_csv(out_path, sep="\t", index=False, encoding = "utf-8")

