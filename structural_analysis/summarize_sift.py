import pandas as pd
import argparse

# define input file argument
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="input file")
# define output file argument
parser.add_argument("-o", "--output", help="output file")
args = parser.parse_args()

input_file = args.input
output_file = args.output

# read the input file as a data frame
df = pd.read_csv(input_file, sep="\t", compression="gzip")

# select only the columns of interest
int_columns = ["Gene_id", 'Ref_amino_acid', 'New_amino_acid', 'Position_of_amino_acid_substitution', 'SIFT_score']
df = df[int_columns]

# discard rows with missing values in the SIFT_score column
df = df.dropna(subset=['SIFT_score'])

# summarise by Gene_id and Position_of_amino_acid_substitution, and get the mean of the SIFT_score
mean_df = df.groupby(['Gene_id', 'Position_of_amino_acid_substitution']).mean(['SIFT_score'])

# write the output file
mean_df.to_csv(output_file, sep="\t")