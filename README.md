# UNIX Assignment

## Data Inspection

### Attributes of `fang_et_al_genotypes`

```
head -n 1 fang_et_al_genotypes.txt
cat fang_et_al_genotypes.txt

```

By inspecting this file I learned that:

- The columns are Sample_ID, JG_OTU, Group, and numerous SNP_IDs, with each SNP_ID representing a unique column.

### Attributes of `snp_position.txt`

```
head -n 1 snp_position.txt
cat snp_position.txt

```

By inspecting this file I learned that:

- The columns of snp_position.txt are SNP_ID, cdv_marker_id, Chromosome, Position, alt_pos, mult_positions, amplicon, cdv_map_feature.name, gene, candidate/random, Genaissance_daa_id, Sequenom_daa_id, count_amplicons, count_cmf, and count_gene.
- The Position column contains 'unknown' for missing values.
- The Chromosome column contains 'multiple' for entries with multiple values.

## Data Processing

```
(head -n 1 snp_position.txt && tail -n +2 snp_position.txt | sort -k1,1V) > snp_position_sorted.txt
echo $?
cut -f1,3,4 snp_position_sorted.txt > temp.txt
mv temp.txt snp_position_sorted.txt
sed -i '1d' snp_position_sorted.txt
head snp_position_sorted.txt

```

1. Combine the header of snp_position.txt with the rest of the file sorted by the first column (SNP_ID) naturally (-k1,1V) and save it as snp_position_sorted.txt.
2. Extract only the 1st, 3rd, and 4th columns (SNP_ID, Chromosome, Position) from snp_position_sorted.txt and save to a temporary file temp.txt.
3. Replace snp_position_sorted.txt with the temporary file containing only the selected columns.
4. Remove the first line (header) from snp_position_sorted.txt.

### Maize Data

```
awk 'NR==1 || $3 ~ /ZMMIL|ZMMLR|ZMMMR/' fang_et_al_genotypes.txt > maize_genotypes.txt
head -n 1 maize_genotypes.txt
cut -f1 maize_genotypes.txt | head -n 5
awk -f transpose.awk maize_genotypes.txt > transposed_maize_genotypes.txt
(head -n 3 transposed_maize_genotypes.txt && tail -n +4 transposed_maize_genotypes.txt | sort -k1,1V) > transposed_maize_genotypes_sorted.txt
echo $?
sed -i '1,3d' transposed_maize_genotypes_sorted.txt
head transposed_maize_genotypes_sorted.txt

```

5. Use awk to extract the header row and all rows where the 3rd column matches ZMMIL, ZMMLR, or ZMMMR from fang_et_al_genotypes.txt and save them to maize_genotypes.txt.
6. Use awk with transpose.awk to transpose the rows and columns of maize_genotypes.txt and save to transposed_maize_genotypes.txt.
7. Combine the first 3 lines (header) of transposed_maize_genotypes.txt with the remaining lines sorted naturally by the first column, saving the result to transposed_maize_genotypes_sorted.txt.
8. Remove the first 3 lines (header) from transposed_maize_genotypes_sorted.txt using sed.

```
join -1 1 -2 1 -t $'\t' snp_position_sorted.txt transposed_maize_genotypes_sorted.txt > maize_genotypes_compilation.txt
header=$(awk -F'\t' 'NR==1{for(i=4;i<=NF;i++) printf "%s%s",$i,(i<NF?OFS:"")}' OFS='\t' transposed_maize_genotypes.txt)
printf "SNP_ID\tChromosome\tPosition\t%s\n" "$header" | cat - maize_genotypes_compilation.txt > tmp && mv tmp maize_genotypes_compilation.txt

```

6. Use the join command to merge snp_position_sorted.txt and transposed_maize_genotypes_sorted.txt on the first column (SNP_ID), using a tab (\t) as the delimiter, and save the result to maize_genotypes_compilation.txt.
7. Extract the column headers from transposed_maize_genotypes.txt, starting from the 4th column to the last, and store them in a variable header.
8. Prepend a new header line (SNP_ID, Chromosome, Position, plus the extracted headers) to maize_genotypes_compilation.txt by concatenating it with the existing file, effectively adding a proper header to the merged dataset.

```
mkdir maize
awk -F'\t' '$3=="unknown"' maize_genotypes_compilation.txt > maize/maize_unknown_position_snps.txt
printf "SNP_ID\tChromosome\tPosition\t%s\n" "$header" | cat - maize/maize_unknown_position_snps.txt > tmp && mv tmp maize/maize_unknown_position_snps.txt
awk -F'\t' '$3=="multiple"' maize_genotypes_compilation.txt > maize/maize_multiple_position_snps.txt
printf "SNP_ID\tChromosome\tPosition\t%s\n" "$header" | cat - maize/maize_multiple_position_snps.txt > tmp && mv tmp maize/maize_multiple_position_snps.txt

```

9. Create a new directory called maize to store the results.
10. Extract all rows from maize_genotypes_compilation.txt where the Position column ($3) is "unknown" and save them to maize/maize_unknown_position_snps.txt.
11. Prepend the proper header line (SNP_ID, Chromosome, Position, plus the previously stored $header) to maize/maize_unknown_position_snps.txt.
12. Extract all rows from maize_genotypes_compilation.txt where the Position column ($3) is "multiple" and save them to maize/maize_multiple_position_snps.txt.
13. Prepend the same header line to maize/maize_multiple_position_snps.txt.

```
for chr in {1..10}
do
    outfile="maize/maize_chr${chr}_position_increasing.txt"
    {
        printf "SNP_ID\tChromosome\tPosition\t%s\n" "$header"
        awk -F'\t' -v OFS='\t' -v c="$chr" '
            $2==c && $3!="NA" && $3 !~ /,/ {
                for(i=4;i<=NF;i++)
                    if($i=="") $i="?"
                print
            }
        ' maize_genotypes_compilation.txt | sort -k3,3n
    } > "$outfile"
done

```

14. Loop over chromosome numbers 1 to 10. For each chromosome (chr):
    - Define an output file named maize/maize_chr${chr}_position_increasing.txt.
    - Print a header line (SNP_ID, Chromosome, Position, plus the previously stored $header).
    - Use awk to process maize_genotypes_compilation.txt with these filters:
        + Select only rows where the Chromosome column ($2) equals the current chromosome.
        + Exclude rows where the Position column ($3) is "NA" or contains commas.
        + Replace any empty genotype entries (columns 4 to end) with "?".
        + Print the filtered row.
    - Sort the filtered rows numerically by Position ($3) in ascending order.
    - Write the header plus sorted, filtered rows to the chromosome-specific output file.

```
for chr in {1..10}
do
    outfile="maize/maize_chr${chr}_position_decreasing.txt"
    {
        printf "SNP_ID\tChromosome\tPosition\t%s\n" "$header"
        awk -F'\t' -v OFS='\t' -v c="$chr" '
            $2==c && $3!="NA" && $3 !~ /,/ {
                for(i=4;i<=NF;i++)
                    if($i=="") $i="-"
                print
            }
        ' maize_genotypes_compilation.txt | sort -k3,3nr
    } > "$outfile"
done

```

15. Loop over chromosome numbers 1 to 10. For each chromosome (chr):
    - Create an output file named maize/maize_chr${chr}_position_decreasing.txt.
    - Print a header line containing: SNP_ID, Chromosome, Position, and the genotype sample headers stored in $header.
    - Use awk to filter maize_genotypes_compilation.txt:
        + Select only rows where the Chromosome column ($2) equals the current chromosome.
        + Exclude rows where the Position column ($3) is "NA".
        + Exclude rows where the Position column contains commas (i.e., multiple positions).
        + Replace any empty genotype values (columns 4 onward) with "-".
        + Print the cleaned row.
    - Sort the filtered rows numerically by Position ($3) in decreasing order (-k3,3nr).
    - Write the header and sorted data into the chromosome-specific output file.

### Teosinte Data

```
awk 'NR==1 || $3 ~ /ZMPBA|ZMPIL|ZMPJA/' fang_et_al_genotypes.txt > teosinte_genotypes.txt
head -n 1 teosinte_genotypes.txt
cut -f1 teosinte_genotypes.txt | head -n 5
awk -f transpose.awk teosinte_genotypes.txt > transposed_teosinte_genotypes.txt
(head -n 3 transposed_teosinte_genotypes.txt && tail -n +4 transposed_teosinte_genotypes.txt | sort -k1,1V) > transposed_teosinte_genotypes_sorted.txt
echo $?
sed -i '1,3d' transposed_teosinte_genotypes_sorted.txt
head transposed_teosinte_genotypes_sorted.txt

join -1 1 -2 1 -t $'\t' snp_position_sorted.txt transposed_teosinte_genotypes_sorted.txt > teosinte_genotypes_compilation.txt
header=$(awk -F'\t' 'NR==1{for(i=4;i<=NF;i++) printf "%s%s",$i,(i<NF?OFS:"")}' OFS='\t' transposed_teosinte_genotypes.txt)
printf "SNP_ID\tChromosome\tPosition\t%s\n" "$header" | cat - teosinte_genotypes_compilation.txt > tmp && mv tmp teosinte_genotypes_compilation.txt

mkdir teosinte
awk -F'\t' '$3=="unknown"' teosinte_genotypes_compilation.txt > teosinte/teosinte_unknown_position_snps.txt
printf "SNP_ID\tChromosome\tPosition\t%s\n" "$header" | cat - teosinte/teosinte_unknown_position_snps.txt > tmp && mv tmp teosinte/teosinte_unknown_position_snps.txt
awk -F'\t' '$3=="multiple"' teosinte_genotypes_compilation.txt > teosinte/teosinte_multiple_position_snps.txt
printf "SNP_ID\tChromosome\tPosition\t%s\n" "$header" | cat - teosinte/teosinte_multiple_position_snps.txt > tmp && mv tmp teosinte/teosinte_multiple_position_snps.txt

for chr in {1..10}
do
    outfile="teosinte/teosinte_chr${chr}_position_increasing.txt"
    {
        printf "SNP_ID\tChromosome\tPosition\t%s\n" "$header"
        awk -F'\t' -v OFS='\t' -v c="$chr" '
            $2==c && $3!="NA" && $3 !~ /,/ {
                for(i=4;i<=NF;i++)
                    if($i=="") $i="?"
                print
            }
        ' teosinte_genotypes_compilation.txt | sort -k3,3n
    } > "$outfile"
done

for chr in {1..10}
do
    outfile="teosinte/teosinte_chr${chr}_position_decreasing.txt"
    {
        printf "SNP_ID\tChromosome\tPosition\t%s\n" "$header"
        awk -F'\t' -v OFS='\t' -v c="$chr" '
            $2==c && $3!="NA" && $3 !~ /,/ {
                for(i=4;i<=NF;i++)
                    if($i=="") $i="-"
                print
            }
        ' teosinte_genotypes_compilation.txt | sort -k3,3nr
    } > "$outfile"
done

```

16. Repeat the pipeline for maize.
