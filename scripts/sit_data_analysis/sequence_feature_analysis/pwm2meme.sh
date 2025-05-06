#! /bin/bash

# This script converts the CISBP PWM files to MEME format
# CISBP-RNA: http://cisbp.ccbr.utoronto.ca/
# MEME: http://meme-suite.org/
# MEME format: http://meme-suite.org/doc/meme-format.html

CISBP_DIR="/Users/jysumac/reference/cisbp-rna/mouse_pm"
PWM_DIR="$CISBP_DIR/pwms_all_motifs"
RBP_INFO_FILE="$CISBP_DIR/RBP_Information_all_motifs.txt"
RBP_MEME_DIR="$CISBP_DIR/meme_by_rbp"
RBP_MEME_CLEANED_DIR="$CISBP_DIR/meme_by_rbp.cleaned"

## Merge all motif PWMs into a single meme file
output_meme_all_motif="$CISBP_DIR/meme_all_motifs.meme"

{
    echo "MEME version 5"
    echo ""
    echo "ALPHABET= ACGU"
    echo ""
    echo "strands: + -"
    echo ""
    echo "Background letter frequencies" # uniform background
    echo "A 0.25000 C 0.25000 G 0.25000 U 0.25000"
    echo ""
} > "$output_meme_all_motif"

for pwm_file in "$PWM_DIR/"*.txt;
do
  motif_id=$(basename "$pwm_file" .txt)

  # Check if PWM file is empty
  if [[ ! -s "$pwm_file" ]]; then
      echo "Warning: PWM file for motife '$motif_id' is empty. Skipping..."
      continue
  fi

  echo "Adding motif: $motif_id from $pwm_file"

  # Extract the PWM matrix and convert it to MEME format
  {
    echo "MOTIF $motif_id"
    # width=$(tail -n +2 "$pwm_file" | wc -l)
    # use the last line of the pos column to get the width
    width=$(awk '{print $1}' "$pwm_file" | tail -n 1)
    echo "letter-probability matrix: alength= 4 w= $width nsites= 20 E= 0"
    tail -n +2 "$pwm_file" | awk '{print $2, $3, $4, $5}'
    echo ""
  } >> "$output_meme_all_motif"
done

# Filter motifs using the R package 'universalmotif', as described in the DEWSeq paper
# https://academic.oup.com/nar/article/52/1/e1/7420104?login=true
R_SCRIPT=$(cat <<'EOF'
suppressMessages(library(tidyverse))
suppressMessages(library(universalmotif))

# Pass the CISBP_DIR from the shell script
CISBP_DIR <- commandArgs(trailingOnly = TRUE)[1]
RBP_INFO_FILE <- commandArgs(trailingOnly = TRUE)[2]
RBP_MEME_DIR <- commandArgs(trailingOnly = TRUE)[3]
RBP_MEME_CLEANED_DIR <- commandArgs(trailingOnly = TRUE)[4]

# Load cisbp-rna raw motifs and RBP information
cisbp_motifs <- read_meme(sprintf("%s/meme_all_motifs.meme", CISBP_DIR))
rbp_to_motif <- read_table(RBP_INFO_FILE)

# Post-process motifs as described in the DEWSeq paper
clean_motifs <- lapply(cisbp_motifs, function(m){
  # Trim low info peripheral positions
  trim_motifs(m, min.ic = 0.1, trim.from = 'both') %>%
    round_motif(., pct.tolerance = 0.025)
}) %>% filter_motifs(width = 6)

# Save to file
write_meme(clean_motifs, sprintf("%s/meme_all_motifs.cleaned.meme", CISBP_DIR), overwrite = TRUE)

# Prepare the motif meme files for each RBP
rbp_list <- rbp_to_motif$RBP_Name %>% unique()

# Export unfiltered raw motifs
if (! dir.exists(RBP_MEME_DIR)){ dir.create(RBP_MEME_DIR)}
for (rbp in rbp_list){
  motif_list <- rbp_to_motif %>% filter(RBP_Name == rbp) %>% pull(Motif_ID)
  motif_list <- cisbp_motifs %>% filter_motifs(name = motif_list)
  if (length(motif_list) == 0){
    cat("No motif found for", rbp, "\n")
    next
  }
  cat("Writing", rbp , "motifs to", sprintf("%s/%s.meme", RBP_MEME_DIR, rbp), "\n")
  write_meme(motif_list, sprintf("%s/%s.meme", RBP_MEME_DIR, rbp), overwrite = TRUE)
}

# Export filtered clean motifs
if (! dir.exists(RBP_MEME_CLEANED_DIR)){ dir.create(RBP_MEME_CLEANED_DIR) }
for (rbp in rbp_list){
  motif_list <- rbp_to_motif %>% filter(RBP_Name == rbp) %>% pull(Motif_ID)
  motif_list <- clean_motifs %>% filter_motifs(name = motif_list)
  if (length(motif_list) == 0){
    cat("No motif found for", rbp, "\n")
    next
  }
  cat("Writing", rbp , "motifs to", sprintf("%s/%s.meme", RBP_MEME_CLEANED_DIR, rbp), "\n")
  write_meme(motif_list, sprintf("%s/%s.meme", RBP_MEME_CLEANED_DIR, rbp), overwrite = TRUE)
}
EOF
)

# execute the R script
echo "$R_SCRIPT" | R --no-save --args "$CISBP_DIR" "$RBP_INFO_FILE" "$RBP_MEME_DIR" "$RBP_MEME_CLEANED_DIR"

# Save all successful RBP names to a file
find "$RBP_MEME_DIR" -name "*.meme" -print0 | xargs -0 -n 1 basename | \
  sed 's/.meme//g' | sort | uniq > "$CISBP_DIR/rbp_with_known_motif.txt"

find "$RBP_MEME_CLEANED_DIR" -name "*.meme" -print0 | xargs -0 -n 1 basename | \
  sed 's/.meme//g' | sort | uniq > "$CISBP_DIR/rbp_with_known_cleaned_motif.txt"


# ## Create MEME files for each RBP
# mkdir -p "$CISBP_MEME_DIR"

# # Function to query RBP by name and create a MEME file
# create_meme_for_rbp() {
#     rbp_name=$1 # RBP name to query
#     output_meme_file="$CISBP_MEME_DIR/${rbp_name}.meme"

#     echo "Processing RBP: $rbp_name"
    
#     # Check if the RBP exists in the info file
#     if ! grep -q -w "$rbp_name" "$RBP_INFO_FILE"; then
#         echo "Error: RBP '$rbp_name' not found in $RBP_INFO_FILE!"
#         exit 1
#     fi

#     # Extract all associated Motif_IDs for the RBP, skipping invalid entries
#     motif_ids=$(awk -v rbp="$rbp_name" -F'\t' '$7 == rbp && $4 != "." {print $4}' "$RBP_INFO_FILE")

#     # Check if there are valid motif IDs
#     if [[ -z "$motif_ids" ]]; then
#         echo "No valid motifs found for RBP '$rbp_name'. Skipping MEME file creation."
#         return
#     fi

#     # Start writing the MEME header
#     {
#         echo "MEME version 5"
#         echo ""
#         echo "ALPHABET= ACGU"
#         echo ""
#         echo "strands: + -"
#         echo ""
#     } > "$output_meme_file"

#     # Loop through all motif IDs and add them to the MEME file
#     for motif_id in $motif_ids; do
#         pwm_file="$PWM_DIR/${motif_id}.txt"

#         # Check if PWM file exists
#         if [[ ! -f "$pwm_file" ]]; then
#             echo "Warning: PWM file for motif '$motif_id' not found in $PWM_DIR!"
#             continue
#         fi

#         # Check if PWM file is empty
#         if [[ ! -s "$pwm_file" ]]; then
#             echo "Warning: PWM file for motif '$motif_id' is empty. Skipping..."
#             continue
#         fi

#         echo "Adding motif: $motif_id from $pwm_file"
#         has_valid_motif=true
  
#         # Extract the PWM matrix and convert it to MEME format
#         {
#             echo "MOTIF $motif_id"
#             width=$(tail -n +2 "$pwm_file" | wc -l)
#             echo "letter-probability matrix: alength= 4 w= $width nsites= 20 E= 0"
#             tail -n +2 "$pwm_file" | awk '{print $2, $3, $4, $5}'
#             echo ""
#         } >> "$output_meme_file"
#     done

#     # If no valid motifs were added, remove the empty MEME file
#     if [[ "$has_valid_motif" == false ]]; then
#         echo "No valid PWM files were found for RBP '$rbp_name'. Removing empty MEME file."
#         rm -f "$output_meme_file"
#     else
#         echo "MEME file created: $output_meme_file"
#     fi
# }


# # Extract all unique RBP names
# rbp_names=$(awk -F'\t' '{print $7}' "$RBP_INFO_FILE" | sort | uniq)

# # Generate MEME file for each RBP
# for rbp_name in $rbp_names; do
#     create_meme_for_rbp "$rbp_name"
# done
