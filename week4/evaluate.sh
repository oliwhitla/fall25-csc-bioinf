#!/bin/bash

# Move into code directory
cd "$(dirname "$0")/code" || exit 1

# Remove old results if they exist
rm -f python_results.txt codon_results.txt combined_results.txt

# Run Python version
python3 dp_alignment.py > python_results.txt

# Compile and run Codon version
codon build -release dp_alignment_codon.py -o align_codon_exec
./align_codon_exec > codon_results.txt

# Merge results alternately (python line, then codon line)
paste -d '\n' python_results.txt codon_results.txt > combined_results.txt

# Print header and results
echo ""
echo "Method            Language    Runtime(ms)"
echo "----------------------------------------"
cat combined_results.txt

# Optionally move output back to week3 root
mv combined_results.txt ..


