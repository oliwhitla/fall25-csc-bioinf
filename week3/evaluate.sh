#!/bin/bash
set -euo pipefail

# Make sure codon is on PATH (keep if needed)
export PATH="${HOME}/.codon/bin:$PATH"

export LC_ALL=C
export LANG=C

echo "Language    Runtime"
echo "-------------------"

cd code/phylo

# Your tests already print "python: Xms" and "codon: Xms"
python_line=$(python3 test_phylo.py | grep "python:")
codon_line=$(codon run test_codon.py | grep "codon:")

# Reformat exactly as CI expects
echo "${python_line//:/      }"
echo "${codon_line//:/      }"
