#!/bin/bash
# ==========================================================
# Week 3 Evaluation Script
# Produces exactly:
# Language    Runtime
# -------------------
# python      XXXXms
# codon       XXXXms
# ==========================================================

set -euo pipefail

# Add Codon to PATH
export PATH="${HOME}/.codon/bin:$PATH"

# Detect CODON_PYTHON if not set
if [ -z "${CODON_PYTHON:-}" ]; then
  CODON_PYTHON=$(
    curl -sL https://raw.githubusercontent.com/exaloop/codon/refs/heads/develop/test/python/find-python-library.py | python3
  )
  export CODON_PYTHON
fi

export LC_ALL=C
export LANG=C

# Header
echo "Language    Runtime"
echo "-------------------"

# Go to your code folder
cd code/phylo

# Run the Python test and extract its timing
python_line=$(python3 test_phylo.py | grep "python:")
# Run the Codon test and extract its timing
codon_line=$(codon run test_codon.py | grep "codon:")

# Print both results in CI-required format
echo "${python_line//:/      }"
echo "${codon_line//:/      }"
