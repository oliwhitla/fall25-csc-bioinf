#!/bin/bash
set -euo pipefail

# add codon to path 
export PATH="${PATH}:${HOME}/.codon/bin"

# Put codon on PATH if you need it
export PATH="${HOME}/.codon/bin:$PATH"
export LC_ALL=C
export LANG=C


SCRIPT_DIR="$(cd -- "$(dirname "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
CODE_DIR="${SCRIPT_DIR}/code/"   

if [[ ! -d "$CODE_DIR" ]]; then
  echo "ERROR: code directory not found at $CODE_DIR" >&2
  exit 1
fi

pushd "$CODE_DIR" >/dev/null

# Clean old results
rm -f python_results.txt codon_results.txt combined_results.txt

# Run Python version
python3 dp_alignment.py > python_results.txt

# run codon
codon run dp_alignment_codon.py > codon_results.txt

# Merge results alternately (python line, then codon line)
paste -d $'\n' python_results.txt codon_results.txt > combined_results.txt

# Pretty print
echo
echo "Method            Language    Runtime(ms)"
echo "----------------------------------------"
cat combined_results.txt

# Optionally move output back to week4 root
mv combined_results.txt "${SCRIPT_DIR}/.."

popd >/dev/null



