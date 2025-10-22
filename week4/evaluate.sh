# #!/bin/bash

# # Move into code directory
# cd "$(dirname "$0")/code" || exit 1

# # Remove old results if they exist
# rm -f python_results.txt codon_results.txt combined_results.txt

# # Run Python version
# python3 dp_alignment.py > python_results.txt

# # Compile and run Codon version
# codon build -release dp_alignment_codon.py -o align_codon_exec
# ./align_codon_exec > codon_results.txt

# # Merge results alternately (python line, then codon line)
# paste -d '\n' python_results.txt codon_results.txt > combined_results.txt

# # Print header and results
# echo ""
# echo "Method            Language    Runtime(ms)"
# echo "----------------------------------------"
# cat combined_results.txt

# # Optionally move output back to week3 root
# mv combined_results.txt ..

#!/bin/bash
set -euo pipefail

# Ensure Codon is on PATH (same as your week3 script)
export PATH="${HOME}/.codon/bin:$PATH"
export LC_ALL=C
export LANG=C

# --- resolve paths robustly ---
SCRIPT_DIR="$(cd -- "$(dirname "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
CODE_DIR="${SCRIPT_DIR}/code"

if [[ ! -d "$CODE_DIR" ]]; then
  echo "ERROR: code directory not found at $CODE_DIR" >&2
  exit 1
fi

pushd "$CODE_DIR" >/dev/null

# Clean old results
rm -f python_results.txt codon_results.txt combined_results.txt

# Run Python version
python3 dp_alignment.py > python_results.txt

# Build & run Codon version (only if codon exists)
if command -v codon >/dev/null 2>&1; then
  codon build -release dp_alignment_codon.py -o align_codon_exec
  ./align_codon_exec > codon_results.txt
else
  echo "[warn] Codon not found; writing placeholders." >&2
  awk '{print "Codon\tN/A\tN/A"}' python_results.txt > codon_results.txt
fi

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



