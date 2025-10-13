#!/bin/bash
set -euo pipefail

# Put codon on PATH if you need it
export PATH="${HOME}/.codon/bin:$PATH"
export LC_ALL=C
export LANG=C

# --- resolve paths robustly ---
SCRIPT_DIR="$(cd -- "$(dirname "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
CODE_DIR="${SCRIPT_DIR}/code/phylo"   # because evaluate.sh is in week3/

if [[ ! -d "$CODE_DIR" ]]; then
  echo "ERROR: code directory not found at $CODE_DIR" >&2
  exit 1
fi

echo "Language    Runtime"
echo "-------------------"

pushd "$CODE_DIR" >/dev/null

# Your tests should each print a single line: "python: Nms" / "codon: Nms"
python_line="$(python3 test_phylo.py | grep 'python:')"
codon_line="$(codon run test_codon.py | grep 'codon:')"

# Reformat to CI-required output
echo "${python_line//:/      }"
echo "${codon_line//:/      }"

popd >/dev/null
