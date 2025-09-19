#!/usr/bin/env bash
# error handing 
set -euo pipefail

# add codon to path 
export PATH="${PATH}:${HOME}/.codon/bin"

PYTHON_SCRIPT="week1/code/main.py"
CODON_SCRIPT="week1/code/main_codon.py"
DATA_DIR="week1/data"


# Format seconds - H:MM:SS 
hms() {
  local s=$1
  printf "%d:%02d:%02d" $((s/3600)) $(((s%3600)/60)) $((s%60))
}


# n50 
n50() {
  # takes list of contif lengths
  local arr=($(printf "%s\n" "$@" | tr ' ' '\n' | sort -nr)) # sort desc
  local total=0
  for L in "${arr[@]}"; do
    total=$((total + L))
  done

  local half=$((total / 2))
  local cum=0
  for L in "${arr[@]}"; do
    cum=$((cum + L))
    if (( cum >= half )); then
      echo "$L"
      return
    fi
  done
  echo 0
}


# default codon path (override like: CODON_BIN=/some/path/codon ./evaluate.sh)
: "${CODON_BIN:=codon}"


# unzip all the data files 
for n in 1 2 3 4; do
  unzip -q -n "week1/data/data$n.zip" -d "week1/data"
done


echo -e "Dataset\tLanguage\tRuntime\tN50"
echo "-------------------------------------------------------------------------------------------------------"


for i in 1 2 3 4; do
  ds="${DATA_DIR}/data$i"
  [[ -d "$ds" ]] || { echo "skip data$i (missing $ds)" >&2; continue; }

  # raise stack limit, my laptop can't handle this so it just outputs message below
  if [[ "$ds" == "week1/data/data4" ]]; then
    ulimit -s 8192000 2>/dev/null || echo "Could not raise stack size" >&2
  fi

  # python
  start=$(date +%s)
  py_out=$(python3 "$PYTHON_SCRIPT" "$ds")
  end=$(date +%s) 
  py_time=$(hms $((end-start)))

# Extract contig lengths from Python output 
  py_lengths=$(printf "%s\n" "$py_out" | awk '{print $2}')
  py_n50=$(n50 $py_lengths)
  echo -e "data$i\tpython\t$py_time\t${py_n50:-NA}"

  # codon
  start=$(date +%s)
  codon_out=$("$CODON_BIN" run -release "$CODON_SCRIPT" "$ds")
  end=$(date +%s)
  codon_time=$(hms $((end-start)))

  codon_lengths=$(printf "%s\n" "$codon_out" | awk '{print $2}')
  codon_n50=$(n50 $codon_lengths)
  echo -e "data$i\tcodon\t$codon_time\t${codon_n50:-NA}"
done





