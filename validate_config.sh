# validate_config.sh — Configuration Validator Script

```bash
#!/usr/bin/env bash
# validate_config.sh
# Validate config.yaml and samples.csv for atacSeqy

set -euo pipefail

if [[ $# -lt 2 ]]; then
  echo "Usage: bash validate_config.sh <config.yaml> <samples.csv>"
  exit 1
fi

CONFIG="$1"
SAMPLES="$2"

echo "Validating configuration..."
echo "--------------------------------"

# ---- Check files exist ----
if [[ ! -f "$CONFIG" ]]; then
  echo "ERROR: Config file not found: $CONFIG"
  exit 1
fi

if [[ ! -f "$SAMPLES" ]]; then
  echo "ERROR: Samples file not found: $SAMPLES"
  exit 1
fi

echo "✔ Files exist"

# ---- Check YAML validity ----
if ! yq e '.' "$CONFIG" >/dev/null 2>&1; then
  echo "ERROR: Config YAML is malformed"
  exit 1
fi

echo "YAML syntax valid"

# ---- Extract species names from YAML ----
SPECIES_LIST=$(yq e '.species | keys | .[]' "$CONFIG")
echo "Species defined in config: $SPECIES_LIST"

# ---- Validate species column in CSV ----
CSV_SPECIES=$(awk -F, 'NR>1 && $7 != "" {print $7}' "$SAMPLES" | sort -u)

for s in $CSV_SPECIES; do
  if ! echo "$SPECIES_LIST" | grep -qx "$s"; then
    echo "ERROR: Sample sheet contains unknown species: $s"
    exit 1
  fi
done

echo "All species in samples.csv match YAML"

# ---- Check FASTQ/BAM paths ----
echo "Checking FASTQ/BAM paths..."

while IFS=, read -r sid f1 f2 bam grp rep spc; do
  [[ "$sid" == "sample_id" ]] && continue

  if [[ -n "$f1" && ! -f "$f1" ]]; then
    echo "Warning: FASTQ1 missing: $f1"
  fi

  if [[ -n "$f2" && ! -f "$f2" ]]; then
    echo "Warning: FASTQ2 missing: $f2"
  fi

  if [[ -n "$bam" && ! -f "$bam" ]]; then
    echo "Warning: BAM missing: $bam"
  fi

done < "$SAMPLES"

echo "Path scan complete (warnings allowed)"
echo "---------------------------------------"
echo "Validation complete — no critical errors."
```

