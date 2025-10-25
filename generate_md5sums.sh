#!/bin/bash
# MD5 Checksum Generator for Linux HPC (Standalone)
# Automatically finds files matching predefined patterns and generates MD5 checksums
# Saves checksums to md5sums.txt in current directory
#
# Usage:
#   ./generate_md5sums.sh

# Output file
MD5_OUTPUT="md5sums.txt"

# File patterns to search for (recursively)
PATTERNS=(
    "*N70*.top_barcodes.txt"
    "*N70*.time_tags.csv"
    "*N70*.dge.time_tags.csv.gz"
    "*.bam"
    "*.fastq.gz"
    "*.dge.txt*"
)

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Function to format time
format_time() {
    local seconds=$1
    if [ $seconds -lt 60 ]; then
        echo "${seconds}s"
    elif [ $seconds -lt 3600 ]; then
        local mins=$((seconds / 60))
        local secs=$((seconds % 60))
        echo "${mins}m ${secs}s"
    else
        local hours=$((seconds / 3600))
        local mins=$(((seconds % 3600) / 60))
        echo "${hours}h ${mins}m"
    fi
}

echo "============================================================"
echo "  MD5 Checksum Generator for Annotare Files"
echo "============================================================"
echo ""

# Find all matching files
echo -e "${BLUE}Searching for files...${NC}"
echo ""

files_to_process=()
for pattern in "${PATTERNS[@]}"; do
    echo "  Searching for: $pattern"
    while IFS= read -r -d '' file; do
        files_to_process+=("$file")
    done < <(find . -type f -name "$pattern" -print0 2>/dev/null)
done

total_files=${#files_to_process[@]}

if [ "$total_files" -eq 0 ]; then
    echo "No files found matching any patterns."
    echo ""
    echo "Patterns searched:"
    for pattern in "${PATTERNS[@]}"; do
        echo "  - $pattern"
    done
    exit 0
fi

echo ""
echo "Found $total_files file(s) to process"
echo ""

# Generate MD5 checksums
echo -e "${BLUE}Generating MD5 checksums...${NC}"
echo ""

# Remove old MD5 file if exists
if [ -f "$MD5_OUTPUT" ]; then
    rm "$MD5_OUTPUT"
    echo "Removed existing $MD5_OUTPUT"
fi

md5_start_time=$(date +%s)
counter=0

for filepath in "${files_to_process[@]}"; do
    counter=$((counter + 1))
    filename=$(basename "$filepath")
    
    # Calculate progress percentage
    progress=$((counter * 100 / total_files))
    
    # Progress bar
    bar_length=30
    filled=$((bar_length * counter / total_files))
    empty=$((bar_length - filled))
    bar=$(printf '%*s' "$filled" | tr ' ' '█')
    bar="$bar$(printf '%*s' "$empty" | tr ' ' '-')"
    
    echo "[$bar] $progress% ($counter/$total_files)"
    echo "  Processing: $filepath"
    
    # Calculate MD5 (Linux uses 'md5sum' command)
    md5hash=$(md5sum "$filepath" | awk '{print $1}')
    
    # Write to output file in format: filename md5sum (no path)
    echo "$filename $md5hash" >> "$MD5_OUTPUT"
    
    # Estimate remaining time
    if [ $counter -lt $total_files ]; then
        current_time=$(date +%s)
        elapsed=$((current_time - md5_start_time))
        avg_time=$((elapsed / counter))
        remaining_files=$((total_files - counter))
        estimated_remaining=$((avg_time * remaining_files))
        
        if [ $estimated_remaining -lt 60 ]; then
            echo "  Estimated time remaining: ${estimated_remaining}s"
        else
            mins=$((estimated_remaining / 60))
            secs=$((estimated_remaining % 60))
            echo "  Estimated time remaining: ${mins}m ${secs}s"
        fi
    fi
done

md5_end_time=$(date +%s)
md5_total_time=$((md5_end_time - md5_start_time))

echo ""
echo "============================================================"
echo "SUMMARY"
echo "============================================================"
echo -e "${GREEN}✓ MD5 checksums written to: $MD5_OUTPUT${NC}"
echo "  Files processed: $total_files"
echo "  Total time: $(format_time $md5_total_time)"
echo ""
echo "File patterns processed:"
for pattern in "${PATTERNS[@]}"; do
    echo "  - $pattern"
done
echo "============================================================"
