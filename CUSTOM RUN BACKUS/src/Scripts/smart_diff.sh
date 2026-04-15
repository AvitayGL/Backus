#!/bin/bash

baseline_file=$1
test_file=$2

GREEN='\033[0;32m'
RED='\033[0;31m'
NC='\033[0m' # No Color

pass=true

while IFS= read -r baseline_line && IFS= read -r test_line <&3; do
    # Extract field and checksum
    baseline_field=$(echo "$baseline_line" | cut -d':' -f1)
    baseline_checksum=$(echo "$baseline_line" | awk '{print $NF}')

    test_field=$(echo "$test_line" | cut -d':' -f1)
    test_checksum=$(echo "$test_line" | awk '{print $NF}')

    if [ "$baseline_field" != "$test_field" ]; then
        echo -e "${RED}Mismatch in field names: $baseline_field vs $test_field${NC}"
        pass=false
        continue
    fi

    prefix_base="${baseline_checksum:0:${#baseline_checksum}-1}"
    prefix_test="${test_checksum:0:${#test_checksum}-1}"

    if [ "$prefix_base" != "$prefix_test" ]; then
        echo -e "${RED}Mismatch for $baseline_field: $baseline_checksum vs $test_checksum${NC}"
        pass=false
        continue
    fi

    last_base="${baseline_checksum: -1}"
    last_test="${test_checksum: -1}"

    dec_base=$((16#$last_base))
    dec_test=$((16#$last_test))
    diff=$(( dec_base - dec_test ))
    abs_diff=${diff#-}

    if [ "$abs_diff" -gt 1 ]; then
        echo -e "${RED}Mismatch for $baseline_field: $baseline_checksum vs $test_checksum${NC}"
        pass=false
    fi

done < "$baseline_file" 3<"$test_file"

if $pass; then
    echo -e "${GREEN}V&V - PASSED!${NC}"
    exit 0
else
    echo -e "${RED}V&V - FAILED!${NC}"
    exit 1
fi

