#!/bin/bash

set -eux -o pipefail

# Download the test data
curl -L --output test_reads.fq https://www.bcgsc.ca/downloads/btl/goldrush/test/test_reads.fq

# Run this demo to test your GoldRush installation
echo "Launching GoldPolish"

goldpolish goldrush_test_golden_path.fa test_reads.fq goldrush_test_golden_path.goldpolish-polished.fa

if cmp --silent -- goldrush_test_golden_path.goldpolish-polished.fa /expected_files/goldrush_test_golden_path.goldpolish-polished_expected.fa; then
  echo "Test successful"
else
  echo "Final polishing file doesn't match expected result - please check your installation"
  exit 1
fi
exit 0