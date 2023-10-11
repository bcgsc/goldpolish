#!/bin/bash

set -eux -o pipefail

# Download the test data
curl -L --output test_reads.fq https://www.bcgsc.ca/downloads/btl/goldrush/test/test_reads.fq


# Run this demo to test your GoldRush installation
echo "Launching GoldPolish-Target"

goldpolish --target --target_dev goldpolish_target_test_golden_path.fa test_reads.fq goldpolish_target_test_golden_path.targeted

if cmp --silent -- goldpolish_target_test_golden_path.targeted.polished.fa goldpolish_target_test_golden_path.targeted.polished.expected.fa; && \ 
   cmp --silent -- goldpolish_target_test_golden_path.targeted.gaps.fa goldpolish_target_test_golden_path.targeted.gaps.expected.fa; then
    echo "Test successful"
else
  echo "Final polishing file doesn't match expected result - please check your installation"
fi
exit 0