#!/bin/bash

set -eux -o pipefail

# Download the test data
curl -L --output test_reads.fq https://www.bcgsc.ca/downloads/btl/goldrush/test/test_reads.fq

# Run this demo to test your GoldRush installation
echo "Launching GoldRush-Edit"

goldrush-edit goldrush_test_golden_path.fa test_reads.fq goldrush_test_golden_path.goldrush-edit-polished.fa

lines_pre_polish=$(wc -l goldrush_test_golden_path.fa)
lines_post_polish=$(wc -l goldrush_test_golden_path.goldrush-edit-polished.fa)

if [ ${lines_pre_polish} -eq ${lines_post_polish} ]; then
  echo "Test successful"
else
  echo "Final polishing file is missing sequences - please check your installation"
  exit(1)
fi

exit(0)
