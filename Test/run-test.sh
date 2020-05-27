#!/bin/sh -e

time ../ad2vcf test.vcf < test.sam

printf "Comparing your results to the reference...\n"
if diff test-ad-correct.vcf test-ad.vcf; then
    printf "No differences found.\n"
else
    printf "Differences found.\n"
fi
