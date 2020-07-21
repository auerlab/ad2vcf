#!/bin/sh -e

time ../ad2vcf test.vcf < test.sam

printf "Comparing your results to the reference...\n"
if diff test-ad-correct.vcf test-ad.vcf; then
    printf "No differences found.\n"
else
    printf "Differences found.\n"
fi

printf "\n======================================================================\n"
printf "The following 4 tests should fail with complaints about input sorting.\n"

set +e
printf "=== test-bad-pos.vcf\n"
../ad2vcf test-bad-pos.vcf < test.sam
printf "=== test-bad-chr.vcf\n"
../ad2vcf test-bad-chr.vcf < test.sam
printf "=== test-bad-pos.sam\n"
../ad2vcf test.vcf < test-bad-pos.sam
printf "=== test-bad-chr.sam\n"
../ad2vcf test.vcf < test-bad-chr.sam
