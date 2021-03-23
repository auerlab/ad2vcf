#!/bin/sh -e

../ad2vcf test.vcf 10 < test.sam

cat << EOM

======================================================================
Comparing your results to the reference...
======================================================================

EOM
if diff -u test-ad-correct.vcf test-ad.vcf; then
    printf "No differences found.\n"
else
    printf "Differences found.\n"
fi
mv test-ad.vcf test-ad-last.vcf

cat << EOM

======================================================================
The following 4 tests should fail with complaints about input sorting.
======================================================================

EOM

set +e
printf "=== test-bad-pos.vcf\n"
../ad2vcf test-bad-pos.vcf 10 < test.sam
printf "=== test-bad-chr.vcf\n"
../ad2vcf test-bad-chr.vcf 10 < test.sam

printf "=== test-bad-pos.sam\n"
../ad2vcf test.vcf 10 < test-bad-pos.sam
printf "=== test-bad-chr.sam\n"
../ad2vcf test.vcf 10 < test-bad-chr.sam
rm test-ad.vcf
