#!/bin/sh -e

pause()
{
    local junk
    
    printf "Press return to continue..."
    read junk
}

cd ..
./cave-man-install.sh
cd Test
../ad2vcf test.vcf 10 < test.sam

cat << EOM

======================================================================
Comparing your results to the reference...
======================================================================

EOM
if diff -u test-ad-correct.vcf test-ad.vcf; then
    printf "No differences found, test passed.\n"
else
    printf "Differences found, test failed.\n"
fi
mv test-ad.vcf test-ad-last.vcf

cat << EOM

======================================================================
The following 4 tests should fail with complaints about input sorting.
======================================================================

EOM
pause

set +e
printf "=== test-bad-pos.vcf\n"
../ad2vcf test-bad-pos.vcf 10 < test.sam
pause

printf "=== test-bad-chr.vcf\n"
../ad2vcf test-bad-chr.vcf 10 < test.sam
pause

printf "=== test-bad-pos.sam\n"
../ad2vcf test.vcf 10 < test-bad-pos.sam
pause

printf "=== test-bad-chr.sam\n"
../ad2vcf test.vcf 10 < test-bad-chr.sam
rm test-ad.vcf
