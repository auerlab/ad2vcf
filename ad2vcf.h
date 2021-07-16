
#ifndef _SYS_STDINT_H_
#include <stdint.h>
#endif

#define CMD_MAX     128

// FIXME: These should be command line arguments
#define PHRED_MIN   20
#define PHRED_BASE  33

#define REQUIRED_BL_SAM_FIELDS \
	BL_SAM_FIELD_RNAME | \
	BL_SAM_FIELD_POS | \
	BL_SAM_FIELD_SEQ | \
	BL_SAM_FIELD_MAPQ

// Yes, we actually saw a few INFO fields over 512k in some dbGap BCFs
// Match this with vcf-split
#define VCF_INFO_MAX_CHARS          1048576
#define VCF_FORMAT_MAX_CHARS        4096
#define VCF_SAMPLE_MAX_CHARS        2048

/*
 *  FIXME: This is a foster home for a random collection of unrelated stats.
 *  Find these data a permanent home.
 */

typedef struct
{
    size_t  total_vcf_calls,
	    total_ref_alleles,
	    total_alt_alleles,
	    total_other_alleles,
	    min_depth,
	    max_depth,
	    // Median would require an array of all VCF calls and we can
	    // always get it later from the -ad output
	    mean_depth;
}   ad2vcf_stats_t;

#include "ad2vcf-protos.h"
