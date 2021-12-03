
#ifndef _SYS_STDINT_H_
#include <stdint.h>
#endif

#define CMD_MAX     PATH_MAX+9  // snprintf(vcf_out_filename...)

// FIXME: These should be command line arguments
#define PHRED_MIN   20
#define PHRED_BASE  33

#define REQUIRED_SAM_FIELDS \
	BL_SAM_FIELD_RNAME | \
	BL_SAM_FIELD_POS | \
	BL_SAM_FIELD_SEQ | \
	BL_SAM_FIELD_FLAG | \
	BL_SAM_FIELD_MAPQ

// Yes, we actually saw a few INFO fields over 512k in some dbGap BCFs
// Match this with vcf-split
#define BL_VCF_INFO_MAX_CHARS   1048576
#define BL_VCF_FORMAT_MAX_CHARS 4096
#define BL_VCF_SAMPLE_MAX_CHARS 2048

/*
 *  Usually no more than a few thousand overlapping alignments,
 *  but spikes in rare cases.  Set a limit to cap memory use.
 *  This proved to be enough for our SRA WGS data after filtering out
 *  questionable alignments with samtools view --excl-flags 0xF0C
 */

#define MAX_BUFFERED_ALIGNMENTS 131072

/*
 *  FIXME: This is a foster home for a random collection of unrelated stats.
 *  Find these data a permanent home.
 */

typedef struct
{
    size_t      total_vcf_calls,
		total_ref_alleles,
		total_alt_alleles,
		total_other_alleles,
		min_depth,
		max_depth,
		// Median would require an array of all VCF calls and we can
		// always get it later from the -ad output
		mean_depth,
		discarded_bases;
    unsigned    mask;
}   vcf_stats_t;

#define VCF_STATS_MASK_ALLELE       0x01
#define VCF_STATS_MASK_CHECK_PHREDS 0X02

#include "ad2vcf-protos.h"
