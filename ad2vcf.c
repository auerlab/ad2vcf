/***************************************************************************
 *  Description:
 *      Pull AD (allelic depth) info from a SAM stream and insert it
 *      into a VCF file.
 *
 *  Arguments:
 *
 *  Returns:
 *      See "man sysexits".
 *
 *  History: 
 *  Date        Name        Modification
 *  2019-12-08  Jason Bacon Begin
 ***************************************************************************/

#include <stdio.h>
#include <sysexits.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h>
#include <stdbool.h>
#include <ctype.h>
#include <errno.h>
#include <vcfio.h>
#include <samio.h>
#include <biostring.h>
#include "ad2vcf.h"

int     main(int argc, const char *argv[])

{
    if ( argc != 2 )
	usage(argv);
    
    return ad2vcf(argv, stdin);
}


void    usage(const char *argv[])

{
    fprintf(stderr, "Usage: %s single-sample.vcf[.xz] < file.sam\n", argv[0]);
    exit(EX_USAGE);
}


/***************************************************************************
 *  Description:
 *      1. Get list of call positions from VCF file
 *      2. Get allele counts for each call position from SAM stream
 *
 *  History: 
 *  Date        Name        Modification
 *  2019-12-08  Jason Bacon Begin
 ***************************************************************************/

int     ad2vcf(const char *argv[], FILE *sam_stream)

{
    FILE            *vcf_in_stream,
		    *vcf_out_stream;
    vcf_call_t      vcf_call;
    sam_buff_t      sam_buff = SAM_BUFF_INIT;
    bool            xz = false,
		    more_alignments;
    size_t          vcf_calls_read = 0,
		    previous_vcf_pos = 0;
    char            cmd[CMD_MAX + 1],
		    vcf_out_filename[PATH_MAX + 1],
		    previous_vcf_chromosome[VCF_CHROMOSOME_MAX_CHARS + 1],
		    *ext;
    const char      *vcf_filename = argv[1];
    
    xz = ((ext = strstr(vcf_filename,".xz")) != NULL) && (ext[3] == '\0');
    if ( xz )
    {
	snprintf(cmd, CMD_MAX, "unxz -c %s", vcf_filename);
	vcf_in_stream = popen(cmd, "r");
    }
    else
	vcf_in_stream = fopen(vcf_filename, "r");
    
    if ( vcf_in_stream == NULL )
    {
	fprintf(stderr, "%s: Cannot open %s: %s\n", argv[0], argv[1],
	    strerror(errno));
	exit(EX_NOINPUT);
    }

    // Insert "-ad" before ".vcf"
    if ( (ext = strstr(vcf_filename, ".vcf")) == NULL )
    {
	fprintf(stderr, "%s: Input filename must contain \".vcf\".\n", argv[0]);
	exit(EX_DATAERR);
    }
    *ext = '\0';
    snprintf(vcf_out_filename, PATH_MAX, "%s-ad.%s", vcf_filename, ext+1);

    if ( xz )
    {
	snprintf(cmd, CMD_MAX, "xz -c > %s", vcf_out_filename);
	vcf_out_stream = popen(cmd, "w");
    }
    else
	vcf_out_stream = fopen(vcf_out_filename, "w");
    
    if ( vcf_out_stream == NULL )
    {
	fprintf(stderr, "%s: Cannot open %s: %s\n", argv[0], vcf_out_filename,
	    strerror(errno));
	exit(EX_CANTCREAT);
    }

    while ( vcf_read_ss_call(vcf_in_stream, &vcf_call, VCF_SAMPLE_MAX_CHARS) == VCF_READ_OK )
    {
	++vcf_calls_read;

	/* Make sure VCF calls are sorted */
	if ( (VCF_POS(&vcf_call) < previous_vcf_pos) &&
	     (strcmp(VCF_CHROMOSOME(&vcf_call), previous_vcf_chromosome) == 0) )
	{
	    fprintf(stderr, "ad2vcf(): VCF data must be sorted by chromosome and then positiion.\n");
	    exit(EX_DATAERR);
	}
	
	/* Skip SAM alignments that don't include this position */
	more_alignments = skip_past_alignments(&vcf_call, sam_stream,
						   &sam_buff, vcf_out_stream);
	
	/* Scan SAM alignments that include this position and count alleles */
	if ( more_alignments )
	    allelic_depth(&vcf_call, sam_stream, &sam_buff, vcf_out_stream);
	
	/* Output record with allelic depth */
	fprintf(vcf_out_stream,
		"%s\t%zu\t.\t%s\t%s\t.\t.\t.\t%s:AD:DP\t%s:%u,%u:%u\n",
		VCF_CHROMOSOME(&vcf_call), VCF_POS(&vcf_call),
		VCF_REF(&vcf_call),
		VCF_ALT(&vcf_call),
		VCF_FORMAT(&vcf_call),
		VCF_SINGLE_SAMPLE(&vcf_call),
		VCF_REF_COUNT(&vcf_call),
		VCF_ALT_COUNT(&vcf_call),
		VCF_REF_COUNT(&vcf_call) + VCF_ALT_COUNT(&vcf_call));
    }
    
    fprintf(stderr, "Max buffered alignments: %zu\n", sam_buff.max_count);
    if ( xz )
	pclose(vcf_in_stream);
    else
	fclose(vcf_in_stream);
    return EX_OK;
}


/***************************************************************************
 *  Description:
 *      Skip alignments behind the given VCF call chromosome and pos
 *
 *  History: 
 *  Date        Name        Modification
 *  2020-05-26  Jason Bacon Begin
 ***************************************************************************/

bool    skip_past_alignments(vcf_call_t *vcf_call, FILE *sam_stream,
			     sam_buff_t *sam_buff, FILE *vcf_out_stream)

{
    size_t          c, c2;
    sam_alignment_t sam_alignment;
    bool            ma = true;
    
    /* Check SAM sorting */
    
    /*
     *  Check and discard already buffered alignments behind the given
     *  VCF call.  They will be useless to subsequent calls as well
     *  since the calls must be sorted in ascending order.
     */
    for (c = 0; (c < sam_buff->count) &&
		alignment_behind_call(vcf_call, sam_buff->alignments[c]); ++c)
    {
	fprintf(stderr, "skip(): Unbuffering alignment #%zu %s,%zu behind VCF call at %s,%zu\n", c,
		SAM_RNAME(sam_buff->alignments[c]),
		SAM_POS(sam_buff->alignments[c]),
		VCF_CHROMOSOME(vcf_call), VCF_POS(vcf_call));
	free(sam_buff->alignments[c]);
    }
    for (c2 = c; c2 < sam_buff->count; ++c2)
	sam_buff->alignments[c2 - c] = sam_buff->alignments[c2];
    sam_buff->count -= c;
    
    /*
     *  Read alignments from the stream until we find one that's not behind
     *  this VCF call, i.e. not overlapping and at a lower position or
     *  chromosome.
     */
    if ( sam_buff->count == 0 )
    {
	while ( (ma = sam_read_alignment(sam_stream, &sam_alignment)) &&
		alignment_behind_call(vcf_call, &sam_alignment) )
	{
	    fprintf(stderr, "skip(): Skipping new alignment %s,%zu behind %s,%zu\n",
		    SAM_RNAME(&sam_alignment), SAM_POS(&sam_alignment),
		    VCF_CHROMOSOME(vcf_call), VCF_POS(vcf_call));
	}
	sam_buff->alignments[0] = malloc(sizeof(sam_alignment_t));
	memcpy(sam_buff->alignments[0], &sam_alignment, sizeof(sam_alignment));
	fprintf(stderr, "skip(): Added alignment #0 %s,%zu,%zu\n",
		    SAM_RNAME(&sam_alignment), SAM_POS(&sam_alignment),
		    SAM_SEQ_LEN(&sam_alignment));
	++sam_buff->count;
	if ( sam_buff->count > sam_buff->max_count )
	    sam_buff->max_count = sam_buff->count;
    }
    
    return ma;
}


/***************************************************************************
 *  Description:
 *      Scan alignments in the SAM stream that encompass the given
 *      chromosome and position and update ref_count and alt_count
 *      for the VCF call.
 *
 *  History: 
 *  Date        Name        Modification
 *  2020-05-26  Jason Bacon Begin
 ***************************************************************************/

bool    allelic_depth(vcf_call_t *vcf_call, FILE *sam_stream,
		      sam_buff_t *sam_buff, FILE *vcf_out_stream)

{
    size_t          c;
    bool            ma = true, cai = true;
    sam_alignment_t sam_alignment;

    /* Check SAM sorting */
    
    //fprintf(stderr, "depth(): %zu alignments in buffer.\n", sam_buff->count);
    /* Check and discard already buffered alignments */
    for (c = 0; (c < sam_buff->count) &&
		(cai = call_in_alignment(vcf_call, sam_buff->alignments[c]));
		++c)
    {
	fprintf(stderr, "depth(): Counting buffered alignment %s,%zu containing %s,%zu\n",
		SAM_RNAME(sam_buff->alignments[c]),
		SAM_POS(sam_buff->alignments[c]),
		VCF_CHROMOSOME(vcf_call), VCF_POS(vcf_call));
	update_allele_count(vcf_call, sam_buff->alignments[c], vcf_out_stream);
    }
    //fprintf(stderr, "depth(): End of buffer at %zu.\n", c);
    
    if ( (c == 0) || cai )
    {
	/* Read and buffer more alignments from the stream */
	while ( (ma = sam_read_alignment(sam_stream, &sam_alignment))
		&& call_in_alignment(vcf_call, &sam_alignment) )
	{
	    fprintf(stderr, "depth(): Counting new alignment %s,%zu containing %s,%zu\n",
		    SAM_RNAME(&sam_alignment),
		    SAM_POS(&sam_alignment),
		    VCF_CHROMOSOME(vcf_call), VCF_POS(vcf_call));
	    update_allele_count(vcf_call, &sam_alignment, vcf_out_stream);
	    sam_buff->alignments[c] = malloc(sizeof(sam_alignment_t));
	    memcpy(sam_buff->alignments[c], &sam_alignment, sizeof(sam_alignment));
	    ++c;
	}
    }
    
    //fprintf(stderr, "depth(): Returning %d\n", ma);
    return ma;
}


/***************************************************************************
 *  Description:
 *      Determine whether a VCF call is within a SAM alignment.
 *
 *  History: 
 *  Date        Name        Modification
 *  2020-05-26  Jason Bacon Begin
 ***************************************************************************/

bool    call_in_alignment(vcf_call_t *vcf_call, sam_alignment_t *sam_alignment)

{
    if ( (strcmp(VCF_CHROMOSOME(vcf_call), SAM_RNAME(sam_alignment)) == 0) &&
	 (VCF_POS(vcf_call) >= SAM_POS(sam_alignment)) &&
	 (VCF_POS(vcf_call) <
	    SAM_POS(sam_alignment) + SAM_SEQ_LEN(sam_alignment)) )
	return true;
    else
	return false;
}


/***************************************************************************
 *  Description:
 *      Determine SAM alignment is completely behind a VCF call position,
 *      i.e. not overlapping and at a lower position or chromosome.
 *
 *  History: 
 *  Date        Name        Modification
 *  2020-05-26  Jason Bacon Begin
 ***************************************************************************/

bool    alignment_behind_call(vcf_call_t *vcf_call, sam_alignment_t *sam_alignment)

{
    if ( (strcmp(VCF_CHROMOSOME(vcf_call), SAM_RNAME(sam_alignment)) > 0) ||
	 (VCF_POS(vcf_call) >=
	  SAM_POS(sam_alignment) + SAM_SEQ_LEN(sam_alignment)) )
	return true;
    else
	return false;
}


/***************************************************************************
 *  Description:
 *      Update allele counts
 *
 *  History: 
 *  Date        Name        Modification
 *  2020-05-26  Jason Bacon Begin
 ***************************************************************************/

void    update_allele_count(vcf_call_t *vcf_call, sam_alignment_t *sam_alignment,
			   FILE *vcf_out_stream)

{
    unsigned char   allele;
    
    allele = SAM_SEQ(sam_alignment)[VCF_POS(vcf_call) - SAM_POS(sam_alignment)];
#ifdef DEBUG
    /*
    fprintf(stderr, "===\n%s pos=%zu len=%zu %s\n",
	    SAM_RNAME(sam_alignment), SAM_POS(sam_alignment),
	    SAM_SEQ_LEN(sam_alignment), SAM_SEQ(sam_alignment));
    */
    fprintf(stderr, "Found allele %c for call %s,%zu in seq %s,%zu.\n",
	    allele, VCF_CHROMOSOME(vcf_call), VCF_POS(vcf_call),
	    SAM_RNAME(sam_alignment), SAM_POS(sam_alignment));
    fputs("===\n", stderr);
    putc(allele, vcf_out_stream);
#endif
    if ( allele == *VCF_REF(vcf_call) )
	++VCF_REF_COUNT(vcf_call);
    else if ( allele == *VCF_ALT(vcf_call) )
	++VCF_ALT_COUNT(vcf_call);
    else
	++VCF_OTHER_COUNT(vcf_call);
}

