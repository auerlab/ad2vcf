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
    static sam_buff_t      sam_buff;    // Watch out for small stack sizes
    FILE            *vcf_in_stream,
		    *vcf_out_stream;
    vcf_call_t      vcf_call;
    bool            xz = false,
		    more_alignments;
    size_t          previous_vcf_pos,
		    total_alleles,
		    depth,
		    depth_sum = 0;
    ad2vcf_stats_t  stats = AD2VCF_STATS_INIT;
    char            cmd[CMD_MAX + 1],
		    vcf_out_filename[PATH_MAX + 1],
		    previous_vcf_chromosome[VCF_CHROMOSOME_MAX_CHARS + 1] = "",
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

    printf("\nProcessing \"%s\":\n\n", vcf_filename);
    
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

    vcf_call_init(&vcf_call);
    sam_buff_init(&sam_buff);
    
    while ( vcf_read_ss_call(vcf_in_stream, &vcf_call, VCF_SAMPLE_MAX_CHARS) == VCF_OK )
    {
	++stats.total_vcf_calls;
	
#ifdef DEBUG
	fprintf(stderr, "\n=========================\n");
	fprintf(stderr, "New VCF call: %s, %zu\n",
		VCF_CHROMOSOME(&vcf_call), VCF_POS(&vcf_call));
	fprintf(stderr, "=========================\n");
#endif

	/* Make sure VCF calls are sorted */
	if ( strcmp(VCF_CHROMOSOME(&vcf_call), previous_vcf_chromosome) == 0 )
	{
	    if ( VCF_POS(&vcf_call) < previous_vcf_pos )
		 vcf_out_of_order(&vcf_call, previous_vcf_chromosome,
				  previous_vcf_pos);
	    else
		previous_vcf_pos = VCF_POS(&vcf_call);
	}
	else if ( chromosome_name_cmp(VCF_CHROMOSOME(&vcf_call),
				      previous_vcf_chromosome) < 0 )
	{
	    vcf_out_of_order(&vcf_call, previous_vcf_chromosome,
			     previous_vcf_pos);
	}
	else
	{
	    printf("Starting VCF chromosome %s.\n",
		    VCF_CHROMOSOME(&vcf_call));
	    fflush(stdout);
	    strlcpy(previous_vcf_chromosome, VCF_CHROMOSOME(&vcf_call),
		    VCF_CHROMOSOME_MAX_CHARS);
	    previous_vcf_pos = VCF_POS(&vcf_call);
	}
	
	/* Skip SAM alignments that don't include this position */
	more_alignments = skip_upstream_alignments(&vcf_call, sam_stream,
					       &sam_buff, vcf_out_stream,
					       &stats);
	
	/* Scan SAM alignments that include this position and count alleles */
	if ( more_alignments )
	    allelic_depth(&vcf_call, sam_stream, &sam_buff, vcf_out_stream, &stats);
	
	depth = VCF_REF_COUNT(&vcf_call) + VCF_ALT_COUNT(&vcf_call);
	depth_sum += depth;
	if ( depth < stats.min_depth )
	    stats.min_depth = depth;
	if ( depth > stats.max_depth )
	    stats.max_depth = depth;
	
	/* Output record with allelic depth */
#ifdef DEBUG
	fputc('\n', vcf_out_stream);
#endif
	/* Compute stats on phred scores */
	/*
	qsort(VCF_PHREDS(&vcf_call), VCF_PHRED_COUNT(&vcf_call), 1,
	      (int (*)(const void *, const void *))uchar_cmp);
	*/
	
	//fprintf(stderr, "%s\n", VCF_PHREDS(&vcf_call));
	fprintf(vcf_out_stream,
		"%s\t%zu\t.\t%s\t%s\t.\t.\t.\t%s:AD:DP\t%s:%u,%u,%u:%u\n",
		VCF_CHROMOSOME(&vcf_call), VCF_POS(&vcf_call),
		VCF_REF(&vcf_call),
		VCF_ALT(&vcf_call),
		VCF_FORMAT(&vcf_call),
		VCF_SINGLE_SAMPLE(&vcf_call),
		VCF_REF_COUNT(&vcf_call),
		VCF_ALT_COUNT(&vcf_call),
		VCF_OTHER_COUNT(&vcf_call),
		VCF_REF_COUNT(&vcf_call) + VCF_ALT_COUNT(&vcf_call));

	// vcf_phred_blank(&vcf_call);
    }

    printf("\nFinal statistics:\n\n");
    printf("%zu VCF calls processed\n", stats.total_vcf_calls);
    printf("%zu SAM alignments processed\n",
	    stats.total_sam_alignments);
    printf("Max buffered alignments: %zu\n", sam_buff.max_count);
    printf("%zu SAM alignments discarded (%zu%%)\n",
	    stats.discarded_sam_alignments, 
	    stats.discarded_sam_alignments * 100 / stats.total_sam_alignments);
    if ( stats.discarded_sam_alignments != 0 )
	printf("MAPQ min = %zu  max = %zu  mean = %f\n",
		stats.min_discarded_score, stats.max_discarded_score,
		(double)stats.discarded_score_sum / stats.discarded_sam_alignments);
    total_alleles = stats.total_ref_alleles + stats.total_alt_alleles +
		    stats.total_other_alleles;
    printf("%zu total REF alleles (%zu%%)\n",
	    stats.total_ref_alleles,
	    stats.total_ref_alleles * 100 / total_alleles);
    printf("%zu total ALT alleles (%zu%%)\n",
	    stats.total_alt_alleles,
	    stats.total_alt_alleles * 100 / total_alleles);
    printf("%zu total OTHER alleles (%zu%%)\n",
	    stats.total_other_alleles,
	    stats.total_other_alleles * 100 / total_alleles);
    printf("Min depth = %zu\n", stats.min_depth);
    printf("Max depth = %zu\n", stats.max_depth);
    printf("Mean depth = %f\n",
	    (double)depth_sum / stats.total_vcf_calls);

    if ( xz )
	pclose(vcf_in_stream);
    else
	fclose(vcf_in_stream);
    
    return EX_OK;
}


/***************************************************************************
 *  Description:
 *      Skip alignments upstream of the given VCF call chromosome and pos
 *
 *  History: 
 *  Date        Name        Modification
 *  2020-05-26  Jason Bacon Begin
 ***************************************************************************/

bool    skip_upstream_alignments(vcf_call_t *vcf_call, FILE *sam_stream,
			     sam_buff_t *sam_buff, FILE *vcf_out_stream,
			     ad2vcf_stats_t *stats)

{
    size_t          c;
    bool            ma = true;
    // static so sam_alignment_read() won't keep reallocating seq
    static sam_alignment_t sam_alignment = SAM_ALIGNMENT_INIT;

    if ( sam_alignment.seq == NULL )
	sam_alignment_init(&sam_alignment, SAM_SEQ_MAX_CHARS);
    
    /*
     *  Check and discard already buffered alignments upstream of the given
     *  VCF call.  They will be useless to subsequent calls as well
     *  since the calls must be sorted in ascending order.
     */
    for (c = 0; (c < sam_buff->count) &&
		alignment_upstream_of_call(vcf_call, sam_buff->alignments[c]); ++c)
    {
#ifdef DEBUG
	fprintf(stderr, "skip(): Unbuffering alignment #%zu %s,%zu upstream of variant %s,%zu\n",
		c, SAM_RNAME(sam_buff->alignments[c]),
		SAM_POS(sam_buff->alignments[c]),
		VCF_CHROMOSOME(vcf_call), VCF_POS(vcf_call));
#endif
    }
    
    /* If anything to unbuffer, shift */
    if ( c > 0 )
	sam_buff_shift(sam_buff, c);
    
    /*
     *  Read alignments from the stream until we find one that's not upstream of
     *  this VCF call, i.e. not overlapping and at a lower position or
     *  chromosome.
     */
    if ( sam_buff->count == 0 )
    {
	while ( (ma = sam_alignment_read(sam_stream, &sam_alignment)) )
	{
	    ++stats->total_sam_alignments;
	    /*
	    fprintf(stderr, "sam_alignment_read(): %s,%zu,%zu,%zu,%u\n",
		    SAM_RNAME(&sam_alignment), SAM_POS(&sam_alignment),
		    SAM_SEQ_LEN(&sam_alignment), SAM_QUAL_LEN(&sam_alignment),
		    SAM_MAPQ(&sam_alignment));
	    */
	    if ( SAM_MAPQ(&sam_alignment) < MAPQ_MIN )
		stats_update_discarded(stats, &sam_alignment);
	    else
	    {
		/*
		 *  We're done when we find an alignment overlapping or after
		 *  the VCF call
		 */
		if ( ! alignment_upstream_of_call(vcf_call, &sam_alignment) )
		    break;
#ifdef DEBUG
		else
		    fprintf(stderr, "skip(): Skipping new alignment %s,%zu upstream of variant %s,%zu\n",
			    SAM_RNAME(&sam_alignment), SAM_POS(&sam_alignment),
			    VCF_CHROMOSOME(vcf_call), VCF_POS(vcf_call));
#endif
	    }
	}
#ifdef DEBUG
	fprintf(stderr, "skip(): Buffering alignment #%zu %s,%zu,%zu\n",
		    sam_buff->count,
		    SAM_RNAME(&sam_alignment), SAM_POS(&sam_alignment),
		    SAM_SEQ_LEN(&sam_alignment));
#endif
	sam_buff_add_alignment(sam_buff, &sam_alignment);
    }

    return ma;
}


/***************************************************************************
 *  Description:
 *      Scan alignments in the SAM stream that encompass the given variant
 *      chromosome and position and update ref_count and alt_count for the
 *      VCF call.
 *
 *  History: 
 *  Date        Name        Modification
 *  2020-05-26  Jason Bacon Begin
 ***************************************************************************/

bool    allelic_depth(vcf_call_t *vcf_call, FILE *sam_stream,
		      sam_buff_t *sam_buff, FILE *vcf_out_stream,
		      ad2vcf_stats_t *stats)

{
    size_t          c;
    bool            ma = true, overlapping = true;
    // static so sam_alignment_read() won't keep reallocating seq
    static sam_alignment_t sam_alignment = SAM_ALIGNMENT_INIT;

    if ( sam_alignment.seq == NULL )
	sam_alignment_init(&sam_alignment, SAM_SEQ_MAX_CHARS);

    /* Check and discard already buffered alignments */
    for (c = 0; (c < sam_buff->count) &&
		(overlapping = call_in_alignment(vcf_call, sam_buff->alignments[c]));
		++c)
    {
#ifdef DEBUG
	fprintf(stderr, "depth(): Counting buffered alignment #%zu %s,%zu containing call %s,%zu\n",
		c, SAM_RNAME(sam_buff->alignments[c]),
		SAM_POS(sam_buff->alignments[c]),
		VCF_CHROMOSOME(vcf_call), VCF_POS(vcf_call));
#endif
	update_allele_count(vcf_call, sam_buff->alignments[c], vcf_out_stream, stats);
    }
    
    if ( (c == 0) || overlapping )
    {
	/* Read and buffer more alignments from the stream */
	while ( (ma = sam_alignment_read(sam_stream, &sam_alignment)) )
	{
	    ++stats->total_sam_alignments;
	    /*
	    fprintf(stderr, "sam_alignment_read(): Read %s,%zu,%zu,%u\n",
		    SAM_RNAME(&sam_alignment), SAM_POS(&sam_alignment),
		    SAM_SEQ_LEN(&sam_alignment), SAM_MAPQ(&sam_alignment));
	    */
	    if ( SAM_MAPQ(&sam_alignment) < MAPQ_MIN )
		stats_update_discarded(stats, &sam_alignment);
	    else
	    {
#ifdef DEBUG
		fprintf(stderr, "depth(): Buffering new alignment #%zu %s,%zu,%zu\n",
			sam_buff->count, SAM_RNAME(&sam_alignment),
			SAM_POS(&sam_alignment), SAM_SEQ_LEN(&sam_alignment));
#endif
		sam_buff_add_alignment(sam_buff, &sam_alignment);
		
		if ( call_in_alignment(vcf_call, &sam_alignment) )
		{
#ifdef DEBUG
		    fprintf(stderr, "depth(): Counting new alignment %s,%zu containing call %s,%zu\n",
			    SAM_RNAME(&sam_alignment), SAM_POS(&sam_alignment),
			    VCF_CHROMOSOME(vcf_call), VCF_POS(vcf_call));
#endif
		    update_allele_count(vcf_call, &sam_alignment, vcf_out_stream, stats);
		}
		else
		{
#ifdef DEBUG
		    fprintf(stderr, "depth(): Does not contain call %s,%zu\n",
			    VCF_CHROMOSOME(vcf_call), VCF_POS(vcf_call));
#endif
		    break;
		}
	    }
	}
    }

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
 *      Determine SAM alignment is completely upstream of a VCF call position,
 *      i.e. not overlapping and at a lower position or chromosome.
 *
 *  History: 
 *  Date        Name        Modification
 *  2020-05-26  Jason Bacon Begin
 ***************************************************************************/

bool    alignment_upstream_of_call(vcf_call_t *vcf_call, sam_alignment_t *alignment)

{
    /*fprintf(stderr, "alignment_upstream_of_call(): %s,%zu,%zu %s,%zu\n",
	    SAM_RNAME(sam_alignment),SAM_POS(sam_alignment),
	    SAM_SEQ_LEN(sam_alignment),
	    VCF_CHROMOSOME(vcf_call),VCF_POS(vcf_call));*/
    if ( (SAM_POS(alignment) + SAM_SEQ_LEN(alignment) <= VCF_POS(vcf_call)) &&
	  (strcmp(SAM_RNAME(alignment), VCF_CHROMOSOME(vcf_call)) == 0) )
	return true;
    else if ( chromosome_name_cmp(SAM_RNAME(alignment), VCF_CHROMOSOME(vcf_call)) < 0 )
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
			   FILE *vcf_out_stream, ad2vcf_stats_t *stats)

{
    unsigned char   allele;
    size_t          position_in_sequence;
    
    position_in_sequence = VCF_POS(vcf_call) - SAM_POS(sam_alignment);
    allele = SAM_SEQ(sam_alignment)[position_in_sequence];
    
    /*fprintf(stderr, "%zu %zu %zu\n", position_in_sequence,
	    SAM_QUAL_LEN(sam_alignment), SAM_SEQ_LEN(sam_alignment));*/
    
#if 0
    if ( SAM_QUAL_LEN(sam_alignment) == SAM_SEQ_LEN(sam_alignment) )
    {
	phred = SAM_QUAL(sam_alignment)[position_in_sequence];
	if ( phred < PHRED_BASE + PHRED_MIN )
	{
	    ++stats->discarded_bases;
#ifdef DEBUG
	    fprintf(stderr,
		    "Discarding low-quality base: %s,%zu,%zu = %u ('%c')\n",
		    SAM_RNAME(sam_alignment), SAM_POS(sam_alignment),
		    position_in_sequence, phred - PHRED_BASE, phred);
#endif
	    return;
	}
	vcf_phred_add(vcf_call, phred);
    }
#endif
    
#ifdef DEBUG
    char            *atype;

    atype = allele == *VCF_REF(vcf_call) ? "ref" :
	allele == *VCF_ALT(vcf_call) ? "alt" : "other";
    fprintf(stderr, "Found \"%s\" allele %c at pos %zu in seq %s,%zu for call %s,%zu\n",
	    atype, allele,
	    VCF_POS(vcf_call) - SAM_POS(sam_alignment) + 1,
	    SAM_RNAME(sam_alignment), SAM_POS(sam_alignment),
	    VCF_CHROMOSOME(vcf_call), VCF_POS(vcf_call));
    fputs("===\n", stderr);
    putc(allele, vcf_out_stream);
#endif
    if ( allele == *VCF_REF(vcf_call) )
    {
	++VCF_REF_COUNT(vcf_call);
	++stats->total_ref_alleles;
    }
    else if ( allele == *VCF_ALT(vcf_call) )
    {
	++VCF_ALT_COUNT(vcf_call);
	++stats->total_alt_alleles;
    }
    else
    {
	++VCF_OTHER_COUNT(vcf_call);
	++stats->total_other_alleles;
    }
}


/***************************************************************************
 *  Description:
 *      Check order of sam alignments being read into buffer.
 *
 *  History: 
 *  Date        Name        Modification
 *  2020-05-27  Jason Bacon Begin
 ***************************************************************************/

void    sam_buff_check_order(sam_buff_t *sam_buff, sam_alignment_t *sam_alignment)

{
    /*fprintf(stderr, "Previous SAM: %s %zu, Current SAM: %s %zu\n",
	    sam_buff->previous_rname, sam_buff->previous_pos,
	    sam_alignment->rname, sam_alignment->pos);*/
    if ( strcmp(sam_alignment->rname, sam_buff->previous_rname) == 0 )
    {
	// Silly to assign when ==, but sillier to add another check
	if (sam_alignment->pos < sam_buff->previous_pos )
	    sam_buff_out_of_order(sam_buff, sam_alignment);
	else
	    sam_buff->previous_pos = sam_alignment->pos;
    }
    else if ( chromosome_name_cmp(sam_alignment->rname, sam_buff->previous_rname) < 0 )
	sam_buff_out_of_order(sam_buff, sam_alignment);
    else
    {
	strlcpy(sam_buff->previous_rname, sam_alignment->rname, SAM_RNAME_MAX_CHARS);
	sam_buff->previous_pos = sam_alignment->pos;
    }
}


/***************************************************************************
 *  Description:
 *      Initialize a SAM alignment buffer
 *
 *  History: 
 *  Date        Name        Modification
 *  2020-05-27  Jason Bacon Begin
 ***************************************************************************/

void    sam_buff_init(sam_buff_t *sam_buff)

{
    size_t  c;
    
    sam_buff->count = 0;
    sam_buff->max_count = 0;
    sam_buff->previous_pos = 0;
    *sam_buff->previous_rname = '\0';
    sam_buff->buff_size = SAM_BUFF_START_SIZE;
    /*
     *  Dynamically allocating the pointers is probably senseless since they
     *  take very little space compared to the alignment data.  By the time
     *  the pointer array takes a significant amount of RAM, you're probably
     *  already thrashing to accomodate the sequence data.  The pointer array
     *  size is capped by SAM_BUFF_MAX_SIZE to prevent memory exhaustion.
     *  We may save a few megabytes with this, though.
     */
    sam_buff->alignments =
	(sam_alignment_t **)malloc(sam_buff->buff_size *
				   sizeof(sam_alignment_t **));
    for (c = 0; c < sam_buff->buff_size; ++c)
	sam_buff->alignments[c] = NULL;
}


/***************************************************************************
 *  Description:
 *      Add a new alignment to the buffer
 *
 *  History: 
 *  Date        Name        Modification
 *  2020-05-27  Jason Bacon Begin
 ***************************************************************************/

void    sam_buff_add_alignment(sam_buff_t *sam_buff,
			       sam_alignment_t *sam_alignment)

{
    size_t  old_buff_size,
	    c;

    sam_buff_check_order(sam_buff, sam_alignment);
    
    // Just allocate the static fields, sam_alignment_copy() does the rest
    if ( sam_buff->alignments[sam_buff->count] == NULL )
    {
	//fprintf(stderr, "Allocating alignment #%zu\n", sam_buff->count);
	sam_buff->alignments[sam_buff->count] = malloc(sizeof(sam_alignment_t));
	if ( sam_buff->alignments[sam_buff->count] == NULL )
	{
	    fprintf(stderr, "sam_buff_add_alignment(): malloc() failed.\n");
	    exit(EX_UNAVAILABLE);
	}
	// Redundant to sam_alignment_copy()
	// sam_alignment_init(sam_buff->alignments[sam_buff->count], 0);
    }
    else
	sam_alignment_free(sam_buff->alignments[sam_buff->count]);
    
    sam_alignment_copy(sam_buff->alignments[sam_buff->count], sam_alignment);
    
    ++sam_buff->count;

    if ( sam_buff->count > sam_buff->max_count )
    {
	sam_buff->max_count = sam_buff->count;
	// fprintf(stderr, "sam_buff->max_count = %zu\n", sam_buff->max_count);
    }
    
    if ( sam_buff->count == SAM_BUFF_MAX_SIZE )
    {
	fprintf(stderr,
		"sam_buff_add_alignment(): Hit SAM_BUFF_MAX_SIZE=%u.\n",
		SAM_BUFF_MAX_SIZE);
	fprintf(stderr, "Terminating to prevent runaway memory use.\n");
	fprintf(stderr, "Check your SAM input.\n");
	exit(EX_DATAERR);
    }
    
    if ( sam_buff->count == sam_buff->buff_size )
    {
	fprintf(stderr,
		"sam_buff_add_alignment(): Hit buff_size=%zu, doubling buffer size.\n",
		sam_buff->buff_size);
	fprintf(stderr, "RNAME: %s  POS: %zu  LEN: %zu\n",
		SAM_RNAME(sam_alignment), SAM_POS(sam_alignment),
		SAM_SEQ_LEN(sam_alignment));
	old_buff_size = sam_buff->buff_size;
	sam_buff->buff_size *= 2;
	sam_buff->alignments =
	    (sam_alignment_t **)realloc(sam_buff->alignments,
					sam_buff->buff_size *
					sizeof(sam_alignment_t **));
	for (c = old_buff_size; c < sam_buff->buff_size; ++c)
	    sam_buff->alignments[c] = NULL;
    }
}


/***************************************************************************
 *  Description:
 *      Explain SAM input sort error and exit.
 *
 *  History: 
 *  Date        Name        Modification
 *  2020-05-27  Jason Bacon Begin
 ***************************************************************************/

void    sam_buff_out_of_order(sam_buff_t *sam_buff, sam_alignment_t *sam_alignment)

{
    fprintf(stderr, "ad2vcf: Error: SAM input must be sorted by chromosome and then position.\n");
    fprintf(stderr, "Found %s,%zu after %s,%zu.\n",
	    SAM_RNAME(sam_alignment), SAM_POS(sam_alignment),
	    sam_buff->previous_rname, sam_buff->previous_pos);
    exit(EX_DATAERR);
}


/***************************************************************************
 *  Description:
 *      Explain VCF input sort error and exit.
 *
 *  History: 
 *  Date        Name        Modification
 *  2020-05-27  Jason Bacon Begin
 ***************************************************************************/

void    vcf_out_of_order(vcf_call_t *vcf_call,
			 char *previous_chromosome, size_t previous_pos)

{
    fprintf(stderr, "ad2vcf: Error: VCF input must be sorted by chromosome and then position.\n");
    fprintf(stderr, "Found %s,%zu after %s,%zu.\n",
	    VCF_CHROMOSOME(vcf_call), VCF_POS(vcf_call),
	    previous_chromosome, previous_pos);
    exit(EX_DATAERR);
}


/***************************************************************************
 *  Description:
 *      Free an element of the SAM alignment array
 *
 *  History: 
 *  Date        Name        Modification
 *  2020-05-29  Jason Bacon Begin
 ***************************************************************************/

void    sam_buff_free_alignment(sam_buff_t *sam_buff, size_t c)

{
    sam_alignment_free(sam_buff->alignments[c]);
    sam_alignment_init(sam_buff->alignments[c], 0);
    if ( sam_buff->alignments[c] != NULL )
    {
	free(sam_buff->alignments[c]);
	sam_buff->alignments[c] = NULL;
    }
}


/***************************************************************************
 *  Description:
 *      Shift SAM array elements up c positions.
 *
 *  History: 
 *  Date        Name        Modification
 *  2020-05-29  Jason Bacon Begin
 ***************************************************************************/

void    sam_buff_shift(sam_buff_t *sam_buff, size_t c)

{
    size_t  c2;

    /*
     *  FIXME: A circular queue would be more efficient, but won't matter
     *  much to the bottom line to avoid shifting a few pointers on top
     *  of everything else going on here
     */
    
    /* Make sure elements to be removed are freed */
    for (c2 = 0; c2 < c; ++c2)
	sam_buff_free_alignment(sam_buff, c2);

    /* Shift elements */
    for (c2 = 0; c2 < sam_buff->count - c; ++c2)
	sam_buff->alignments[c2] = sam_buff->alignments[c2 + c];
    
    /* Clear vacated elements */
    while ( c2 < sam_buff->count )
	sam_buff->alignments[c2++] = NULL;
    
    sam_buff->count -= c;
}


int     uchar_cmp(unsigned char *c1, unsigned char *c2)

{
    return *c1 - *c2;
}


void    stats_update_discarded(ad2vcf_stats_t *stats,
			       sam_alignment_t *sam_alignment)

{
    ++stats->discarded_sam_alignments;
    stats->discarded_score_sum += SAM_MAPQ(sam_alignment);
    if ( SAM_MAPQ(sam_alignment) < stats->min_discarded_score )
	stats->min_discarded_score = SAM_MAPQ(sam_alignment);
    if ( SAM_MAPQ(sam_alignment) > stats->max_discarded_score )
	stats->max_discarded_score = SAM_MAPQ(sam_alignment);

    fprintf(stderr, "Discarding low quality read: %s,%zu MAPQ=%u\n",
	    SAM_RNAME(sam_alignment), SAM_POS(sam_alignment),
	    SAM_MAPQ(sam_alignment));
}
