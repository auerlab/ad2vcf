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
#include <sys/param.h>  // MIN()
#include <xtend.h>
#include <vcfio.h>
#include <sam-buff.h>
#include <biostring.h>
#include "ad2vcf.h"

int     main(int argc, const char *argv[])

{
    if ( argc != 3 )
	usage(argv);
    
    return ad2vcf(argv, stdin);
}


void    usage(const char *argv[])

{
    fprintf(stderr, "Usage: %s single-sample.vcf[.xz] minimum-MAPQ < file.sam\n", argv[0]);
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
    sam_buff_t      sam_buff;
    FILE            *vcf_in_stream,
		    *vcf_out_stream;
    vcf_call_t      vcf_call = VCF_CALL_INIT;
    bool            xz = false,
		    more_alignments;
    size_t          previous_vcf_pos,
		    total_alleles,
		    depth,
		    depth_sum = 0;
    ad2vcf_stats_t  stats;
    char            cmd[CMD_MAX + 1],
		    vcf_out_filename[PATH_MAX + 1],
		    previous_vcf_chromosome[BIO_CHROMOSOME_MAX_CHARS + 1] = "",
		    *ext,
		    *end;
    const char      *vcf_filename = argv[1];
    unsigned int    mapq_min;

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
    
    mapq_min = strtol(argv[2], &end, 10);
    if ( *end != '\0' )
    {
	fprintf(stderr, "%s: Invalid MAPQ minimum: %s\n", argv[0], argv[2]);
	exit(EX_USAGE);
    }
    
    sam_buff_init(&sam_buff, mapq_min);
    ad2vcf_stats_init(&stats);
    
    printf("\nProcessing \"%s\", MAPQ min = %u:\n\n", vcf_filename, mapq_min);
    
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

    vcf_call_init(&vcf_call, VCF_INFO_MAX_CHARS, VCF_FORMAT_MAX_CHARS,
		  VCF_SAMPLE_MAX_CHARS);
    
    while ( vcf_read_ss_call(vcf_in_stream, &vcf_call) == BIO_READ_OK )
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
		    BIO_CHROMOSOME_MAX_CHARS);
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

#ifdef DEBUG
    static sam_alignment_t sam_alignment = SAM_ALIGNMENT_INIT;
    if ( sam_alignment.seq == NULL )
	sam_alignment_init(&sam_alignment, SAM_SEQ_MAX_CHARS);

    // Debug discarded count
    puts("Gathering stats on trailing alignments...");
    while ( sam_alignment_read(sam_stream, &sam_alignment) )
    {
	SAM_BUF_INC_TOTAL_ALIGNMENTS(&sam_buff);
	SAM_BUF_INC_TRAILING_ALIGNMENTS(&sam_buff);
	if ( !sam_buff_alignment_ok(&sam_buff, &sam_alignment) )
	    SAM_BUF_INC_DISCARDED_TRAILING(&sam_buff);
    }
#endif

    printf("\nFinal statistics:\n\n");
    printf("%zu VCF calls processed\n", stats.total_vcf_calls);
    printf("%zu SAM alignments processed\n",
	    SAM_BUFF_TOTAL_ALIGNMENTS(&sam_buff));
    printf("Max buffered alignments: %zu\n",
	    SAM_BUFF_MAX_COUNT(&sam_buff));
    printf("%zu low MAPQ alignments discarded (%zu%%)\n",
	    SAM_BUFF_DISCARDED_ALIGNMENTS(&sam_buff),
	    SAM_BUFF_DISCARDED_ALIGNMENTS(&sam_buff) * 100 /
		SAM_BUFF_TOTAL_ALIGNMENTS(&sam_buff));
    printf("%zu unmapped alignments discarded (%zu%%)\n",
	    SAM_BUFF_UNMAPPED_ALIGNMENTS(&sam_buff),
	    SAM_BUFF_UNMAPPED_ALIGNMENTS(&sam_buff) * 100 /
		SAM_BUFF_TOTAL_ALIGNMENTS(&sam_buff));
#ifdef DEBUG
    printf("%zu SAM alignments beyond last call.\n",
	    SAM_BUFF_TRAILING_ALIGNMENTS(&sam_buff));
    printf("%zu trailing SAM alignments discarded (%zu%%)\n",
	    SAM_BUFF_DISCARDED_TRAILING(&sam_buff),
	    SAM_BUFF_DISCARDED_TRAILING(&sam_buff) * 100 /
		SAM_BUFF_TRAILING_ALIGNMENTS(&sam_buff));
#endif
    if ( SAM_BUFF_DISCARDED_ALIGNMENTS(&sam_buff) != 0 )
	printf("MAPQ min discarded = %zu  max discarded = %zu  mean = %f\n",
		SAM_BUFF_MIN_DISCARDED_SCORE(&sam_buff),
		SAM_BUFF_MAX_DISCARDED_SCORE(&sam_buff),
		(double)SAM_BUFF_DISCARDED_SCORE_SUM(&sam_buff) /
		    SAM_BUFF_DISCARDED_ALIGNMENTS(&sam_buff));
    printf("MAPQ min used = %zu  max used = %zu  mean = %f\n",
	    SAM_BUFF_MAPQ_LOW(&sam_buff), SAM_BUFF_MAPQ_HIGH(&sam_buff),
	    (double)SAM_BUFF_MAPQ_SUM(&sam_buff) / SAM_BUFF_READS_USED(&sam_buff));
    
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
    for (c = 0; (c < SAM_BUFF_BUFFERED_COUNT(sam_buff)) &&
	 vcf_call_downstream_of_alignment(vcf_call, SAM_BUFF_ALIGNMENTS(sam_buff, c));
	 ++c)
    {
#ifdef DEBUG
	fprintf(stderr, "skip(): Unbuffering alignment #%zu %s,%zu upstream of variant %s,%zu\n",
		c, SAM_RNAME(SAM_BUFF_ALIGNMENTS(sam_buff,c)),
		SAM_POS(SAM_BUFF_ALIGNMENTS(sam_buff,c)),
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
    if ( SAM_BUFF_BUFFERED_COUNT(sam_buff) == 0 )
    {
	while ( (ma = sam_alignment_read(sam_stream, &sam_alignment)) )
	{
	    SAM_BUFF_INC_TOTAL_ALIGNMENTS(sam_buff);
	    /*
	    fprintf(stderr, "sam_alignment_read(): %s,%zu,%zu,%zu,%u\n",
		    SAM_RNAME(&sam_alignment), SAM_POS(&sam_alignment),
		    SAM_SEQ_LEN(&sam_alignment), SAM_QUAL_LEN(&sam_alignment),
		    SAM_MAPQ(&sam_alignment));
	    */
	    if ( sam_buff_alignment_ok(sam_buff, &sam_alignment) )
	    {
		/*
		 *  We're done when we find an alignment overlapping or after
		 *  the VCF call
		 */
		if ( ! vcf_call_downstream_of_alignment(vcf_call, &sam_alignment) )
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
		    SAM_BUFF_BUFFERED_COUNT(sam_buff),
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
    for (c = 0; (c < SAM_BUFF_BUFFERED_COUNT(sam_buff)) &&
		(overlapping = vcf_call_in_alignment(vcf_call,
		    SAM_BUFF_ALIGNMENTS(sam_buff,c)));
		++c)
    {
#ifdef DEBUG
	fprintf(stderr, "depth(): Counting buffered alignment #%zu %s,%zu containing call %s,%zu\n",
		c, SAM_RNAME(SAM_BUFF_ALIGNMENTS(sam_buff,c)),
		SAM_POS(SAM_BUFF_ALIGNMENTS(sam_buff,c),
		VCF_CHROMOSOME(vcf_call), VCF_POS(vcf_call));
#endif
	update_allele_count(vcf_call, SAM_BUFF_ALIGNMENTS(sam_buff,c),
	    vcf_out_stream, stats);
    }
    
    if ( (c == 0) || overlapping )
    {
	/* Read and buffer more alignments from the stream */
	while ( (ma = sam_alignment_read(sam_stream, &sam_alignment)) )
	{
	    SAM_BUFF_INC_TOTAL_ALIGNMENTS(sam_buff);
	    /*
	    fprintf(stderr, "sam_alignment_read(): Read %s,%zu,%zu,%u\n",
		    SAM_RNAME(&sam_alignment), SAM_POS(&sam_alignment),
		    SAM_SEQ_LEN(&sam_alignment), SAM_MAPQ(&sam_alignment));
	    */
	    if ( sam_buff_alignment_ok(sam_buff, &sam_alignment) )
	    {
#ifdef DEBUG
		fprintf(stderr, "depth(): Buffering new alignment #%zu %s,%zu,%zu\n",
			SAM_BUFF_BUFFERED_COUNT(sam_buff),
			SAM_RNAME(&sam_alignment),
			SAM_POS(&sam_alignment), SAM_SEQ_LEN(&sam_alignment));
#endif
		sam_buff_add_alignment(sam_buff, &sam_alignment);
		
		if ( vcf_call_in_alignment(vcf_call, &sam_alignment) )
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


int     uchar_cmp(unsigned char *c1, unsigned char *c2)

{
    return *c1 - *c2;
}


void    ad2vcf_stats_init(ad2vcf_stats_t *stats)

{
    stats->total_vcf_calls = 0;
    stats->total_ref_alleles = 0;
    stats->total_alt_alleles = 0;
    stats->total_other_alleles = 0;
    stats->min_depth = SIZE_MAX;
    stats->max_depth = 0;
    stats->mean_depth = 0;
}
