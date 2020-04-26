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
    FILE            *vcf_stream,
		    *allele_stream;
    sam_alignment_t sam_alignment;
    int             more_alignments,
		    allele,
		    c;
    bool            xz = false;
    size_t          vcf_pos = 0,
		    vcf_calls_read = 0,
		    alignments_read = 1,
		    previous_vcf_pos = 0,
		    previous_alignment_pos = 0;
    char            cmd[CMD_MAX + 1],
		    allele_filename[PATH_MAX + 1],
		    *vcf_chromosome,
		    previous_vcf_chromosome[VCF_CHROMOSOME_MAX_CHARS + 1] = "",
		    previous_sam_rname[SAM_RNAME_MAX + 1] = "",
		    *ext;
    const char      *vcf_filename = argv[1];
    /*
     *  Linux default stack size is 8 MiB and VCF_INFO_MAX_CHARS has to be
     *  large to accommodate dpGap BCFs.  Make this static to avoid stack
     *  overflow (manifests as a seg fault on Linux).
     */
    static vcf_calls_for_position_t    vcf_calls_for_position;
    
    xz = ((ext = strstr(vcf_filename,".xz")) != NULL) && (ext[3] == '\0');
    if ( xz )
    {
	snprintf(cmd, CMD_MAX, "unxz -c %s", vcf_filename);
	vcf_stream = popen(cmd, "r");
    }
    else
	vcf_stream = fopen(vcf_filename, "r");
    
    if ( vcf_stream == NULL )
    {
	fprintf(stderr, "%s: Cannot open %s: %s\n", argv[0], argv[1],
	    strerror(errno));
	exit(EX_NOINPUT);
    }

    if ( (ext = strstr(vcf_filename, ".vcf")) == NULL )
    {
	fprintf(stderr, "%s: Input filename must contain \".vcf\".\n", argv[0]);
	exit(EX_DATAERR);
    }
    
    // Insert "-ad" before ".vcf"
    *ext = '\0';
    snprintf(allele_filename, PATH_MAX, "%s-ad.%s", vcf_filename, ext+1);

    if ( xz )
    {
	snprintf(cmd, CMD_MAX, "xz -c > %s", allele_filename);
	allele_stream = popen(cmd, "w");
    }
    else
	allele_stream = fopen(allele_filename, "w");
    
    if ( allele_stream == NULL )
    {
	fprintf(stderr, "%s: Cannot open %s: %s\n", argv[0], allele_filename,
	    strerror(errno));
	exit(EX_CANTCREAT);
    }
    
    /*
     *  Read in VCF fields
     */
    
    // Prime loop by reading first SAM alignment
    more_alignments = sam_read_alignment(argv, sam_stream, &sam_alignment);
    
    /*
     *  Check each SAM alignment against all ALT alleles.  There may be
     *  consecutive VCF calls with the same position, so we read in all of
     *  them at once and match each SAM read to every one.
     */
    while ( more_alignments &&
	    ((vcf_read_calls_for_position(vcf_stream,
				      &vcf_calls_for_position)) > 0) )
    {
	// chromosome and position are the same for all calls
	vcf_pos = vcf_calls_for_position.call[0].pos;
	vcf_chromosome = vcf_calls_for_position.call[0].chromosome;
	// FIXME: Can this be done before the loop?
	// This if check is a tiny performance hit
	if ( *previous_vcf_chromosome == '\0' )
	    strlcpy(previous_vcf_chromosome, vcf_chromosome,
		    VCF_CHROMOSOME_MAX_CHARS);
	
	/*
	 *  VCF input must be sorted by chromosome, then position.
	 *  If current position < previous and chromosome is the same,
	 *  then the input is not sorted.
	 */
	if ( vcf_pos < previous_vcf_pos )
	{
	    if ( strcmp(vcf_chromosome, previous_vcf_chromosome) == 0 )
	    {
		fprintf(stderr, "%s: VCF input must be sorted first by chromosome then position.\n", argv[0]);
		exit(EX_DATAERR);
	    }
	    else
	    {
		fprintf(stderr, "INFO: Finished chromosome %s: %zu calls processed.\n",
			previous_vcf_chromosome, vcf_calls_read);
		// Begin next chromosome, reset pos
		strlcpy(previous_vcf_chromosome, vcf_chromosome,
			VCF_CHROMOSOME_MAX_CHARS);
		previous_vcf_pos = vcf_pos;
		vcf_calls_read = 0;
	    }
	}
	else
	    previous_vcf_pos = vcf_pos;
	
	// Debug
	// fprintf(stderr, "VCF call %s %zu\n", vcf_chromosome, vcf_pos);
	vcf_calls_read += vcf_calls_for_position.count;
	if ( vcf_calls_for_position.count > 1 )
	    fprintf(stderr, "INFO: %zu calls at chromosome %s pos %zu.\n",
		    vcf_calls_for_position.count, vcf_chromosome, vcf_pos);

	// Skip remaining alignments for previous chromosome after VCF
	// chromosome changes
	while ( more_alignments &&
		(chromosome_name_cmp(SAM_RNAME(&sam_alignment), vcf_chromosome) < 0) )
	{
	    // Debug
	    /*
	    fprintf(stderr, "Skipping %s %zu on the way to VCF %s %zu\n",
		    SAM_RNAME(&sam_alignment), SAM_POS(&sam_alignment),
		    vcf_chromosome, vcf_pos);
	    */
	    more_alignments =
		sam_read_alignment(argv, sam_stream, &sam_alignment);
	    ++alignments_read;
	}
	
#ifdef DEBUG
	fprintf(allele_stream, "# %s %zu ", vcf_chromosome, vcf_pos);
#endif

	/*
	 *  Now check all SAM alignments for the same chromosome with sequence
	 *  starting positions <= the VCF call position. Both SAM and VCF
	 *  should be sorted by chromosome first and then position, so when we
	 *  encounter a SAM alignment with a different chromosome or a position
	 *  beyond the VCF position, we're done with this VCF call.
	 *
	 *  Do integer position compare before strcmp().  It's
	 *  less intuitive to check the second sort key first, but faster
	 *  since the strcmp() is unnecessary when the integer
	 *  comparison is false.
	 */
	while ( more_alignments &&
		(SAM_POS(&sam_alignment) <= vcf_pos) &&
		(strcmp(SAM_RNAME(&sam_alignment), vcf_chromosome) == 0) )
	{
	    /*
	     *  We know at this point that the VCF call position is downstream
	     *  from the start of the SAM alignment sequence.  Is it also
	     *  upstream of the end?  If so, record the exact position and
	     *  allele.
	     */
	    if ( vcf_pos < SAM_POS(&sam_alignment) + SAM_SEQ_LEN(&sam_alignment) )
	    {
		allele = SAM_SEQ(&sam_alignment)[vcf_pos - SAM_POS(&sam_alignment)];
#ifdef DEBUG
		fprintf(stderr, "===\n%s pos=%zu len=%zu %s\n",
			SAM_RNAME(&sam_alignment), SAM_POS(&sam_alignment),
			SAM_SEQ(&sam_alignment)_len, SAM_SEQ(&sam_alignment));
		fprintf(stderr, "Found allele %c (%d) for call pos %zu on %s aligned seq starting at %zu.\n",
			allele, allele, vcf_pos, vcf_chromosome,
			SAM_POS(&sam_alignment));
		fprintf(stderr, "Calls: %zu  Alignments: %zu\n",
			vcf_calls_read, alignments_read);
		putc(allele, allele_stream);
#endif
		
		for (c = 0; c < vcf_calls_for_position.count; ++c)
		{
		    if ( allele == *VCF_REF(&vcf_calls_for_position.call[c]) )
			++VCF_REF_COUNT(&vcf_calls_for_position.call[c]);
		    else if ( allele == *VCF_ALT(&vcf_calls_for_position.call[c]) )
			++VCF_ALT_COUNT(&vcf_calls_for_position.call[c]);
		    else
			++VCF_OTHER_COUNT(&vcf_calls_for_position.call[c]);
		}
	    }
	    
	    more_alignments = sam_read_alignment(argv, sam_stream, &sam_alignment);
	    ++alignments_read;
	    
	    /*
	     *  SAM input must be sorted by rname (chromosome), then position.
	     *  If current position < previous and chromosome is the same,
	     *  then the SAM input is not sorted.
	     */
	    if ( SAM_POS(&sam_alignment) < previous_alignment_pos )
	    {
		if ( strcmp(SAM_RNAME(&sam_alignment), previous_sam_rname) == 0 )
		{
		    fprintf(stderr, "%s: SAM input is not sorted.\n", argv[0]);
		    exit(EX_DATAERR);
		}
		else
		{
		    // Begin next chromosome, reset pos
		    strlcpy(previous_sam_rname, SAM_RNAME(&sam_alignment), SAM_RNAME_MAX);
		    previous_alignment_pos = SAM_POS(&sam_alignment);
		}
	    }
	    else
		previous_alignment_pos = SAM_POS(&sam_alignment);
	}

#ifdef DEBUG
	putc('\n', allele_stream);
	/*
	fprintf(stderr, "===\nOut of alignments for VCF %s %zu\n",
		vcf_chromosome, vcf_pos);
	fprintf(stderr, "rname = %s pos=%zu\n",
		SAM_RNAME(&sam_alignment), SAM_POS(&sam_alignment));
	*/
#endif
	
	for (c = 0; c < vcf_calls_for_position.count; ++c)
	{
	    // FIXME: Use vcf_write_call()
	    // Haplohseq expects DP to be sum of AD values (P. Auer)
	    fprintf(allele_stream,
		    "%s\t%zu\t.\t%s\t%s\t.\t.\t.\t%s:AD:DP\t%s:%u,%u:%u\n",
		    vcf_chromosome, vcf_pos,
		    VCF_REF(&vcf_calls_for_position.call[c]),
		    VCF_ALT(&vcf_calls_for_position.call[c]),
		    VCF_FORMAT(&vcf_calls_for_position.call[c]),
		    VCF_SINGLE_SAMPLE(&vcf_calls_for_position.call[c]),
		    VCF_REF_COUNT(&vcf_calls_for_position.call[c]),
		    VCF_ALT_COUNT(&vcf_calls_for_position.call[c]),
		    VCF_REF_COUNT(&vcf_calls_for_position.call[c]) +
			VCF_ALT_COUNT(&vcf_calls_for_position.call[c]));
		    // vcf_calls_for_position.call[c].other_count);
	}
    }
    
    fprintf(stderr, "INFO: Finished chromosome %s: %zu calls processed.\n",
	    vcf_chromosome, vcf_calls_read);
    
    // Debug
    fprintf(stderr, "Loop terminated with more_alignments = %d, new calls = %zu\n",
	    more_alignments, vcf_calls_for_position.count);
    fprintf(stderr, "===\nrname = %s pos=%zu len=%zu %s\n",
	    SAM_RNAME(&sam_alignment), SAM_POS(&sam_alignment),
	    SAM_SEQ_LEN(&sam_alignment), SAM_SEQ(&sam_alignment));
    fprintf(stderr, "vcf_pos = %zu, vcf_chromosome = %s\n",
	    vcf_pos, vcf_chromosome);
    
    if ( xz )
	pclose(vcf_stream);
    else
	fclose(vcf_stream);
    return EX_OK;
}
