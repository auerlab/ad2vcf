# ad2vcf

ad2vcf extracts allelic depth info from a SAM stream and adds it to a
corresponding single-sample VCF file.

SAM input is read via stdin and the VCF input file is taken as a command-line
argument.  This allows expensive BAM/CRAM decoding to occur in-parallel using
a pipe:

```sh
samtools view -@ 2 --input-fmt-option required_fields=0x218 \
	../SRR6990379/NWD102903.b38.irc.v1.cram \
	| ./ad2vcf file.vcf
```

ad2vcf is written entirely in C and attempts to optimize CPU, memory,
and disk access.  It does not inhale large amounts of data into RAM, so memory
use is trivial and it runs mostly from cache, making it very fast.

Processing a large CRAM file with human genome alignments against from the SRA
project against dbGaP VCF data takes about 20 minutes on a modern workstation
or server with ad2vcf averaging 4 to 5 MB (not GB) resident memory use.
Memory use will spike briefly due to alignment buffering when processing
regions where many alignments overlap multiple variant calls.

Both SAM and VCF inputs must be sorted first by chromosome and then by
read/call position.

ad2vcf is intended to build cleanly in any POSIX environment.  Please
don't hesitate to open an issue if you encounter problems on any
Unix-like system.

Primary development is done on FreeBSD with clang, but the code is frequently
tested on CentOS, Mac OS X, and NetBSD as well.  MS Windows is not supported,
unless using a POSIX environment such as Cygwin or Windows Subsystem for Linux.

Building and installing:

ad2vcf depends on [biolibc](https://github.com/auerlab/biolibc).

Set LOCALBASE to the prefix of lib/libbiolibc.a.  Default is ../local.
(See Makefile).

Set PREFIX to the prefix where you would like to install.  Default is
$LOCALBASE.

Then simply run

```sh
make install
```
