# ad2vcf

## Description

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

Both SAM and VCF inputs must be sorted first by chromosome and then by
read/call position.

## Design and Implementation

The code is organized following basic object-oriented design principals, but
implemented in C to minimize overhead and keep the source code accessible to
scientists who don't have time to master the complexities of C++.

Structures are treated as classes, with accessor macros and mutator functions
provided, so dependent applications and libraries need not access
structure members directly.  Since the C language cannot enforce this, it's
up to application programmers to exercise self-discipline.

ad2vcf is written entirely in C and attempts to optimize CPU, memory,
and disk access.  It does not inhale large amounts of data into RAM, so memory
use is trivial and it runs mostly from cache, making it very fast.

Processing a large CRAM file with human genome alignments against from the SRA
project against dbGaP VCF data takes about 20 minutes on a modern workstation
or server with ad2vcf averaging 4 to 5 MB (not GB) resident memory use.
Memory use will spike briefly due to alignment buffering when processing
regions where many alignments overlap multiple variant calls.

## Building and installing

ad2vcf is intended to build cleanly in any POSIX environment on
any CPU architecture.  Please
don't hesitate to open an issue if you encounter problems on any
Unix-like system.

Primary development is done on FreeBSD with clang, but the code is frequently
tested on Linux, MacOS, NetBSD, and OpenIndiana as well.  MS Windows is not supported,
unless using a POSIX environment such as Cygwin or Windows Subsystem for Linux.

The Makefile is designed to be friendly to package managers, such as
[Debian packages](https://www.debian.org/distrib/packages),
[FreeBSD ports](https://www.freebsd.org/ports/),
[MacPorts](https://www.macports.org/), [dreckly](https://github.com/drecklypkg/dreckly), etc.

End users should install using a package manager, to ensure that
dependencies are properly managed.

I maintain a FreeBSD port and a dreckly package, which is sufficient to install
cleanly on virtually any POSIX platform.  If you would like to see a
package in another package manager, please consider creating a package
yourself.  This will be one of the easiest packages in the collection and
hence a good vehicle to learn how to create packages.

Note that dreckly can be used by anyone, on virtually any POSIX operating
system, with or without administrator privileges.

For an overview of available package managers, see the
[Repology website](https://repology.org/).

### Installing ad2vcf on FreeBSD:

FreeBSD is a highly underrated platform for scientific computing, with over
2,000 scientific libraries and applications in the FreeBSD ports collection
(of more than 30,000 total), modern clang compiler, fully-integrated ZFS
filesystem, and renowned security, performance, and reliability.
FreeBSD has a somewhat well-earned reputation for being difficult to set up
and manage compared to user-friendly systems like [Ubuntu](https://ubuntu.com/).
However, if you're a little bit Unix-savvy, you can very quickly set up a
workstation, laptop, or VM using
[desktop-installer](http://www.acadix.biz/desktop-installer.php).
[GhostBSD](https://ghostbsd.org/) offers an experience very similar
to Ubuntu, but is built on FreeBSD rather than Debian Linux.  GhostBSD
packages lag behind FreeBSD ports slightly, but this is not generally
an issue and there are workarounds.

To install the binary package on FreeBSD:

```
pkg install ad2vcf
```

You can just as easily build and install from source.  This is useful for
FreeBSD ports with special build options, for building with non-portable
optimizations such as -march=native, and for 
[work-in-progress ports](https://github.com/outpaddling/freebsd-ports-wip),
for which binary packages are not yet maintained.

```
cd /usr/ports/biology/ad2vcf && env CFLAGS='-march=native -O2' make install
cd /usr/ports/wip/ad2vcf && make install
```

### Installing via dreckly

[Dreckly](https://github.com/drecklypkg/dreckly) is a cross-platform package manager that works on any Unix-like
platform. It is derived from pkgsrc, which is part of [NetBSD](https://www.netbsd.org/),, and well-supported
on [Illumos](https://illumos.org/), [MacOS](https://www.apple.com/macos/),
[RHEL](https://www.redhat.com)/[CentOS](https://www.centos.org/), and
many other Linux distributions.
Unlike most package managers, using dreckly does not require admin privileges.  You can install a dreckly
tree in any directory to which you have write access and easily install any
of the nearly 20,000 packages in the collection.

The
[auto-dreckly-setup](https://github.com/outpaddling/auto-admin/blob/master/User-scripts/auto-dreckly-setup)
script will help you install dreckly in about 10 minutes.  Just download it
and run

```
sh auto-dreckly-setup
```

Then, assuming you selected current packages and the default prefix

```
source ~/Dreckly/pkg/etc/dreckly.sh   # Or dreckly.csh for csh or tcsh
cd ~/Dreckly/dreckly/biology/ad2vcf
sbmake install clean clean-depends
```

### Other package managers

Packages for libxtend are known to exist in the following package managers.
These are maintained by third parties and not directly supported here.

[BioArchLinux](https://github.com/BioArchLinux/Packages)

## Instructions for packagers

### Prerequisites

* **[biolibc](https://github.com/auerlab/biolibc)**
* **[libxtend](https://github.com/outpaddling/libxtend)**

If you would like to add this project to another package manager
rather than use FreeBSD ports or dreckly, basic manual build instructions
for package can be found
[here](https://github.com/outpaddling/Coding-Standards/blob/main/package.md).
Your contribution is greatly appreciated!
