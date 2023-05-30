"""
signet snakefile (C) 2023, Joel Bader and Zeinab Mousavi, Johns Hopkins University
"""

import os

"""
Everything is underneath this directory for now
Snakemake complains if a file path starts with ./
"""
rootdir = '.'
progname = 'signet.py'
dummyfile = 'dummy.txt'
if (rootdir != '.'):
    progname = os.path.join(rootdir, progname)
    dummyfile = os.path.join(rootdir, dummyfile)

rule all:
    input:
        dummyfile

rule run_signet:
    """
    Read raw database dump with column-sensitive fields defined by layout
    Write tab-delimited text and xlsx
    """
    input:
        progname
    output:
        dummyfile
    shell:
        """
        /usr/bin/env python {progname} -s0 0 -s1 2 -d 125000
        """
