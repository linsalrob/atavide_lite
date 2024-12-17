import gzip

from .atavide_error import FastqFormatError


__author__ = 'Rob Edwards'


def stream_fastq(fqfile):
    """Read a fastq file and provide an iterable of the sequence ID, the
    full header, the sequence, and the quaity scores.

    Note that the sequence ID is the header up until the first space,
    while the header is the whole header.
    """

    if fqfile.endswith('.gz'):
        qin = gzip.open(fqfile, 'rt')
    else:
        qin = open(fqfile, 'r')

    linecounter = 0
    while True:
        header = qin.readline()
        linecounter += 1
        if not header:
            break
        if not header.startswith("@"):
            raise FastqFormatError(f"The file {fqfile} does not appear to be a four-line fastq "
                                   f"file at line {linecounter}")
        header = header.strip()
        seqidparts = header.split(' ')
        seqid = seqidparts[0]
        seqid = seqid.replace('@', '')
        seq = qin.readline().strip()
        linecounter += 1
        qualheader = qin.readline()
        if not qualheader.startswith("+"):
            raise FastqFormatError(f"The file does not appear to be a four-line fastq file at line {linecounter}")
        linecounter += 1
        qualscores = qin.readline().strip()
        linecounter += 1
        header = header.replace('@', '', 1)
        if len(qualscores) != len(seq):
            raise FastqFormatError(f"The sequence and qual scores are not the same length at line {linecounter}")
        yield seqid, header, seq, qualscores
