import os
import sys

__author__ = 'Rob Edwards'
from .colours import colours, colors, message
from .sequences import stream_fastq, stream_fasta
from .read_definitions import read_definitions

__all__ = [
    'colours', 'colors', 'message', 'stream_fastq', 'stream_fasta', 'read_definitions'
    ]
