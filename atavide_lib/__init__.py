import os
import sys

__author__ = 'Rob Edwards'
from .colours import colours, colors, message
from .sequences import stream_fastq
from .read_definitions import read_definitions

__all__ = [
    'colours', 'colors', 'message', 'stream_fastq', 'read_definitions'
    ]
