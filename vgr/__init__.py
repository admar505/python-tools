#!/usr/bin/env python
"""
A VCFv4.0 and 4.1 parser for Python.

Online version of PyVCF documentation is available at http://pyvcf.rtfd.org/
"""


from vgr.parser import Reader, Writer
from vgr.parser import VGReader, VGWriter
from vgr.filters import Base as Filter
from vgr.parser import RESERVED_INFO, RESERVED_FORMAT
from vgr.sample_filter import SampleFilter

VERSION = '0.6.8'
