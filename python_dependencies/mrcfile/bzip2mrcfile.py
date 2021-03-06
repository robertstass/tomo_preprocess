# Copyright (c) 2016, Science and Technology Facilities Council
# This software is distributed under a BSD licence. See LICENSE.txt.
"""
bzip2mrcfile
------------

Module which exports the :class:`Bzip2MrcFile` class.

Classes:
    :class:`Bzip2MrcFile`: An object which represents a bzip2 MRC file.

"""

# Import Python 3 features for future-proofing
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)


import bz2
import os

from .mrcfile import MrcFile


class Bzip2MrcFile(MrcFile):
    
    """:class:`~mrcfile.mrcfile.MrcFile` subclass for handling bzip2 files.
    
    Usage is the same as for :class:`~mrcfile.mrcfile.MrcFile`.
    
    """
    
    def __repr__(self):
        return "Bzip2MrcFile('{0}', mode='{1}')".format(self._fname,
                                                        self._mode)
    
    def _open_file(self, name):
        """Override _open_file() to open a bzip2 file."""
        self._fname = name
        if 'w' in self._mode and not os.path.exists(name):
            open(name, mode='w').close()
        self._iostream = bz2.BZ2File(name, mode='r')
    
    def _read(self):
        """Override _read() to ensure bzip2 file is in read mode."""
        self._ensure_readable_stream()
        super(Bzip2MrcFile, self)._read()
    
    def _ensure_readable_stream(self):
        """Make sure _iostream is a bzip2 stream that can be read."""
        self._iostream.close()
        self._iostream = bz2.BZ2File(self._fname, mode='r')
    
    def _get_file_size(self):
        """Override _get_file_size() to avoid seeking from end."""
        self._ensure_readable_stream()
        return super(Bzip2MrcFile, self)._get_file_size()
    
    def flush(self):
        """Override flush() since BZ2File objects need special handling."""
        if not self._read_only:
            self._iostream.close()
            self._iostream = bz2.BZ2File(self._fname, mode='w')
            
            # Arrays converted to bytes so gzip can calculate sizes correctly
            self._iostream.write(self.header.tobytes())
            self._iostream.write(self.extended_header.tobytes())
            self._iostream.write(self.data.tobytes())
            # no equivalent for flush() with BZ2File
