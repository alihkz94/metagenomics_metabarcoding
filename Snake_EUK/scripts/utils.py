#!/usr/bin/env python3
"""
Utility module for robust, memory-efficient FASTA file handling.
This module defines a file wrapper class that streams a FASTA file,
reading it with 'latin-1' decoding and converting each line to ASCII on the fly.
The wrapper implements the full file-like interface (read, readline, seek, tell, etc.)
so that BioPython’s SeqIO can use it without error.
"""

class AsciiFileWrapper:
    def __init__(self, filename):
        # Open file using latin-1 encoding with error replacement.
        self.handle = open(filename, 'r', encoding='latin-1', errors='replace')
    
    def read(self, size=-1):
        # Read raw data and convert to ASCII.
        data = self.handle.read(size)
        return data.encode('ascii', errors='replace').decode('ascii')
    
    def readline(self):
        # Read a single line and convert to ASCII.
        line = self.handle.readline()
        return line.encode('ascii', errors='replace').decode('ascii')
    
    def __iter__(self):
        return self
    
    def __next__(self):
        line = self.readline()
        if line == "":
            raise StopIteration
        return line
    
    def seek(self, pos, whence=0):
        return self.handle.seek(pos, whence)
    
    def tell(self):
        return self.handle.tell()
    
    def close(self):
        return self.handle.close()
    
    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

def get_ascii_fasta_stream(filename):
    """
    Returns a file-like object that streams the given FASTA file,
    converting its contents from latin-1 to ASCII on the fly.
    """
    return AsciiFileWrapper(filename)
