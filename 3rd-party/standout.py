# standout.py
# v3.0.0
# A flexible object to redirect standard output and standard error
# Allows logging to a file and to set a level of verbosity

# Copyright (C) 2006 Michael Foord
# E-mail: fuzzyman AT voidspace DOT org DOT uk

# StandOut 3
# http://www.voidspace.org.uk/python/standout.html

# Released subject to the BSD License
# Please see http://www.voidspace.org.uk/python/license.shtml

# Scripts maintained at http://www.voidspace.org.uk/python/index.shtml
# For information about bugfixes, updates and support, please join the
# Pythonutils mailing list:
# http://groups.google.com/group/pythonutils/
# Comments, suggestions and bug reports welcome.

import sys

__version__ = '3.0.0'
__docformat__ = "restructuredtext en"

__all__ = (
    'StandOut',
)

class StandOut(object):
    
    def __init__(self, logfile=None, logmode="w", stdout=True, stderr=True,
        errPrefix='[err] ', priority=5, threshold=None, errThreshold=None,
        outThreshold=None, errLogfileThreshold=None, outLogfileThreshold=None,
        outStreamThreshold=None, errStreamThreshold=None, unbuffered=False):
        
        self.__stdout = None
        self.__stderr = None
        self.logfile = logfile
        self._logfile = None
        self._logmode = logmode
        self.__outputStream = None
        self.__errorStream = None
        self.__priority = priority
        
        _outStreamThreshold = 5
        _errStreamThreshold = 0
        _outLogfileThreshold = 5
        _errLogfileThreshold = 0
        
        if threshold is not None:
            _outStreamThreshold = _errStreamThreshold = _outLogfileThreshold = _errLogfileThreshold = threshold
        if outThreshold is not None:
            _outStreamThreshold = _outLogfileThreshold = outThreshold
        if errThreshold is not None:
            errStreamThreshold = errLogfileThreshold = errThreshold
        if outStreamThreshold is not None:
            _outStreamThreshold = outStreamThreshold
        if errStreamThreshold is not None:
            _errStreamThreshold = errStreamThreshold
        if outLogfileThreshold is not None:
            _outLogfileThreshold = outLogfileThreshold
        if errLogfileThreshold is not None:
            _errLogfileThreshold = errLogfileThreshold
        
        if logfile:
            self._logfile = open(self.logfile, self._logmode)
        
        if stdout:
            self.__stdout = sys.stdout
            self.__outputStream = _Stream(sys.stdout, self._logfile,
                priority=priority, streamThreshold=_outStreamThreshold,
                logfileThreshold=_outLogfileThreshold, unbuffered=unbuffered)
            sys.stdout = self.__outputStream
        if stderr:
            self.__errorStream = _Stream(sys.stderr, self._logfile,
                prefix=errPrefix, priority=priority,
                streamThreshold=_errStreamThreshold,
                logfileThreshold=_errLogfileThreshold, unbuffered=unbuffered)
            self.__stderr = sys.stderr
            sys.stderr = self.__errorStream


    def __setPriority(self, value):
        self.__priority = value
        if self.__errorStream:
            self.__errorStream._priority = value
        if self.__outputStream:
            self.__outputStream._priority = value
            
    priority = property(lambda self: self.__priority, __setPriority)
            
            
    def __getThreshold(self):
        if self.errThreshold == self.outThreshold:
            return self.errThreshold
        return -1
    
    def __setThreshold(self, value):
        if self.__errorStream:
            self.__errorStream._streamThreshold = value
            self.__errorStream._logfileThreshold = value
        if self.__outputStream:
            self.__outputStream._streamThreshold = value
            self.__outputStream._logfileThreshold = value
    
    threshold = property(__getThreshold, __setThreshold)


    def __setErrThreshold(self, value):
        if self.__errorStream:
            self.__errorStream._streamThreshold = value
            self.__errorStream._logfileThreshold = value

    def __getErrThreshold(self):
        if self.__errorStream and (self.__errorStream._streamThreshold
                                   == self.__errorStream._logfileThreshold):
            return self.__errorStream._streamThreshold
        return -1

    errThreshold = property(__getErrThreshold, __setErrThreshold)


    def __setOutThreshold(self, value):
        if self.__outputStream:
            self.__outputStream._streamThreshold = value
            self.__outputStream._logfileThreshold = value

    def __getOutThreshold(self):
        if self.__outputStream and (self.__outputStream._streamThreshold
                                   == self.__outputStream._logfileThreshold):
            return self.__outputStream._streamThreshold
        return -1

    outThreshold = property(__getOutThreshold, __setOutThreshold)


    def __setErrStreamThreshold(self, value):
        if self.__errorStream:
            self.__errorStream._streamThreshold = value

    def __getErrStreamThreshold(self):
        if self.__errorStream:
            return self.__errorStream._streamThreshold
        return -1

    errStreamThreshold = property(__getErrStreamThreshold, __setErrStreamThreshold)


    def __setErrLogfileThreshold(self, value):
        if self.__errorStream:
            self.__errorStream._logfileThreshold = value

    def __getErrLogfileThreshold(self):
        if self._logfile and self.__errorStream:
            return self.__errorStream._logfileThreshold
        return -1

    errLogfileThreshold = property(__getErrLogfileThreshold, __setErrLogfileThreshold)


    def __setOutStreamThreshold(self, value):
        if self.__outputStream:
            self.__outputStream._streamThreshold = value

    def __getOutStreamThreshold(self):
        if self.__outputStream:
            return self.__outputStream._streamThreshold
        return -1

    outStreamThreshold = property(__getOutStreamThreshold, __setOutStreamThreshold)


    def __setOutLogfileThreshold(self, value):
        if self.__outputStream:
            self.__outputStream._logfileThreshold = value

    def __getOutLogfileThreshold(self):
        if self._logfile and self.__outputStream:
            return self.__outputStream._logfileThreshold
        return -1

    outLogfileThreshold = property(__getOutLogfileThreshold, __setOutLogfileThreshold)
    
    
    def close(self):
        if self.__stdout:
            sys.stdout = self.__stdout
        if self.__stderr:
            sys.stderr = self.__stderr
        if self._logfile:
            self._logfile.close()


class _Stream(object):

    def __init__(self, stream, outfile=None, prefix=None, priority=5,
        streamThreshold=5, logfileThreshold=5, unbuffered=False):
        self._stream = stream
        self._outfile = outfile
        self._prefix = prefix
        self._done_newline = True
        self._priority = priority
        self._streamThreshold = streamThreshold
        self._logfileThreshold = logfileThreshold
        self._unbuffered = unbuffered

    def write(self, data, priority=None):
        if not data:
            # avoid writing the prefix for null input
            return
        if priority is None:
            priority = self._priority
        
        if  self._prefix is not None:
            # Need every newline (including the first) to start with the prefix
            # but not print the prefix for a newline if it is the last character
            # instead do it on the next print
            if self._done_newline:
                data = self._prefix + data
            self._done_newline = False
            terminated = False
            if data[-1] == '\n':
                terminated = True
                data = data[:-1]
                self._done_newline = True
            data = data.replace('\n', '\n' + self._prefix)
            if terminated:
                data = data + '\n'
        
        if priority >= self._streamThreshold:
            self._stream.write(data)
            if self._unbuffered:
                self._stream.flush()
        if self._outfile is not None and priority >= self._logfileThreshold:
            self._outfile.write(data)

    def writelines(self, inLines):
        self.write(''.join(inLines))

    def __getattr__(self, attribute):
        # doesn't proxy attributes shread by this class and the underlying stream
        # this affects next and some double underscore attributes
        try:
            return self.__dict__[attribute]
        except KeyError:
            return getattr(self._stream, attribute)
