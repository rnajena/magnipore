# author: Jannes Spangenberg
# e-mail: jannes.spangenberg@uni-jena.de
# github: https://github.com/JannesSP
# website: https://jannessp.github.io
# Logger repository: https://github.com/JannesSP/Logger

import datetime
import sys
from io import TextIOWrapper
from src.Helper import ANSI
import psutil

def get_memory_usage():
    """Returns current memory usage in MB."""
    process = psutil.Process()
    return process.memory_info().rss / (1024 * 1024)  # Convert bytes to MB

class Logger():
    '''
    Logger class to print and write logs to stdout, stderr and a logfile.
    '''
    
    def __init__(self, logfilepointer: TextIOWrapper = None):
        self.lp = logfilepointer

    def _writeLog(self, string):
        '''
        Write string to logfile if logfilepointer is set.
        
        @param string: String to write to logfile.
        '''
        if self.lp is not None:
            self.lp.write(string)

    def error(self, string, error_type : str = '1'):
        '''
        Write error string to stderr and logfile if logfilepointer is set.
        Exit program with given error code.
        
        @param string: String to write to stderr and logfile.
        @param error_type: used error_type, default 1.
        '''
        sys.stderr.write(f'{datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")}, MEM: {get_memory_usage():.2f} MB, {ANSI.RED}ERROR: {string}\nMagnipore Error Code: {error_type}\n{ANSI.END}\n')
        self._writeLog(f'{datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")}, MEM: {get_memory_usage():.2f} MB, ERROR: {string}\nMagnipore Error Code: {error_type}\n')
        sys.exit(error_type)

    def warning(self, string):
        '''
        Write warning string to stderr and logfile if logfilepointer is set.
        
        @param string: String to write to stderr and logfile.
        '''
        sys.stderr.write(f'{datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")}, MEM: {get_memory_usage():.2f} MB, {ANSI.RED}WARNING: {string}{ANSI.END}\n')
        self._writeLog(f'{datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")}, MEM: {get_memory_usage():.2f} MB, WARNING: {string}\n')

    def printLog(self, string, newline_before=False, newline_after=True):
        '''
        Write datetime and string to stdout and logfile if logfilepointer is set.
        
        @param string: String to write to stdout and logfile.
        @param newline_before: Add newline before string, default False.
        @param newline_after: Add newline after string, default True.
        '''
        if newline_before:
            sys.stdout.write('\n')
            self._writeLog('\n')
            
        sys.stdout.write(f'{datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")}, MEM: {get_memory_usage():.2f} MB, LOG: {string}')
        self._writeLog(f'{datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")}, MEM: {get_memory_usage():.2f} MB, LOG: {string}')
        
        if newline_after:
            sys.stdout.write('\n')
            self._writeLog('\n')