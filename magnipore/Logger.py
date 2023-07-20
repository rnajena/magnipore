# author: Jannes Spangenberg
# e-mail: jannes.spangenberg@uni-jena.de
# github: https://github.com/JannesSP
# website: https://jannessp.github.io
# Logger repository: https://github.com/JannesSP/Logger

import datetime
import sys
from io import TextIOWrapper
from magnipore.Helper import ANSI

class Logger():
    '''
    Logger class to print and write logs to stdout, stderr and a logfile.
    '''
    
    def __init__(self, logfilepointer: TextIOWrapper = None):
        self.lp = logfilepointer

    def writeLog(self, string):
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
        sys.stderr.write(f'{ANSI.RED}ERROR: {string}\nMagnipore Error Code: {error_type}\n{ANSI.END}\n')
        self.writeLog(f'ERROR: {string}\nMagnipore Error Code: {error_type}\n')
        sys.exit(error_type)

    def warning(self, string):
        '''
        Write warning string to stderr and logfile if logfilepointer is set.
        
        @param string: String to write to stderr and logfile.
        '''
        sys.stderr.write(f'{ANSI.RED}WARNING: {string}{ANSI.END}\n')
        self.writeLog(f'WARNING: {string}\n')

    def printLog(self, string, newline_before=False, newline_after=True):
        '''
        Write datetime and string to stdout and logfile if logfilepointer is set.
        
        @param string: String to write to stdout and logfile.
        @param newline_before: Add newline before string, default False.
        @param newline_after: Add newline after string, default True.
        '''
        if newline_before:
            sys.stdout.write('\n')
            self.writeLog('\n')
            
        sys.stdout.write(f'{datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")} LOG: {string}')
        self.writeLog(f'{datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")} LOG: {string}')
        
        if newline_after:
            sys.stdout.write('\n')
            self.writeLog('\n')