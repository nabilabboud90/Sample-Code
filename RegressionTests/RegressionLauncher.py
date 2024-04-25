#!/usr/bin/env python

import os
import glob
import shutil
import sys
import subprocess
from   subprocess import call
import argparse
from RunTest import RunTest

class Regression_test:
     def __init__( self, args ):
          self.exe = args.exe
          self.compare_exe = args.compare_exe
          self.keep_log = args.keep_log
          self.keep_output = args.keep_output
          self.tolerance = str( args.tolerance )
          self.np = args.num_procs
          self.pwd = os.getcwd()
     
     def run_tests(self):
          # Initialize counters
          self.total_tests = 0
          self.failed = 0
          self.failed_tests = []
          
          # Loop over all subfolders
          for root, subdirs, files in os.walk(self.pwd):
               # Check if datafile is in the subfolder
               val = [i for i, x in enumerate(subdirs) if x == "datafile"]
               if len(val) > 0:
                    # Go to the subfolder (typically a directory corresponding to a specific test)
                    os.chdir( root )
                    
                    # Execute test
                    num_tests, num_failed = RunTest(self.exe, self.compare_exe, root, self.np, self.keep_log, self.keep_output, self.tolerance)
                    self.total_tests += num_tests
                    self.failed += num_failed
                    if( num_failed > 0 ):
                         self.failed_tests.append(root)

     def output_summary(self):
          print( "Summary: failed " + str( self.failed ) + " tests out of " + str( self.total_tests ) )
          if( self.failed > 0 ):
               print( "List of test failed:" )
               for test in self.failed_tests:
                    print( "  - " + test )
                    

def main():
     # Build parser
     parser = argparse.ArgumentParser()
     parser.add_argument( "exe", help = "path to the executable" )
     parser.add_argument( "compare_exe", help = "path to the executable that compares the vtu files" )
     parser.add_argument( "-klog", "--keep_log", help = "Keep the screen output of the simulation. It will be stored in a .log file", default = False )
     parser.add_argument( "-kout", "--keep_output", help = "Keep the the output of the simulations. It will be stored in a .pvd file", default = False )
     parser.add_argument( "-tol", "--tolerance", help = "Tolerance of the difference between solutions", default = 1.0e-10 )
     parser.add_argument( "-np", "--num_procs", help = "Number of processors", default = '2' )
 
     # Parse
     args = parser.parse_args()

     # Create, run and output regression test
     regression_test = Regression_test(args)
     regression_test.run_tests()
     regression_test.output_summary()

if __name__ == "__main__": main()        
