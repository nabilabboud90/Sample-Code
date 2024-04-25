import os
import glob
import sys
import shutil
import subprocess
from   subprocess import call
import time

def RunTest( exe, compare_exe, current_dir, np, keep_log = False, keep_results = True, tolerance = '1e-6' ):
    # Go to the directory specified by "current_dir"
    os.chdir(current_dir)
    
    number_of_tests = 0;
    failed = 0;

    # Save the path to some files
    datafile_dir = os.path.join(current_dir, 'datafile')
    results_dir  = os.path.join(current_dir, 'results' )
    output_dir = os.path.join(current_dir, 'output' )
    
    # Print the name of the test
    test_name = os.path.basename( current_dir )
    print( "Running " + test_name + ": " )

    # Delete any files with extension "log"
    files_to_remove=glob.glob('*.log')
    for fl in files_to_remove:
        os.unlink(fl)

    # Delete any files with extension "exo"
    files_to_remove=glob.glob('*.exo*')
    for fl in files_to_remove:
        os.unlink(fl)

    # Delete any files with extension "*~"
    files_to_remove=glob.glob('*~')
    for fl in files_to_remove:
        os.unlink(fl)
        
    # loop on the datafiles
    for datafile in glob.glob( os.path.join( datafile_dir, '*' ) ):
        # Get the name of the datafile without path
        datafile_name = os.path.basename( datafile )
        screenoutput =  ", run with the datafile:  " + datafile_name
    
        # Define output file where to store the solution and the corresponding log file
        out_file = datafile_name + "-np" + np + ".log"
        out_errors = datafile_name + "-np" + np + ".errors"
        log_file = open( out_file, 'a' )
        
        # Define the flags for the code
        flags = [datafile]
        
        # Run the code with mpirun
        code = [ "mpirun", "-np", np, exe ] + flags
        log_file.flush()
            
        # Run code
        start = time.time()
        p = subprocess.Popen(code,
                             stdout=log_file,
                             stderr=log_file)
        end = time.time()
        p.wait()
        log_file.flush()
        number_of_tests += 1

        pc = p.communicate()[0]
        hasSimFailed = 0
        if p.returncode != 1:
            hasSimFailed = 1
            
        # Get all the files in the results directory
        if hasSimFailed == 1:
            hasFailed = 1
        else:
            # Compare the new results with the reference results
            for name  in os.listdir(results_dir):
                if name == "PROCESSORRANK" or name == "NODEGID" or name == "ELEMENTGID":
                    continue
                
                results_file_dir = results_dir + "/" + name
                os.chdir(results_file_dir)
                results_file = glob.glob( '*.pvd' )
                results_file = str( results_file ).strip().split("['")[1].strip("]'")
                
                new_results_file_dir = current_dir + "/output/" + name
                os.chdir(new_results_file_dir)
                new_results_file = glob.glob( '*.pvd' )
                new_results_file = str( new_results_file ).strip().split("['")[1].strip("]'")
                
                os.chdir(current_dir)
                
                r = subprocess.Popen([compare_exe, results_file_dir, results_file, new_results_file_dir, new_results_file],
                                     stdout=log_file,
                                     stderr=log_file)
                
                rc = r.communicate()[0]
                if r.returncode == 1:
                    hasFailed = 0
                else:
                    hasFailed = 1
                    break
                
                r.wait()    
                log_file.flush()

            if hasFailed == 0:
                print("\t\tPASSED" + screenoutput + ", ( " + str( end - start ) + " s )")
            else:
                print("\t\tFAILED"  + screenoutput)
                failed += 1
            
        if hasSimFailed == 1:
            print("\t\tFAILED"  + screenoutput + ", (status = " + str( p.returncode ) + ")")
            failed += 1
        
        if not keep_log:
            files_to_remove=glob.glob('*.log')
            for fl in files_to_remove:
                os.unlink(fl)
                
        if not keep_results:
            shutil.rmtree(output_dir)
                    
    return ( number_of_tests, failed );
