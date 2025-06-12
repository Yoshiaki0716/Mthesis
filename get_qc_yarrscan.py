#!/usr/bin/env python3

import pprint
import json
import threading
import sys
import os
import io
import re
import datetime
import subprocess
from optparse import OptionParser
from pathlib import Path
import logging
import glob

sys.path.append(os.environ['ITK_QC_ANALYSIS_DIR'] + '/localdb')
path = "/nas/daq/QCtest/data/repicdaq*"

from localdb_api import ldb_api

logging.basicConfig(format='%(asctime)s - %(name)s - %(levelname)s : %(message)s', level=logging.INFO)
logger = logging.getLogger(__name__)

# Set output directory
temp_zip_dir = 'temp_unzip_directory__to_remove'
outdir_zip = 'temp_zip_download_directory__to_remove'

if __name__ == "__main__":

    parser = OptionParser()
    parser.add_option("--host", dest="hostname", default="192.168.100.104",
                      help="LocalDB Hostname (default: 192.168.100.104)")
    parser.add_option("--port", dest="port", default=27017,
                      help="LocalDB Port Number (default: 27017)")
    parser.add_option("--serialNumber", dest="module_serial_number",
                      help="Module Serial Number (e.g. 20UPGM22601020)")
    parser.add_option("--stage", dest="stage",
                      help="Stage (e.g. MODULE/INITIAL_WARM)")
    parser.add_option("--test", dest="test",
                      help="test (MIN_HEALTH_TEST or PIXEL_FAILURE_ANALYSIS or TUNING)")

    (options, args) = parser.parse_args()

    # Check if required options are given.
    if not options.module_serial_number:
        print("Error: --serialNumber option is required.")
        exit(-1)
    if not options.stage:
        print("Error: --stage option is required.")
        exit(-1)

    if not options.test:
        print("Error: --test option is required.")
        exit(-1)

    if not options.test in ["MIN_HEALTH_TEST", "PIXEL_FAILURE_ANALYSIS", "TUNING"]:
        print("Error:invalid test name. " + options.test )
        exit(-1)
        
    QC_results = {}

    if not os.path.exists(outdir_zip):
        os.system('mkdir -p ' + outdir_zip)
    if not os.path.exists(temp_zip_dir):
        os.system('mkdir -p ' + temp_zip_dir)
    
    ldb = ldb_api( host = options.hostname, port = options.port )

    # Get the list of children FE chip serial numbers
    fe_serial_numbers = ldb.get_child_FEs( options.module_serial_number )
    if fe_serial_numbers == None:
        print("No FEs are found.")
        exit(0)
    
    FElist = "--FElist"
    for fe_serial_number in fe_serial_numbers:
        fe_serial_number_int = int(fe_serial_number[7:])
        fe_serial_number_hex = hex(fe_serial_number_int)
        FElist += f" {fe_serial_number_hex}"

    # Get the corresponding QC result document ID
    fe_test_run_ids = [ ldb.get_qc_status( fe_serial_number, options.stage, options.test ) for fe_serial_number in fe_serial_numbers ]
    
    # Get the corresponding QC result document, limiting the info to the binary attachment file IDs
    fe_test_run_docs = [ ldb.get_test_run_doc( test_run_id, vars = ["gridfs_attachments"] ) for test_run_id in fe_test_run_ids ]

    output_log_fs_data = []
    zip_names = []
    for fe_test_run_doc in fe_test_run_docs:
        if fe_test_run_doc == None:
            print("No log file found for this FE!")
            continue
        for key in fe_test_run_doc["gridfs_attachments"].keys():
            print("attachment: " + key)
            if key.endswith('.zip'):
                output_log_fs_data.append(ldb.get_gridfs_data(fe_test_run_doc.get("gridfs_attachments").get(key)))
                zip_names.append(key)
        #break

    if len(zip_names) == 0:
        print("No log file was found for this module; probably tests are not registered yet.")
        exit(-1)
        
    # Dump to local storage
    runNumbers = {}
    os.system(f"rm -rfv {temp_zip_dir}/*")
    os.system(f"rm -rfv {outdir_zip}/*")
    
    for index, (fe_serial_number, data) in enumerate( zip( fe_serial_numbers, output_log_fs_data ) ):

        if not data:
            logger.warning( f'{options.module_serial_number} FE{index+1} output.log was not possible to acquire' )
            continuec

        output_path = Path( f'{outdir_zip}/{zip_names[index]}' )
        output_path.parent.mkdir( parents = True, exist_ok = True )

        logger.info( f'Dumped: {output_path}' )
        
        with output_path.open('wb') as f:
            f.write( data )

        os.system("unzip " + outdir_zip + "/" + zip_names[index] + " -d " + temp_zip_dir)
        dirnames = glob.glob(temp_zip_dir + "/std_*")

        for dirname in dirnames:
            with open(dirname + "/testRun.json") as f:
                data = json.load(f)
                runNumbers[dirname.split('/')[-1]] = data['runNumber']

        fe_serial_number_int = int(fe_serial_number[7:])
        fe_serial_number_hex = hex(fe_serial_number_int)

        os.system(f"mv -v {temp_zip_dir}/scans.json {temp_zip_dir}/scans_{fe_serial_number_hex}.json")
        os.system(f"mv -v {temp_zip_dir}/info.json {temp_zip_dir}/info_{fe_serial_number_hex}.json")

        # Rename duplicated file names
        scan_names_tmp = subprocess.Popen(f"ls -1d {temp_zip_dir}/std*", stdout = subprocess.PIPE, text = True, shell = True).stdout
        scan_names = []
        for scan_name_tmp in scan_names_tmp:
            scan_name = scan_name_tmp.split('/')[-1].replace('¥n', '')
            os.system(f"mv -v {temp_zip_dir}/{scan_name}/testRun.json {temp_zip_dir}/{scan_name}/{fe_serial_number_hex}_testRun.json")
            os.system(f"mv -v {temp_zip_dir}/{scan_name}/{scan_name}.json {temp_zip_dir}/{scan_name}/{fe_serial_number_hex}_{scan_name}.json")
            os.system(f"mv -v {temp_zip_dir}/{scan_name}/ctrlCfg.json {temp_zip_dir}/{scan_name}/{fe_serial_number_hex}_ctrlCfg.json")
            os.system(f"mv -v {temp_zip_dir}/{scan_name}/componentTestRun.json {temp_zip_dir}/{scan_name}/{fe_serial_number_hex}_componentTestRun.json")

            if scan_name in runNumbers.keys():
                if os.path.exists(f"{temp_zip_dir}/{runNumbers[scan_name]:06}_{scan_name}"):
                    os.system(f"mv -v {temp_zip_dir}/{scan_name}/* {temp_zip_dir}/{runNumbers[scan_name]:06}_{scan_name}")
                    os.system(f"rm -rfv {temp_zip_dir}/{scan_name}")
                else:
                    os.system(f"mv -v {temp_zip_dir}/{scan_name} {temp_zip_dir}/{runNumbers[scan_name]:06}_{scan_name}")

        os.system(f"mv -v {temp_zip_dir}/output {temp_zip_dir}/output_{fe_serial_number_hex}")
            
    #os.system("rm -rfv " + outdir_zip + "/*")
    #os.system("rm -rfv " + temp_zip_dir + "/*")

    # Convert json to ROOT
    scan_paths = glob.glob(f"{temp_zip_dir}/*_std_*")
    for scan_path in scan_paths:
        command = f"python3 ../analysis/YARRscan_hist_saver.py -d {scan_path} -s {options.module_serial_number} {FElist} -y"
        os.system(command)

    stage_short = options.stage.split('/')[-1]
    final_data_dir = f'../../outputs/{options.module_serial_number}/{stage_short}'

    if not os.path.exists(final_data_dir):
        os.system(f"mkdir -p {final_data_dir}")
        
    os.system(f"mv -v ../../outputs/*.root {final_data_dir}")

    #if not os.path.exists(final_data_dir):
    #    os.system(f"mkdir -p {final_data_dir}")
    #temp_unzip_directory__to_remove/output_0x21253/PIXEL_FAILURE_ANALYSIS/2025-02-13_131836/output.log
        
    lognames = glob.glob(f"{temp_zip_dir}/output_*/*/*/output.log")
    for logname in lognames:
        fe_name = logname.split('/')[-4]
        os.system(f"mv -v {logname} {final_data_dir}/{fe_name}.log")
    
    """
    print("")
    print("Scan numbers for " + options.test + ": " + options.module_serial_number + " (" + options.stage + ")")
    print("-----------------------------------")
    for key in runNumbers.keys():
        print(key.ljust(26) + " | " + "%5d"%runNumbers[key])
    print("-----------------------------------")
    print("")

    path += "/" + options.module_serial_number
    #if options.stage.endswith('COLD'):
    #    path += "/cold"
    #else:
    #    path += "/warm"
    path += "/*"
    for key in runNumbers:
        rn = "%06d"%runNumbers[key]
        test = glob.glob(path + "/YARRscans/" + rn + "*")
        print(test[0])
    
    """
exit(0