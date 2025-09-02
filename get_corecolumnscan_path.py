#!/usr/bin/env python3

import os
import sys
import glob
from optparse import OptionParser

sys.path.append(os.environ['ITK_QC_ANALYSIS_DIR'] + '/localdb')

from localdb_api import ldb_api

if __name__ == "__main__":

    parser = OptionParser()
    parser.add_option("--host", dest="hostname", default="192.168.100.104",
                      help="LocalDB Hostname (default: 192.168.100.104)")
    parser.add_option("--port", dest="port", default=27017,
                      help="LocalDB Port Number (default: 27017)")
    parser.add_option("--stage", dest="stage",
                      help="Stage (e.g. MODULE/INITIAL_WARM)")
    #parser.add_option("--test", dest="test",
                    #help="test (MIN_HEALTH_TEST or PIXEL_FAILURE_ANALYSIS or TUNING)")

    (options, args) = parser.parse_args()

    # Check if required options are given.

    if not options.stage:
        print("Error: --stage option is required.")
        exit(-1)
   # if not options.test:
       # print("Error: --test option is required.")
       # exit(-1)    

        
    ldb = ldb_api( host = options.hostname, port = options.port )

    #test = ldb.get_test_run_all_doc( {"testType": "corecolumnscan", "stage": options.stage} )
    #test = ldb.get_list_testRun( {"testType": "corecolumnscan", "stage": options.stage} )
    test = ldb.get_list_testRun( {"testType": "analogscan", "stage": options.stage} )

    paths = {}
    
    for t in test:
        print(t["exec"])
        run = t["runNumber"]
        serial = t["exec"].split()[-1].split("/")[-3]
        print(serial)
        #serial = "20UPGM23610244"
        repicdaq = t["exec"].split()[-1].split("/")[-4]
        print(repiqdaq)

        #path = f"/nas/daq/QCtest/data/{repicdaq}/{serial}/*/YARRscans/{run:06}_corecolumnscan"
        path = f"/nas/daq/QCtest/data/{repicdaq}/{serial}/*/YARRscans/{run:06}_std_analogscan"

        if not serial in paths.keys():
            paths[serial] = []

        paths[serial].append(path)

    for serial in paths.keys():
        path = paths[serial][-1] 
        print(path)

        try:
            dirname = glob.glob(f"{path}")[0]
        except:
            print(f"Error: corecolumnscan wasn't found for {serial}.")
            continue

        stage_short = options.stage.split("/")[-1]
        outputpath = f"../Output/{serial}/{stage_short}"

        os.system(f"mkdir -p {outputpath}")

        os.system(f"python3 ./analysis/YARRscan_hist_saver.py --serial {serial} -d {path} -o {outputpath} --yes")

        #exit(-1)

