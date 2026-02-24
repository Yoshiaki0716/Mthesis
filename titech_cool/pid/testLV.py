import os
import sys
import socket
import time
import math
import numpy as np
from statistics import mean
from functools import reduce
from functools import partial
from operator import add
import json
import traceback
import pprint
import multiprocessing
from Kalman import KalmanSmoother
from DeviceInvoker import *

import smtplib
import datetime
from email.mime.text import MIMEText
import subprocess
import threading

def RIGOLcmd(ip, channel, RIGOLCommand):
        """Call the RIGOL helper and return cleaned stdout string.

        Args:
            ip: device IP string
            channel: channel identifier (not used by helper but kept for signature parity)
            RIGOLCommand: method call string, e.g. 'measureVoltage()'
        """
        RIGOLDir = '/nas/dcs/wienermpod_ivi'
        cmd = f"cd {RIGOLDir}; python3 -c \"from RIGOL_DP821 import RIGOL_DP821_PYVISA; print(RIGOL_DP821_PYVISA().{RIGOLCommand})\" {ip}"

        # run the helper and grab stdout
        proc = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        out = proc.stdout or ''
        # return stripped stdout (single-line). If multiple lines, return the last non-empty line
        lines = [l.strip() for l in out.splitlines() if l.strip()]
        return lines[-1] if lines else out.strip()

def main():
    ip = "192.168.10.92"   
    ch = "2"
    try:
        # RIGOLcmd expects (ip, channel, command)
        vol = float(RIGOLcmd(ip, ch, "measureVoltage()"))
        curr = float(RIGOLcmd(ip, ch, "measureCurrent()"))
        if vol > 1.0 and curr > 1.0:
            is_on = True
        else:
            is_on = False
    except Exception as e:
        # デバッグのため例外情報を出す（本番ではログ化）
        import traceback
        print("error:", e)
        traceback.print_exc()
        # ここは設計次第。判定不能なら None にするのが安全
        is_on = None

    print("is_on:", is_on)
    return is_on

if __name__ == '__main__':
    main()