#!/usr/bin/env python3

import sys
import math
import time
import traceback
from datetime import datetime
from ArduinoSHT85Control import *
from influxdb import InfluxDBClient

# --- 設定 ---
USB_PORT = "/dev/ttyACM0"
BAUD_RATE = 9600
DB_HOST = 'localhost'
DB_PORT = 8086
DB_NAME = 'pidDB'
MEASUREMENT = 'EnvMon'
LOG_FILE_PATH = '/home/admin/titech_cool/pid/test_REPIC_envlog.txt'

try:
    # 初期化
    sht85 = ArduinoSHT85Control({"port": USB_PORT, "baud": BAUD_RATE})
    client = InfluxDBClient(host=DB_HOST, port=DB_PORT, database=DB_NAME)

    print("Waiting 3 seconds for Arduino to boot...")
    time.sleep(3)

    now = datetime.now()
    now_ts = now.timestamp()
    
    # データの取得
    t = sht85.getT()
    rh = sht85.getRH()
    
    # 露点の計算
    b = 18.678
    c = 257.14
    g = math.log(rh / 100.0) + b * t / (c + t)
    dp = c * g / (b - g)

    # 1. ログファイルへ書き込み
    log_line = f"{now_ts:.3f} {now} {t:.3f}[C] {rh:.3f}[RH%] {dp:.3f}[C_DP]\n"
    with open(LOG_FILE_PATH, 'a') as f:
        f.write(log_line)

    # 2. InfluxDBへデータ送信
    json_body = [
        {"measurement": MEASUREMENT, "tags": {"sensor": "sht85_temp"}, "fields": {"SHT85 Temp": float(t)}},
        {"measurement": MEASUREMENT, "tags": {"sensor": "sht85_humi"}, "fields": {"SHT85 RH": float(rh)}},
        {"measurement": MEASUREMENT, "tags": {"sensor": "sht85_dp"},   "fields": {"SHT85 DP": float(dp)}}
    ]
    client.write_points(json_body)
    print(f"[{now}] Successfully sent to InfluxDB: Temp={t:.1f}, RH={rh:.1f}, DP={dp:.1f}")

except Exception as e:
    print(f"[{datetime.now()}] Error: {e}")
    traceback.print_exc()

# ループせずにそのまま終了する
