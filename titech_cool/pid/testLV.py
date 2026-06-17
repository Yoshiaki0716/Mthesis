#!/usr/bin/env python3
import sys
import time
import serial

# ==========================================
# ⚙️ 環境設定（ご自身の環境に合わせて修正してください）
# ==========================================
# TEXIO PFRの絶対パス（/dev/serial/by-id/ 等の固定パスを推奨）
USB_PORT = '/dev/LV'
BAUD_RATE = 9600

# ノイズ丸め込みのしきい値（この値未満は 0 として出力）
V_TOLERANCE = 1.0  # 0.05V 未満は 0V とみなす
I_TOLERANCE = 1.0  # 0.01A (10mA) 未満は 0A とみなす
# ==========================================

def send_texio_cmd(cmd, expect_response=False):
    """TEXIOへコマンドを送信し、必要なら返答を読み取る"""
    try:
        with serial.Serial(USB_PORT, BAUD_RATE, timeout=1.0) as ser:
            time.sleep(0.05) # 接続安定化
            
            # リモートモード（PC制御）へ移行
            ser.write(b':SYSTem:COMMunicate:RLSTate REMote\n')
            time.sleep(0.05)

            # メインのコマンドを送信
            ser.write(f'{cmd}\n'.encode())
            time.sleep(0.05)

            response = None
            if expect_response:
                response = ser.readline().decode('utf-8', errors='ignore').strip()

            # ローカルモード（手動操作可能）へ戻す
            ser.write(b':SYSTem:COMMunicate:RLSTate LOCal\n')
            time.sleep(0.05)

            return response
    except Exception as e:
        # pidServerの誤判定を防ぐため、エラー時は標準エラー出力に逃がす
        print(f"Serial Error: {e}", file=sys.stderr)
        return None

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("使い方: python3 testLV.py [turnOn | turnOff | measureVoltage | measureCurrent]", file=sys.stderr)
        sys.exit(1)

    # 柔軟に対応するため、小文字にして判定
    command = sys.argv[1].lower()

    if command == "on":
        send_texio_cmd('CURR 6.25')
        send_texio_cmd('OUTP ON')
        

    elif command == "off":
        send_texio_cmd('OUTP OFF')

    elif command == "status":
        vol_str = send_texio_cmd('MEAS:VOLT?', expect_response=True)
        curr_str = send_texio_cmd('MEAS:CURR?', expect_response=True)

        if vol_str and curr_str:
            try:
                vol = float(vol_str)
                curr = float(curr_str)
                # 電圧と電流が両方 1.0 以上の時のみ "1" (ON) と判定
                if vol > 1.0 and curr > 1.0:
                    print("1")
                else:
                    print("0")
            except ValueError:
                print("0")
        else:
            print("0")

    else:
        print(f"Error: Unknown command '{sys.argv[1]}'", file=sys.stderr)
        sys.exit(1)