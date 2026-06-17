#!/usr/bin/env python3
import sys
import time
import serial

# ==========================================
# ⚙️ 環境設定（ご自身の環境に合わせて修正してください）
# ==========================================
# Keithley 2400 が繋がっているシリアルポートのパス
# ※ ls /dev/serial/by-id/ 等で調べて書き換えてください
USB_PORT = '/dev/ttyUSB0' 
BAUD_RATE = 9600

# 高電圧(HV)の設定値
TARGET_HV = -113.0  # ONにしたときの最終目標電圧(V)
STEP_V    = 1.0     # 1回に変化させる電圧幅(V)
DELAY_S   = 0.1     # ステップ間の待機時間(秒)
# ==========================================

def send_cmd(ser, cmd, expect_response=False):
    """KeithleyへSCPIコマンドを送信し、必要なら返答を読み取る"""
    ser.write(f"{cmd}\n".encode())
    if expect_response:
        return ser.readline().decode('utf-8').strip()
    return None

def read_voltage(ser):
    """現在の出力電圧を読み取る"""
    # Keithley 2400 の :MEAS:VOLT? は "電圧,電流,抵抗,時間,ステータス" のカンマ区切りで返る
    resp = send_cmd(ser, ":MEAS:VOLT?", expect_response=True)
    if resp:
        try:
            # 最初の要素が電圧値
            vol_str = resp.split(',')[0]
            return float(vol_str)
        except Exception as e:
            print(f"読み取りエラー: {resp} ({e})", file=sys.stderr)
            return None
    return None

def set_voltage(ser, voltage):
    """電圧をセットする"""
    send_cmd(ser, f":SOUR:VOLT:LEV {voltage}")

def ramp_voltage(ser, target_voltage):
    """安全に電圧を目標値までランプ（徐々に変化）させる"""
    current_voltage = read_voltage(ser)
    if current_voltage is None:
        print("現在の電圧が読み取れません。安全のため処理を中断します。", file=sys.stderr)
        sys.exit(1)

    print(f"[Keithley 2400] 現在の電圧: {current_voltage:.3f} V -> 目標電圧: {target_voltage:.3f} V")
    
    # もし現在出力がOFFなら、いきなり電圧をかけないように0VにセットしてからONにする
    outp_state = send_cmd(ser, ":OUTP?", expect_response=True)
    if outp_state == "0":
        set_voltage(ser, 0.0)
        send_cmd(ser, ":OUTP ON")
        current_voltage = 0.0
        time.sleep(0.5)

    # 電圧を上げるか下げるかの方向を決定 (1: プラス方向へ, -1: マイナス方向へ)
    direction = 1 if target_voltage > current_voltage else -1
    v = current_voltage

    while True:
        # 目標に到達（または通り過ぎた）場合の処理
        if (direction == 1 and v >= target_voltage) or (direction == -1 and v <= target_voltage):
            v = target_voltage
            set_voltage(ser, v)
            print(f"[Keithley 2400] 到達: {v:.3f} V")
            break

        # 電圧を1ステップ進める
        v += direction * STEP_V
        
        # 行き過ぎの補正
        if (direction == 1 and v > target_voltage) or (direction == -1 and v < target_voltage):
            v = target_voltage

        set_voltage(ser, round(v, 4))
        print(f"[Keithley 2400] 設定中... {v:.3f} V")
        time.sleep(DELAY_S)

    # 目標が0V（OFFコマンド）の場合は、0Vまで下がりきった後に完全に出力を遮断する
    if target_voltage == 0.0:
        send_cmd(ser, ":OUTP OFF")
        print("出力を完全にOFFにしました。")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("使い方: python3 testHV.py [on | off | status]", file=sys.stderr)
        sys.exit(1)

    command = sys.argv[1].lower()

    try:
        # ランピング中は通信を何度も行うため、1つのセッションを開きっぱなしにする（高速化・安定化）
        with serial.Serial(USB_PORT, BAUD_RATE, timeout=1.0) as ser:
            
            # 手動操作を防ぐためリモートモードへ
            send_cmd(ser, ":SYST:REM")

            if command == "on":
                print("--- 高電圧(HV)の印加を開始します (Ramp UP) ---")
                ramp_voltage(ser, TARGET_HV)
                print("--- 完了 ---")

            elif command == "off":
                print("--- 高電圧(HV)の遮断を開始します (Ramp DOWN) ---")
                ramp_voltage(ser, 0.0)
                print("--- 完了 ---")

            elif command == "status":
                vol = read_voltage(ser)
                if vol is not None:
                    # -1.0V 〜 +1.0V の間は完全に OFF（0）とみなす
                    if abs(vol) < 1.0:
                        print("0")
                    else:
                        print(f"{vol:.3f}")
                else:
                    print("0")

            else:
                print(f"エラー: 不明なコマンド '{command}' です。", file=sys.stderr)
                # エラー時もロックを解除して終了
                send_cmd(ser, ":SYST:LOC")
                sys.exit(1)

            # 処理が終わったらフロントパネルのロックを解除（手動操作可能にする）
            send_cmd(ser, ":SYST:LOC")

    except serial.SerialException as e:
        print(f"シリアル通信エラー: {e}", file=sys.stderr)
        print("※USB_PORTのパスが間違っていないか、ケーブルが抜けていないか確認してください。", file=sys.stderr)
        sys.exit(1)