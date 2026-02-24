import pyvisa
import time
import sys

# --- 設定項目 ---
# 前回の確認結果に基づきポートを指定
RESOURCE_NAME = "ASRL/dev/ttyUSB0::INSTR" 
TARGET_VOLTAGE = 120.0  # 目標電圧 (V)
COMPLIANCE_CURR = 0.0001 # 電流制限 (100uA)
STEP_VOLTAGE = 10.0      # 1ステップで上げる電圧 (V)
STEP_DELAY = 0.5         # ステップごとの待ち時間 (秒)

def main():
    if len(sys.argv) < 2:
        print("Usage: python3 hv_control.py [on/off]")
        return

    action = sys.argv[1].lower()
    rm = pyvisa.ResourceManager('@py')

    try:
        # 接続
        inst = rm.open_resource(RESOURCE_NAME)
        inst.baud_rate = 9600
        inst.write_termination = '\r'
        inst.read_termination = '\r'
        inst.timeout = 5000

        if action == "on":
            print(f"--- HV ON: {TARGET_VOLTAGE}V へ昇圧を開始します ---")
            inst.write("*RST")                  # 初期化
            inst.write(":SOUR:FUNC VOLT")       # 電圧源モード
            inst.write(":SENS:CURR:PROT " + str(COMPLIANCE_CURR)) # 電流制限
            inst.write(":SOUR:VOLT:RANG 200")   # 200Vレンジに固定 (120V出力のため)
            inst.write(":OUTP ON")              # 出力開始 (最初は0V)

            # 0Vから目標電圧まで段階的に上げる (Ramp-up)
            current_v = 0.0
            while current_v < TARGET_VOLTAGE:
                current_v += STEP_VOLTAGE
                if current_v > TARGET_VOLTAGE: current_v = TARGET_VOLTAGE
                inst.write(f":SOUR:VOLT {current_v}")
                print(f"Current Set: {current_v} V")
                time.sleep(STEP_DELAY)
            print("120V 到達。出力を維持します。")

        elif action == "off":
            print("--- HV OFF: 降圧して停止します ---")
            # 安全のため段階的に下げる (Ramp-down)
            # 現在の電圧を取得（簡易的に120Vから下げると想定）
            for v in range(int(TARGET_VOLTAGE), -1, -int(STEP_VOLTAGE)):
                inst.write(f":SOUR:VOLT {v}")
                time.sleep(STEP_DELAY / 2)
            
            inst.write(":SOUR:VOLT 0")
            inst.write(":OUTP OFF")
            print("出力を停止しました。")

    except Exception as e:
        print(f"エラーが発生しました: {e}")
    finally:
        if 'inst' in locals():
            inst.close()

if __name__ == "__main__":
    main()
