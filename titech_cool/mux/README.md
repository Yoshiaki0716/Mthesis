# configMonitor.py定値を調べてinfluxDBに入れるスクリプト
過去のスキャンを遡り,その時点でのMUXの設定値を調べてinfluxDBに入れるスクリプト.



# MUX-reader.py
MUXの値を変えながら測定するスクリプト。

設定ファイル内で(MUX-reader.json),yarr のディレクトリ, connectivity config, controller configを設定してください. 調べたいMUXの設定値のリストはmuxListです。各設定値の意味はhttps://cds.cern.ch/record/2665301?ln=ja の表27: Voltage multiplexer (V_mux) assignments for ATLAS chip.(P110)を参照のこと. 
設定値はリアルタイムにinfluxDBに保存されるほか,結果がjsonファイルとして"MUX result.json"に保存されます.

