* support to run on Raspberry-pi Model 3 B+, with Raspbian version 10.3 or later
* python 2.7.16 or later, 
* need to install the following packages via `pip install`

```bash
pip install numpy pyserial influxdb
```


## How to run

```bash
cd titech_cool/pid
source setup.sh # setup PYTHONPATH to lib/
./bin/pidControl.py config/TakasagoKX.json test.log
```
