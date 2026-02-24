#!/bin/bash

sleep 60

source ${EQC_OPERATOR_DIR}/lib/tools.sh
PROCESS_NAME="watchdog_rigol"

config_path=$( cat ${PIDPATH}/config/default.json | jq -r ".config" )
boxid=$( cat ${PIDPATH}/${config_path} | jq -r ".id" | sed -e "s/CoolingBox//g" | awk '{printf "%d", $0;}' )
mpodIP=$( cat ${PIDPATH}/${config_path} | jq -r ".pidControl.mpodIP" )
mpodHVch=$( cat ${PIDPATH}/${config_path} | jq -r ".pidControl.HVch" )

while :
do
    dt=$(python3 ${PIDPATH}/bin/last_pid_timestamp.py | awk '{printf "%.0f", $0;}')
    if [ ${dt} -ge 60 ]; then
	cmd="/nas/dcs/rigol-dp821/RIGOL_single.py ${boxid} turnOff"
	eval ${cmd}
	snmpset -OqvUp +14.12 -v 2c -m +WIENER-CRATE-MIB -c guru ${mpodIP} outputSwitch.${mpodHVch} i 0
	d=$( date "+%Y-%m-%d %H:%M:%S" )
	echo "$( date +%s ) ${d} executed ${cmd}"
	log_record "ERROR: last pid timestamp is more than 60s. Turned off RIGOL LV and HV for safety"
	exit 1
    fi
done

