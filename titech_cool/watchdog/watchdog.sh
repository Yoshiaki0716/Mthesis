#!/bin/bash

cd /home/atlasj/titech_cool/pid
source setup.sh
cd ../watchdog

lastTS=`influx -host atlastit01.kek.jp -database 'dcsDB' -execute 'select * from TitechCoolingBox02 order by time desc limit 1' -format 'csv' | tail -n1 | cut -d',' -f2`
now=`date "+%s%9N"`

isRun=`ps aux| grep -i pidserver | grep python3 | wc -l`
if [ "$isRun" == "0" ]; then
    exit 0
fi

lastTS_form=`date -d @${lastTS::-9} +%Y-%m-%d\ %H:%M:%S`
now_form=`date -d @${now::-9} +%Y-%m-%d\ %H:%M:%S`
# show in sec
diff=`expr $now - $lastTS`

logfile=/home/atlasj/titech_cool/watchdog/watchdog.log

cat ${logfile} | tail -n1000 > tmp; mv tmp ${logfile}

influx -host atlastit01.kek.jp -database 'dcsDB' -execute 'select * from TitechCoolingBox02 order by time desc limit 1' -format 'json' > /tmp/watch.json
opt=`python3 /home/atlasj/titech_cool/watchdog/getLastOpt.py`

echo "now is $now_form ${now::-6}, last is $lastTS_form  ${lastTS::-6}  | diff = ${diff::-6}"  >> ${logfile}

embargo=/home/atlasj/titech_cool/watchdog/error

if [ "${diff::-6}" -gt 120000 ]; then
    if [ ! -e ${embargo} ]; then
	
	ps aux | grep -i pidClient | grep python3 | awk '{print $2;}' | while read c; do kill -9 $c; done
	ps aux | grep -i pidWidget | grep python3 | awk '{print $2;}' | while read c; do kill -9 $c; done
	echo "killed all PID clients" >> ${logfile}
	
	email_to="To: hideude@gmail.com, asuzuki@hep.phys.titech.ac.jp"
	pid=`ps aux|grep "python3" | grep "pidServer.py" | awk '{print $2;}'`
	
	pidClient.py <<EOF >> ${logfile} &
relaunch_server ${opt}
y
quit
EOF
	
	sleep 90
	
	pid_new=`ps aux|grep "python3" | grep "pidServer.py" | awk '{print $2;}'`
	
	if [ "${pid}" == "${pid_new}" ]; then
	    kill -9 ${pid}
	    echo "${now_form} force killing pidServer.py" >> ${logfile}
	    /home/atlasj/titech_cool/pid
	    ./reset.sh
	    sleep 10
	    pidServer.py ${opt}
	    
	    msg="No influxdb records since ${lastTS_form}"
	    
	    {
		echo "From: itkqc@atlaspc28.kek.jp"
		echo ${email_to}
		echo "Subject: [${now_form} -- ERROR] ITkPixV1 CoolingBox Watchdog -- FORCE RESTART"
		echo
		echo ${msg}
		echo 
		echo "Relaunch procedure failed; force killing pidServer.py"
		echo 
		echo "Check Grafana: http://atlastit01.kek.jp:3000/d/9T_DEUJ7z/coolingbox-itkpix-01?orgId=1&refresh=10s"
	    } | /usr/sbin/sendmail -i -t
	    touch ${embargo}
	    
	    
	else
	    
	    msg="No influxdb records since ${lastTS_form}"
	    
	    {
		echo "From: itkqc@atlaspc28.kek.jp"
		echo ${email_to}
		echo "Subject: [${now_form} -- ERROR] ITkPixV1 CoolingBox Watchdog -- restarting the server"
		echo
		echo ${msg}
		echo "Check Grafana: http://atlastit01.kek.jp:3000/d/9T_DEUJ7z/coolingbox-itkpix-01?orgId=1&refresh=10s"
	    } | /usr/sbin/sendmail -i -t
	    touch ${embargo}
	    
	fi
	
    else
	
	echo "${now_form} embargo is present, skipping" >> ${logfile}
	
    fi
    
    rm -f ${embargo}
    
fi



