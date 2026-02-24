#!/bin/sh

ps aux| grep pidServer.py | awk '{print $2;} '| while read pid
do
    kill -9 ${pid}
done
killall -9 pidClient.py

ps aux | grep pidClient.py | grep -v grep | awk '{print $2;}' | while read pid
do
    kill -9 $pid
    echo "killed pidServer process (pid $pid)"
done

ps aux | grep testRead.py | grep -v grep | awk '{print $2;}' | while read pid
do
    kill -9 $pid
    echo "killed testRead.py process (pid $pid)"
done

if [ "$1" == "0" ]; then
    echo "voltage reset to zero"
    ./bin/reset.py --zero
else
    ./bin/reset.py
fi

if [ $? -eq 0 ]; then
    ./test/testPFR.py
fi

