#!/bin/sh

echo "titech_cool: PATH=${PATH}"

if ! [[ "$PATH" == *"${PWD}/bin"* ]]
then
    export PATH=${PWD}/bin:${PATH}
fi
export PYTHONPATH=${PYTHONPATH}:${PWD}/lib
export PIDPATH=${PWD}
