import subprocess

subout = subprocess.run('snmpget -OqvU -v 2c -m +WIENER-CRATE-MIB -c guru 10.0.0.2 outputStatus.u301', capture_output=True, shell=True, check=True )

print( subout.stdout.decode().strip('\n') )
