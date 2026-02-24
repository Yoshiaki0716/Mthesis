import traceback

def invokeDevices( config, argName=None ):
    devices = {}
    
    for name, cfg in config.items():
        
        if argName != None:
            if name != argName:
                continue
        
        try:
            classname = cfg['class']
            modulename = classname
            if 'module' in cfg:
                modulename = cfg['module']

            module = __import__( modulename )
            c = getattr( module, classname )
            instance = c( cfg )
            devices.update( { name : instance } )
        except:
            print( traceback.format_exc() )
            exit(1)

    return devices



if __name__ == '__main__':
    
    config = { 	    "peltierControl" : { "module"   : "PeltierControl",
				 "class"    : "TakasagoKXControl",
			         "port"     : "/dev/ttyUSB1", 
			         "address"  : 1,
			         "baud"     : 9600
			       },
                    
                    "DcsController" : { "module" : "ArduinoSingleBoxDcsController",
				        "class"  : "ArduinoSingleBoxDcsController",
				        "port"   : "/dev/ttyUSB0",
				        "baud"   : 9600 }
                    }
    
    devs = invokeDevices( config )

    print( devs['peltierControl'].getV() )
    devs['DcsController'].readAll()
    devs['DcsController'].print()

