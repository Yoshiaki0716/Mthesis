import time
import optparse
#from importlib import reload

#import ArduinoSingleBoxDcsController
#reload(ArduinoSingleBoxDcsController)
#print(ArduinoSingleBoxDcsController.__file__)

from ArduinoSingleBoxDcsController import *

parser = optparse.OptionParser()
parser.add_option("-m","--mode", dest="mode",type="int",default=None)
(options, args) = parser.parse_args()


try:
    ctrl = ArduinoSingleBoxDcsController({"port" : "/dev/ttyUSB_ArduinoSingleBoxController", "baud":9600})

#    print(dir(ctrl))

    if options.mode is not None:
        if 0 <= options.mode <= 7:
            mode = options.mode
            ctrl.reset()
            time.sleep(1)
            ctrl.interlockTest(mode)
        else:
            print("mode must be [0-7]")

    else:
        print("please use the option [-m]")
        
except Exception as e:
    print('Exception raised in checking ArduinoDcsController:' ,e)
    
