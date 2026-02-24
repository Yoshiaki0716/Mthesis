#!/usr/bin/env python3

from PidGUI import *
from PidClient import *
from PidMonitor import *
from Logger import *

parser = optparse.OptionParser()
parser.add_option( '-b', action = 'store_true', default = True, dest = 'isBatch' )
parser.add_option( '-v', '--verbose', action = 'store', default = 'info', dest = 'verbosity' )
parser.add_option( '--hostname', action = 'store', default = '127.0.0.1', dest = 'host' )
parser.add_option( '-p', '--port', action = 'store', type = 'int', default = 50007, dest = 'port' )
options, remainder = parser.parse_args( sys.argv )

class PidWidget(QWidget):

    def __init__(self, client):
        super().__init__()
        self.client = client
        self.initUI()

    def initUI(self):
        self.setWindowFlag(Qt.CustomizeWindowHint, True)
        self.setWindowFlag(Qt.WindowCloseButtonHint, False)
        self.setFixedSize(1000, 700)
        self.move(20, 20)
        self.setWindowTitle('CoolingBox Controller')

        self.monControl = VBoxLabel("Control")
        self.button_interlock = ToggleSwitch( 'Cover', ['Free', 'Unlock'], func = self.switchInterlock )
        
        self.monLabel = VBoxLabel("Monitor")
        self.tweaker_setT = SlidingTweaker( 'Set Temperature', '[\u00b0C]', -35, 30, initValue = 25, digits = 1, steps = 65,
                                            func = self.client.do_setT )
        self.tweaker_pelV = SlidingTweaker( 'Peltier Voltage', '[V]', minValue = -12.0, maxValue = 12.0, initValue = 0.0, digits = 2, steps = 240,
                                            func = self.client.do_resetV )

        self.relayLabel = VBoxLabel("Hardware Status")
        
        self.indicators = {
            'Interlock'       : StateIndicator( 'SW Interlock',
                                                states = ['Off', 'Safe', 'Interlocked'],
                                                colors = { 'Off':'gray', 'Safe':'#00aa00', 'Interlocked':'#770000' },
                                                initState = 'Off' ),
            'isStable'        : StateIndicator( '',
                                                states = ['Off', 'Chasing', 'Converging', 'Stable'],
                                                colors = { 'Off':'gray', 'Chasing':'#ffff00', 'Converging':'#aaaaff', 'Stable':'#00aa00' },
                                                initState = 'Off', nameWidth = 0 ),
            'Current'         : ValueIndicator( 'Peltier Current', '[A]', 0, indent=4, digits=2, nameWidth=130 ),
            'T_set_regulated' : ValueIndicator( 'Regulated Set.Temp',        '[\u00b0C]',  initValue = -9999, indent=4, nameWidth=150 ),
            'T_module'        : ValueIndicator( 'Module',          '[\u00b0C]',  initValue = -9999, indent=4 ),
            'T_chuck'         : ValueIndicator( 'Vac.Chuck',       '[\u00b0C]',  initValue = -9999, indent=4 ),
            'T_carrier'       : ValueIndicator( 'Carrier (in)',    '[\u00b0C]',  initValue = -9999, indent=4 ),
            'T_case'          : ValueIndicator( 'Case',            '[\u00b0C]',  initValue = -9999, indent=4 ),
            'T_sink'          : ValueIndicator( 'Heat Sink',       '[\u00b0C]',  initValue = -9999, indent=4 ),
            'T_chiller'       : ValueIndicator( 'Chiller',         '[\u00b0C]',  initValue = -9999, indent=4 ),
            'T_room'          : ValueIndicator( 'Room',            '[\u00b0C]',  initValue = -9999, indent=4 ),
            'DP_carrier'      : ValueIndicator( 'DP Carrier',      '[\u00b0C]',  initValue = -9999, indent=4 ),
            'RH_carrier'      : ValueIndicator( 'RH Carrier',      '[%H]', initValue = -9999, indent=4 ),
            
            'HWInterlock'     : StateIndicator( 'HW Interlock',
                                                states = ['Safe', 'ERR1', 'ERR2', 'ERR3', 'ERR4' ],
                                                colors = { 'Safe':'#00aa00', 'ERR1':'#770000', 'ERR2':'#AA0000', 'ERR3':'#770000', 'ERR4':'#AA0000' },
                                                initState = 'Safe' ),
            
            'Heater'           : StateIndicator( 'Heater', states = ['ON', 'OFF' ], colors = { 'ON':'#ffff00', 'OFF':'gray' }, initState = 'OFF' ),
            'PeltierInterlock' : StateIndicator( 'Peltier Interlock', states = ['Safe', 'Interlocked' ], colors = { 'Safe':'#00aa00', 'Interlocked':'#aa0000' }, initState = 'Safe' ),
            'LVInterlock' : StateIndicator( 'LV Interlock', states = ['Safe', 'Interlocked' ], colors = { 'Safe':'#00aa00', 'Interlocked':'#aa0000' }, initState = 'Safe' ),
            'HVInterlock' : StateIndicator( 'HV Interlock', states = ['Safe', 'Interlocked' ], colors = { 'Safe':'#00aa00', 'Interlocked':'#aa0000' }, initState = 'Safe' ),
            'ChillerInterlock' : StateIndicator( 'Chiller Interlock', states = ['Safe', 'Interlocked' ], colors = { 'Safe':'#00aa00', 'Interlocked':'#aa0000' }, initState = 'Safe' ),
            'PeltierPolPlus'   : StateIndicator( 'Peltier(+)', states = ['Normal', 'Inverted' ], colors = { 'Normal':'gray', 'Inverted':'#ffff00' }, initState = 'Normal' ),
            'PeltierPolMinus'   : StateIndicator( 'Peltier(-)', states = ['Normal', 'Inverted' ], colors = { 'Normal':'gray', 'Inverted':'#ffff00' }, initState = 'Normal' ),
            'Cover1'   : StateIndicator( 'Cover (Front)', states = ['Closed', 'Open' ], colors = { 'Closed':'#00aa00', 'Open':'#ffff00' }, initState = 'Open' ),
            'Cover2'   : StateIndicator( 'Cover (Rear)', states = ['Closed', 'Open' ], colors = { 'Closed':'#00aa00', 'Open':'#ffff00' }, initState = 'Open' )
        }

        self.button_pid = ToggleSwitch( 'PID Control', ['Static Voltage', 'Active'], func = self.switchPID )
        
        self.pidLabel = VBoxLabel("PID Tuning")

        self.tweaker_Kp = SlidingTweaker( 'Kp', '', minValue = 0.0, maxValue = 0.1, initValue = 0.015, digits = 3, steps = 100, func = self.client.do_setKp )
        self.tweaker_Ki = SlidingTweaker( 'Ki', '', minValue = 0.0, maxValue = 2.0, initValue = 1.0, digits = 2, steps = 100, func = self.client.do_setKi )
        self.tweaker_Kd = SlidingTweaker( 'Kd', '', minValue = 0.0, maxValue = 5.0, initValue = 3.5, digits = 3, steps = 100, func = self.client.do_setKd )
        self.tweaker_maxP = SlidingTweaker( 'max (P)', '', 0, 2.0, 0.4, 3, steps = 100, func = self.client.do_setMaxP )
        self.tweaker_maxI = SlidingTweaker( 'max (I)', '', 0, 0.3, 0.03, 3, steps = 300, func = self.client.do_setMaxI )
        self.tweaker_maxD = SlidingTweaker( 'max (D)', '', 0, 2.0, 0.6, 3, steps = 100, func = self.client.do_setMaxD )
        self.tweaker_Idur = SlidingTweaker( 'Integral Duration', '[ticks]', minValue = 0, maxValue = 30, initValue = 7, digits = 0, steps = 30,
                                            func = self.client.do_setIduration )
        self.tweaker_tOffsetP = SlidingTweaker( 'Time Offset (P)', '[s]', minValue = 0, maxValue = 600, initValue = 60, digits = 0, steps = 60,
                                               func = self.client.do_setTimeOffsetP )
        self.tweaker_tOffsetI = SlidingTweaker( 'Time Offset (I)', '[s]', minValue = 0, maxValue = 600, initValue = 60, digits = 0, steps = 60,
                                               func = self.client.do_setTimeOffsetI )
        self.tweaker_tOffsetD = SlidingTweaker( 'Time Offset (D)', '[s]', minValue = 0, maxValue = 600, initValue = 60, digits = 0, steps = 60,
                                               func = self.client.do_setTimeOffsetD )

        # Quit button
        self.hbox_top1 = QHBoxLayout()
        self.hbox_top2 = QHBoxLayout()
        self.hbox_bottom = QHBoxLayout()

        global options
        self.host = StringIndicator('Server Host', '{}:{}'.format( options.host, options.port ) )
        self.hbox_top1.addLayout( self.host.hbox )
        
        self.hbox_top2.addLayout( self.indicators['Interlock'].hbox )
        self.hbox_top2.addLayout( self.button_interlock.hbox )
        self.hbox_top2.addLayout( self.indicators['HWInterlock'].hbox )
        self.hbox_top2.addStretch(1)
        
        self.button = QPushButton( 'Quit', self )
        self.button.clicked.connect( self.clickQuit )

        self.hbox_bottom.addStretch(1)
        self.hbox_bottom.addWidget( self.button )

        # Tabs
        self.tabs = QTabWidget()
        childTabs = { 'Basic': QWidget(), 'Advanced' : QWidget() }
        self.tabs.resize(400, 450)

        for name, tab in childTabs.items():
            self.tabs.addTab( tab, name )

        # Global layout
        vbox1 = QVBoxLayout()
        vbox2 = QVBoxLayout()
        vbox3 = QVBoxLayout()
        vbox4 = QVBoxLayout()

        self.pid_hbox = QHBoxLayout()
        self.pid_hbox.addLayout( self.button_pid.hbox )
        self.pid_hbox.addLayout( self.indicators['isStable'].hbox )
        
        vbox1.addLayout( self.monControl.hbox )
        vbox1.addLayout( self.pid_hbox )
        vbox1.addLayout( self.tweaker_setT.hbox )
        vbox1.addLayout( self.indicators['T_set_regulated'].hbox )
        vbox1.addLayout( self.tweaker_pelV.hbox )
        vbox1.addLayout( self.indicators['Current'].hbox )
        vbox1.addStretch(0)
        
        vbox2.addLayout( self.monLabel.hbox )
        vbox2.addLayout( self.indicators['T_module'].hbox )
        vbox2.addLayout( self.indicators['T_chuck'].hbox )
        vbox2.addLayout( self.indicators['T_carrier'].hbox )
        vbox2.addLayout( self.indicators['T_case'].hbox )
        vbox2.addLayout( self.indicators['T_sink'].hbox )
        vbox2.addLayout( self.indicators['T_chiller'].hbox )
        vbox2.addLayout( self.indicators['T_room'].hbox )
        vbox2.addLayout( self.indicators['DP_carrier'].hbox )
        vbox2.addLayout( self.indicators['RH_carrier'].hbox )
        vbox2.addStretch(0)
        
        vbox4.addLayout( self.relayLabel.hbox )
        vbox4.addLayout( self.indicators['PeltierPolPlus'].hbox )
        vbox4.addLayout( self.indicators['PeltierPolMinus'].hbox )
        vbox4.addLayout( self.indicators['Heater'].hbox )
        vbox4.addLayout( self.indicators['Cover1'].hbox )
        vbox4.addLayout( self.indicators['Cover2'].hbox )
        vbox4.addLayout( self.indicators['ChillerInterlock'].hbox )
        vbox4.addLayout( self.indicators['PeltierInterlock'].hbox )
        vbox4.addLayout( self.indicators['LVInterlock'].hbox )
        vbox4.addLayout( self.indicators['HVInterlock'].hbox )
        vbox4.addStretch(0)
        
        vbox3.addLayout( self.pidLabel.hbox )
        vbox3.addLayout( self.tweaker_Kp.hbox )
        vbox3.addLayout( self.tweaker_Ki.hbox )
        vbox3.addLayout( self.tweaker_Kd.hbox )
        vbox3.addLayout( self.tweaker_maxP.hbox )
        vbox3.addLayout( self.tweaker_maxI.hbox )
        vbox3.addLayout( self.tweaker_maxD.hbox )
        vbox3.addLayout( self.tweaker_Idur.hbox )
        vbox3.addLayout( self.tweaker_tOffsetP.hbox )
        vbox3.addLayout( self.tweaker_tOffsetI.hbox )
        vbox3.addLayout( self.tweaker_tOffsetD.hbox )
        vbox3.addStretch(3)

        hboxBasic = QHBoxLayout()
        hboxBasic.addLayout( vbox1 )
        hboxBasic.addLayout( vbox2 )
        hboxBasic.addLayout( vbox4 )

        childTabs['Basic'].setLayout( hboxBasic )
        childTabs['Advanced'].setLayout( vbox3 )

        self.layout = QVBoxLayout( self )
        self.layout.addLayout( self.hbox_top1 )
        self.layout.addLayout( self.hbox_top2 )
        self.layout.addWidget(self.tabs)
        self.layout.addLayout( self.hbox_bottom )
        self.setLayout( self.layout )

        self.show()

    def switchInterlock( self, flag ):
        print( 'switchInterlock: flag = {}'.format( flag ) )
        
        if flag:
            pass
        else:
            d = AskingDialog("Interlock is going to be disabled.\nAre you sure opening the Cooling Box?").exec_()
            if d == QDialog.Accepted :
                msg = self.client.do_unlock( '' )
                if 'accepted' in msg:
                    self.indicators['Interlock'].changeValue( 'Off' )
                    self.button_interlock.button.setDown( False )
                    self.button_interlock.button.setText( self.button_interlock.stateNames[0] )
                    self.button_interlock.state = False
                    self.tweaker_setT.disable( False )
                    self.tweaker_setT.updateValue( 25.0 )
                    OKDialog(msg).exec_()
                else:
                    self.button_interlock.button.setDown( True )
                    self.button_interlock.button.setText( self.button_interlock.stateNames[1] )
                    self.button_interlock.state = True
                    OKDialog(msg).exec_()
            else:
                self.button_interlock.button.setDown( True )
                self.button_interlock.button.setText( self.button_interlock.stateNames[1] )
                self.button_interlock.state = True
                
    def switchPID( self, flag ):
        
        if flag:
            d = AskingDialog("Are you sure starting the PID feedback?").exec_()
            if d == QDialog.Accepted:
                msg = self.client.do_restartPID( '' )
                if 'accepted' in msg:
                    self.button_pid.button.setDown( True )
                    self.button_pid.button.setText( self.button_pid.stateNames[1] )
                    self.button_pid.state = True
                    self.tweaker_pelV.disable( True )
                    self.tweaker_setT.disable( False )
                    OKDialog(msg).exec_()
                else:
                    self.indicators['isStable'].changeValue( 'Off' )
                    self.tweaker_pelV.disable( False )
                    self.tweaker_setT.disable( True )
                    self.button_pid.button.setDown( False )
                    self.button_pid.button.setText( self.button_pid.stateNames[0] )
                    self.button_pid.state = False
                    self.tweaker_pelV.disable( False )
                    self.tweaker_setT.disable( True )
                    OKDialog(msg).exec_()
            else:
                self.indicators['isStable'].changeValue( 'Off' )
                self.tweaker_pelV.disable( False )
                self.tweaker_setT.disable( True )
                self.button_pid.button.setDown( False )
                self.button_pid.button.setText( self.button_pid.stateNames[0] )
                self.button_pid.state = False
                self.tweaker_pelV.disable( False )
                self.tweaker_setT.disable( True )
        else:
            d = AskingDialog("Are you sure stopping the PID feedback?").exec_()
            if d == QDialog.Accepted:
                msg = self.client.do_stopPID( '' )
                if 'accepted' in msg:
                    self.indicators['isStable'].changeValue( 'Off' )
                    self.tweaker_pelV.disable( False )
                    self.tweaker_setT.disable( True )
                    self.button_pid.button.setDown( False )
                    self.button_pid.button.setText( self.button_pid.stateNames[0] )
                    self.button_pid.state = False
                    self.tweaker_pelV.disable( False )
                    self.tweaker_setT.disable( True )
                    OKDialog(msg).exec_()
                else:
                    self.button_pid.button.setDown( True )
                    self.button_pid.button.setText( self.button_pid.stateNames[1] )
                    self.button_pid.state = True
                    self.tweaker_pelV.disable( True )
                    self.tweaker_setT.disable( False )
                    OKDialog(msg).exec_()
            else:
                self.button_pid.button.setDown( True )
                self.button_pid.button.setText( self.button_pid.stateNames[1] )
                self.button_pid.state = True
                self.tweaker_pelV.disable( True )
                self.tweaker_setT.disable( False )


    def clickQuit( self ):
        self.client.do_quit( '' )
        app.quit()

    def updateParams( self ):
        global params
        global mon

        while True:
            if self.client.kill: break

                            
            try:
                # print( 'fetching log data...')
                command = 'sync all' if len( mon.data ) == 0 else ('sync ' + str( mon.data['elapsed'][-1] ) )
                self.client.logger.debug( 'sending command ' + command )
                
                sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
                sock.connect( (options.host, options.port) )
                sock.sendall( command.encode() )
        
                
                self.client.logger.debug( 'sync command sent...' )
                
                buf = b''
                while True:
                    data = sock.recv(1024)
                    buf += data
                    
                    self.client.logger.verbose( 'received data size: ' + str( len(data) ) )
                    
                    if repr(data).find( '\\n' ) < 0:
                        continue
                    else:
                        break
                    
                self.client.logger.debug( 'completed data reception. size = : ' + str( len(buf) ) )
                
                sock.close()
                
                j = json.loads( buf.decode("utf-8") )
                
                if len(j) == 0:
                    pass
                else:
                    mon.t0 = j[0]
                    for k, v in j[1].items():
                        if k in mon.data:
                            mon.data[k] += v
                        else:
                            mon.data.update( { k:v } )
                            mon.status = j[2]

                for name,log in mon.data.items():
                    
                    time.sleep(0.05)

                    if len( log ) == 0:
                        continue
                        
                    self.client.logger.debug( "{} : {}".format( name, log[-1] ) )
                    
                    if name in self.indicators:
                        if name == 'isStable':
                            pass
                        else:
                            self.indicators[name].changeValue( log[-1] )
                            
                    else:
                        if name == 'T_carrier':
                            self.indicators['T_carrier'].changeValue( log[-1] )
                        elif name == 'RH_carrier':
                            self.indicators['RH_carrier'].changeValue( log[-1] )
                        elif name == 'DP_carrier':
                            self.indicators['DP_carrier'].changeValue( log[-1] )
                        elif name == 'T_chuck_ref':
                            self.indicators['T_room'].changeValue( log[-1] )
                        elif name == 'I':
                            self.indicators['Current'].changeValue( log[-1] )
                        elif name == 'V':
                            self.tweaker_pelV.updateValue( log[-1] )
                        elif name == 'interlockStatus': 
                            self.indicators['HWInterlock'].changeValue( 'ERR'+str(log[-1]) if log[-1]>0 else 'Safe' )
                        elif name == 'RS_heater': 
                            self.indicators['Heater'].changeValue( 'ON' if log[-1]==1 else 'OFF' )
                        elif name == 'RS_peltier': 
                            self.indicators['PeltierInterlock'].changeValue( 'Safe' if log[-1]>0 else 'Interlocked' )
                        elif name == 'RS_chiller': 
                            self.indicators['ChillerInterlock'].changeValue( 'Safe' if log[-1]>0 else 'Interlocked' )
                        elif name == 'RS_lowVoltage': 
                            self.indicators['LVInterlock'].changeValue( 'Safe' if log[-1]>0 else 'Interlocked' )
                        elif name == 'RS_highVoltage': 
                            self.indicators['HVInterlock'].changeValue( 'Safe' if log[-1]>0 else 'Interlocked' )
                        elif name == 'RS_pelPlus': 
                            self.indicators['PeltierPolPlus'].changeValue( 'Inverted' if log[-1]>0 else 'Normal' )
                        elif name == 'RS_pelMinus': 
                            self.indicators['PeltierPolMinus'].changeValue( 'Inverted' if log[-1]>0 else 'Normal' )
                        elif name == 'RS_lockState1': 
                            self.indicators['Cover1'].changeValue( 'Open' if log[-1]==1 else 'Closed' )
                        elif name == 'RS_lockState2': 
                            self.indicators['Cover2'].changeValue( 'Open' if log[-1]==1 else 'Closed' )
                        else:
                            pass #print( '  logs unmatched: ' + name )

            except Exception as e:
                
                if self.client.kill: return
                
                self.client.logger.error( 'getMon() Exception: ' + str(e) )
                print( str(e) )
                print(traceback.format_exc()) 
                

                        
            msg = self.client.do_params( 'quiet' )

            for name,value in self.client.params.items():
                
                self.client.logger.debug( '{} {}'.format( name, value ) )
                
                if name in self.indicators:
                    if name == 'Interlock':
                        if self.button_interlock.state:
                            self.indicators[name].changeValue( value )
                    elif name == 'isStable':
                        if self.button_pid.state:
                            self.indicators[name].changeValue( value )
                    else:
                        self.indicators[name].changeValue( value )
                else:
                    if name == 'T_set':
                        if self.tweaker_setT.isChanging:
                            pass
                        else:
                            self.tweaker_setT.updateValue( value )
                    elif name == 'T_set_regulated':
                        self.indicators['T_set_regulated'].changeValue( value )
                    elif name == 'Kp':
                        self.tweaker_Kp.updateValue( value )
                    elif name == 'Ki':
                        self.tweaker_Ki.updateValue( value )
                    elif name == 'Kd':
                        self.tweaker_Kd.updateValue( value )
                    elif name == 'maxPfeedback':
                        self.tweaker_maxP.updateValue( value )
                    elif name == 'maxIfeedback':
                        self.tweaker_maxI.updateValue( value )
                    elif name == 'maxDfeedback':
                        self.tweaker_maxD.updateValue( value )
                    elif name == 'integralDuration':
                        self.tweaker_Idur.updateValue( value )
                    elif name == 'tOffsetP':
                        self.tweaker_tOffsetP.updateValue( value )
                    elif name == 'tOffsetI':
                        self.tweaker_tOffsetI.updateValue( value )
                    elif name == 'tOffsetD':
                        self.tweaker_tOffsetD.updateValue( value )
                    elif name == 'isInterlockEnabled':
                        interlockStatus = self.client.params['isInterlocked']
                        if value == False:
                            self.button_interlock.update( value )
                            self.button_interlock.disable()
                            self.indicators['Interlock'].changeValue( 'Off' )
                        elif value == True and interlockStatus == False :
                            self.button_interlock.update( value )
                            self.button_interlock.button.setDown( True )
                            self.button_interlock.enable()
                            self.indicators['Interlock'].changeValue( 'Safe' )
                        else:
                            self.button_interlock.update( value )
                            self.button_interlock.button.setDown( True )
                            self.indicators['Interlock'].changeValue( 'Interlocked' )
                    elif name == 'isPidActivated':
                        try:
                            isConverging = mon.data['isConverging'][-1]
                            isStable     = mon.data['isStable'][-1]
                            
                            if isStable:
                                self.indicators['isStable'].changeValue( 'Stable' )
                            elif isConverging:
                                self.indicators['isStable'].changeValue( 'Converging' )
                            else:
                                self.indicators['isStable'].changeValue( 'Chasing' )
                                
                            if value == 0:
                                self.button_pid.update( False )
                                self.button_pid.button.setDown( False )
                                self.indicators['isStable'].changeValue( 'Off' )
                                self.tweaker_pelV.disable( False )
                                self.tweaker_setT.disable( True )
                                
                            elif value == 1:
                                self.button_pid.update( True )
                                self.button_pid.button.setDown( True )
                                self.tweaker_pelV.disable( True )
                                self.tweaker_setT.disable( False )
                                
                            else:
                                pass #print( '  params unmatched : ' + name )
                        except:
                            pass

            time.sleep(2)



if __name__ == '__main__':

    logger = Logger( '{}/logs/client.log'.format( os.getenv('PIDPATH') ), verbosity = options.verbosity )
    mon = PidMonitor(500)
    
    client = PidClient( logger, mon, options )
    
    app = QApplication(sys.argv)
    wid = PidWidget( client )
    
    tUpdate = threading.Thread( target = wid.updateParams )
    tUpdate.start()

    app.exec_()

    tUpdate.join()
    
    sys.exit()

