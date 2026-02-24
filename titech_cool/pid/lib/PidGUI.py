import sys
import time
import functools as ft
from PyQt5.QtWidgets import *
from PyQt5.QtGui import *
from PyQt5.QtCore import *


class ToggleSwitch(QWidget):
    
    def __init__(self, name, stateNames = ['Inactive', 'Activated'], indent=4, func = None ):
        super().__init__()

        self.name = name
        self.isEnabled = True
        self.state = False
        self.stateNames = stateNames
        self.func = func
        self.lastCommandTime = float( time.time() ) - 60
        
        self.hbox = QHBoxLayout()
        
        indent = ''.join( [ ' ' for k in range(0, indent) ] )
        self.label = QLabel( indent + '{}: '.format(name) )

        self.button = QPushButton( self.stateNames[0], self)
        self.button.setFixedWidth(140)
        self.button.setCheckable( True )
        self.button.toggled.connect( self.onClick )

        self.hbox.addWidget( self.label )
        self.hbox.addWidget( self.button )
        self.hbox.addStretch(1)

    def enable(self):
        self.isEnabled = True
        self.button.setEnabled( True )

    def disable(self):
        self.isEnabled = False
        self.button.setEnabled( False )

    def onClick( self, checked ):
        if not self.isEnabled:
            print( 'button is disabled' )
            self.button.setDown( self.state )
            if self.state:
                self.button.setText( self.stateNames[1] )
            else:
                self.button.setText( self.stateNames[0] )
            return

        self.state = not self.state
        print( self.name + ': going to change the state to ' + self.stateNames[self.state] )
        self.lastCommandTime = float( time.time() ) - 60
        
        if self.state:
            self.button.setDown( True )
            self.button.setText( self.stateNames[1] )
        else:
            self.button.setDown( False )
            self.button.setText( self.stateNames[0] )

        if self.func != None:
            self.func( self.state )
    
    def update( self, flag ):
        self.state = flag
        
        if self.state:
            self.button.setDown( True )
            self.button.setText( self.stateNames[1] )
        else:
            self.button.setDown( False )
            self.button.setText( self.stateNames[0] )


class SlidingTweaker(QWidget):
    def __init__(self, name, unit, minValue, maxValue,
                 initValue = 0.0,
                 digits    = 2,
                 steps     = 100,
                 indent    = 4,
                 nameWidth = 130,
                 confirm   = True,
                 func      = None ):

        self.name = name
        self.unit = unit
        self.isEnabled = True
        self.value = initValue
        self.isChanging = False
        self.lastCommandTime = float( time.time() )

        self.digits = digits
        self.digitsForm = '{:.' + str(digits) + 'f}'
        self.steps = steps
        self.func = func
        
        super().__init__()
        
        self.hbox = QHBoxLayout()

        self.minValue = minValue
        self.maxValue = maxValue
        self.value    = initValue

        indent = ''.join( [ ' ' for k in range(0, indent) ] )
        self.label = QLabel( indent + '{}: '.format(name) )
        self.label.setFixedWidth( nameWidth )
        
        self.sld = QSlider(Qt.Horizontal, self)
        self.sld.move( minValue, maxValue )
        self.sld.setMinimum(0)
        self.sld.setMaximum(steps)
        self.sld.setFocusPolicy(Qt.NoFocus)
        self.sld.setGeometry(0, 0, 400, 30)
        self.sld.setValue( int( (initValue - minValue)*self.steps/(maxValue-minValue) ) )

        self.ind = QLineEdit( self )
        self.ind.resize(20,10)
        self.ind.setFixedWidth(60)
        self.ind.setAlignment(Qt.AlignRight)
        self.ind.setReadOnly( False )
        self.ind.setText( self.digitsForm.format( initValue ) )
        
        self.unitLabel = QLabel( unit )
        
        self.hbox.addWidget( self.label )
        self.hbox.addWidget( self.sld )
        self.hbox.addWidget( self.ind )
        self.hbox.addWidget( self.unitLabel )
        #self.hbox.addStretch(1)

        # change the indicator's value while sliding the gauge
        self.sld.valueChanged.connect(self.changeIndicator)

        # change the indicator's value while sliding the gauge
        self.ind.returnPressed.connect(self.changeViaLineEdit)

        # ask actuallyg reflect the change when releasing the slide
        self.sld.sliderReleased.connect(self.changeViaSlider)
        
    def changeIndicator(self):

        self.isChanging = True
        
        if self.isEnabled:
            value_tmp = self.minValue + self.sld.value()/float(self.steps)*(self.maxValue-self.minValue)
            self.ind.setText( self.digitsForm.format(value_tmp) )
        else:
            self.sld.setValue( int( (self.value - self.minValue)*self.steps/(self.maxValue-self.minValue) ) )

    def changeViaSlider(self):
        value_prev = self.value
        value = round( self.minValue + self.sld.value()/float(self.steps)*(self.maxValue-self.minValue), self.digits )

        if self.isEnabled:
            
            d = AskingDialog("Are you sure changing the parameter {} to {} {}?".format( self.name, self.digitsForm.format(value), self.unit ) ).exec_()
            if d == QDialog.Accepted:
                self.value = value
                self.ind.setText( self.digitsForm.format(self.value) )
                self.func( self.digitsForm.format( self.value ) )
                OKDialog("The parameter {} changed to {} {}".format( self.name, self.digitsForm.format(self.value), self.unit ) )
                self.lastCommandTime = float( time.time() )
            else:
                self.ind.setText( self.digitsForm.format(self.value) )
                self.sld.setValue( int( (self.value - self.minValue)*self.steps/(self.maxValue-self.minValue) ) )
        else:
            self.value = value
            self.sld.setValue( int( (self.value - self.minValue)*self.steps/(self.maxValue-self.minValue) ) )
            self.ind.setText( self.digitsForm.format(self.value) )

        self.isChanging = False

    def changeViaLineEdit(self):
        print('changeViaLineEdit() called' )
        value_prev = self.value

        try:
            value = float( self.ind.text() )

            if self.isEnabled:
                
                d = AskingDialog("Are you sure changing the parameter {} to {} {}?".format( self.name, self.digitsForm.format(value), self.unit ) ).exec_()
                if d == QDialog.Accepted:
                    self.value = value
                    self.ind.setText( self.digitsForm.format(self.value) )
                    self.func( self.digitsForm.format( self.value ) )
                    OKDialog("The parameter {} changed to {} {}".format( self.name, self.digitsForm.format(self.value), self.unit ) )
                    self.lastCommandTime = float( time.time() )
                else:
                    self.ind.setText( self.digitsForm.format(self.value) )
                    self.sld.setValue( int( (self.value - self.minValue)*self.steps/(self.maxValue-self.minValue) ) )
            else:
                self.value = value
                self.sld.setValue( int( (self.value - self.minValue)*self.steps/(self.maxValue-self.minValue) ) )
                self.ind.setText( self.digitsForm.format(self.value) )

        except:
            self.ind.setText( self.digitsForm.format(value_prev) )

    def updateValue(self, value):
        self.sld.setValue( int( (self.value - self.minValue)*self.steps/(self.maxValue-self.minValue) ) )
        self.ind.setText( self.digitsForm.format(self.value) )
        self.value = value
        
    def enable( self, flag ):
        self.isEnabled = flag
        self.ind.setReadOnly( not flag )

    def disable( self, flag ):
        self.isEnabled = not flag
        self.ind.setReadOnly( flag )



class ValueIndicator(QWidget):
    def __init__(self, name, unit, initValue = 0.0, digits=2, indent=4, nameWidth=130 ):

        self.digitsForm = '{:.' + str(digits) + 'f}'
        
        super().__init__()
        
        self.hbox = QHBoxLayout()

        self.value    = initValue
        
        indent = ''.join( [ ' ' for k in range(0, indent) ] )
        self.label = QLabel( indent + '{}: '.format(name) )
        self.label.setFixedWidth( nameWidth )
        
        self.ind = QLineEdit( self )
        self.ind.resize(20,10)
        self.ind.setFixedWidth(80)
        self.ind.setAlignment(Qt.AlignRight)
        self.ind.setReadOnly( True )
        self.ind.setText( self.digitsForm.format( initValue ) )
        
        self.unitLabel = QLabel( unit )
        
        self.hbox.addWidget( self.label )
        self.hbox.addWidget( self.ind )
        self.hbox.addWidget( self.unitLabel )
        self.hbox.addStretch(2)

    def changeValue(self, value):
        self.value = value
        self.ind.setText( self.digitsForm.format(self.value) )


class StringIndicator(QWidget):
    def __init__(self, name, initValue = '', indent=4, nameWidth=100 ):

        super().__init__()
        
        self.hbox = QHBoxLayout()

        self.value    = initValue
        
        indent = ''.join( [ ' ' for k in range(0, indent) ] )
        self.label = QLabel( indent + '{}: '.format(name) )
        self.label.setFixedWidth( nameWidth )
        
        self.ind = QLineEdit( self )
        self.ind.setFixedWidth(300)
        self.ind.setAlignment(Qt.AlignLeft)
        self.ind.setReadOnly( True )
        self.ind.setText( self.value )
        
        self.hbox.addWidget( self.label )
        self.hbox.addWidget( self.ind )
        self.hbox.addStretch(2)

    def changeValue(self, value):
        self.value = value
        self.ind.setText( value )


class StateIndicator(QWidget):
    def __init__(self, name,
                 states = ['OFF', 'ON'],
                 colors = { 'OFF':'gray', 'ON': 'green' },
                 initState = 'OFF', indent=4, nameWidth=130 ):

        super().__init__()
        
        self.hbox = QHBoxLayout()
        self.states = states
        self.colors = colors
        self.state = initState

        indent = ''.join( [ ' ' for k in range(0, indent) ] )
        self.label = QLabel( indent + '{}: '.format(name) )
        self.label.setFixedWidth( nameWidth )
        
        self.ind = QLineEdit( self )
        self.ind.resize(20,10)
        self.ind.setFixedWidth(100)
        self.ind.setAlignment(Qt.AlignCenter)
        self.ind.setReadOnly( True )
        self.ind.setText( self.state )
        
        self.hbox.addWidget( self.label )
        self.hbox.addWidget( self.ind )
        self.hbox.addStretch(2)

        self.changeValue( self.state )

    def changeValue(self, state):
        self.state = state
        self.ind.setText( self.state )
        self.ind.setStyleSheet( 'QLineEdit { background: ' + self.colors[state] + '}' )


class VBoxLabel(QWidget):

    def __init__(self, title):
        super().__init__()
        
        self.label = QLabel( title )
        font = self.label.font()
        font.setPointSize(20)
        font.setBold( True )
        self.label.setFont( font )
        self.hbox = QHBoxLayout()
        self.hbox.addWidget( self.label )
        self.hbox.addStretch(1)



class AskingDialog( QDialog ):
    def __init__(self, msg):
        super().__init__()
        
        d_vbox = QVBoxLayout()
        d_hbox = QHBoxLayout()
        msgLabel = QLabel( msg )
        
        proceed = QPushButton("Proceed", self)
        proceed.clicked.connect( self.accept )
        
        cancel  = QPushButton("Cancel", self)
        cancel.clicked.connect( self.reject )
        
        d_hbox.addWidget( proceed )
        d_hbox.addWidget( cancel )
        
        d_vbox.addWidget( msgLabel )
        d_vbox.addLayout( d_hbox )
        
        self.setLayout( d_vbox )
        self.move(100, 200)
        self.setWindowTitle("Confirmation")
        self.setWindowModality(Qt.ApplicationModal)

class OKDialog( QDialog ):
    def __init__(self, msg):
        super().__init__()
        
        d_vbox = QVBoxLayout()
        d_hbox = QHBoxLayout()
        msgLabel = QLabel( msg )
        
        ok = QPushButton("OK", self)
        ok.clicked.connect( self.close )
        
        d_hbox.addWidget( ok )
        
        d_vbox.addWidget( msgLabel )
        d_vbox.addLayout( d_hbox )
        
        self.setLayout( d_vbox )
        self.move(100, 200)
        self.setWindowTitle("Confirmation")
        self.setWindowModality(Qt.ApplicationModal)
