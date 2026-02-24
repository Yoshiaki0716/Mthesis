#!/usr/bin/env python3

from PidGUI import *
from PidClient import *
from PidMonitor import *
from Logger import *
import subprocess

class SerialNumberSwitcher(QWidget):

    def __init__(self):
        super().__init__()
        self.initUI()

    def initUI(self):
        
        self.setWindowTitle('Module SerialNumber Switcher')
        self.setFixedSize(800, 200)
        self.move(50, 50)
        self.label = QLabel(self)
        self.label.setText("Input Module Serial Number to Load")
        self.inputField = QLineEdit( self )
        self.inputField.resize(200, 32)
        
        self.button = QPushButton('Switch', self)
        self.button.clicked.connect( self.register )
        
        self.layout = QVBoxLayout( self )
        self.layout.addStretch(0)
        self.layout.addWidget( self.label )
        self.layout.addWidget( self.inputField )
        self.layout.addWidget( self.button )
        self.layout.addStretch(0)
        
        self.show()
        
    def register(self):
        serial_number = self.inputField.text().upper()
        
        if len( serial_number ) != 14:
            msgBox = QMessageBox.warning(
                self,
                "Warning",
                "[SN_LENGTH_ERROR] Invalid Serial Number!",
                QMessageBox.Ok,
            )
            return
            
        if serial_number.find('20UPGM2') == -1:
            msgBox = QMessageBox.warning(
                self,
                "Warning",
                "[SN_HEADER_ERROR] Invalid Serial Number!",
                QMessageBox.Ok,
            )
            return
            
        try:
            int( serial_number[7:] )
        except:
            msgBox = QMessageBox.warning(
                self,
                "Warning",
                "[SN_NUMBER_ERROR] Invalid Serial Number!",
                QMessageBox.Ok,
            )
            return
            
        msgBox = QMessageBox.information(
            self,
            "Info",
            f"Accepted Switching Serial Number: {serial_number}",
            QMessageBox.Ok,
        )
        
        EQC_OPERATOR_DIR = os.environ.get('EQC_OPERATOR_DIR', '')
        HOSTNAME = os.environ.get('HOSTNAME', '')

        cmd=f"{EQC_OPERATOR_DIR}/setup.sh -c {EQC_OPERATOR_DIR}/configs/setup/{HOSTNAME}_setup.json -m {serial_number}"
        
        out = subprocess.run( cmd, shell=True, capture_output = True, text = True )
        print( out.stdout )
        
        app.quit()



if __name__ == '__main__':

    app = QApplication(sys.argv)
    wid = SerialNumberSwitcher()
    
    app.exec_()

    sys.exit(0)

