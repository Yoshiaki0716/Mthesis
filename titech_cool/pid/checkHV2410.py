import pyvisa


rm = pyvisa.ResourceManager()
print(rm.list_resources())
resource = "ASRL/dev/ttyUSB0::INSTR"
inst = rm.open_resource(resource)
inst.write_termination = '\r\n'
inst.read_termination = '\n'
print(inst.query("*IDN?"))
inst.write("SOUR:FUNC VOLT")
inst.write("SOUR:VOLT 100")
inst.write("SOUR:CURR:PROT 0.01")
inst.write("OUTP ON")
inst.write("MEAS:VOLT?")
print("V:",inst.read())
inst.write("MEAS:CURR?")
print("I:",inst.read())
inst.write("OUTP OFF")
"""
resource = "ASRL/dev/ttyUSB1::INSTR"
inst = rm.open_resource(resource)
inst.write_termination = '\n'
inst.read_termination = '\n'
print(inst.query("*IDN?"))

exit
inst.write("SOUR:FUNC VOLT")
inst.write("SOUR:VOLT 100")
inst.write("SOUR:CURR:PROT 0.01")
inst.write("OUTP ON")
inst.write("OUTP OFF")

inst.write("MEAS:VOLT?")
print("V:",inst.read())
inst.write("MEAS:CURR?")
print("I:",inst.read())

inst.write("OUTP OFF")
"""
