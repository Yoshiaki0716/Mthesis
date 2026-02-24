#include <Wire.h>
#include <SPI.h>
#include "Adafruit_MAX31855.h"
#include "SHTSensor.h"


#define Resi1_NTC 10000.0
#define Vin_NTC   5.15
#define B_NTC     3435.0
#define Resi0_NTC 10000
#define Abs_Temp  273.15
#define Base_Temp 25.0

#define SteinHart_A 0.8676453371787721e-3
#define SteinHart_B 2.541035850140508e-4
#define SteinHart_C 1.868520310774293e-7


#define USE_HYT271 0


SHTSensor sht;

// Default connection is using software SPI, but comment and uncomment one of
// the two examples below to switch between software SPI and hardware SPI:

// Example creating a thermocouple instance with software SPI on any three
// digital IO pins.

#define MAXDO         14
#define MAXCLK        15
#define MAXCS_Chiller 16
#define MAXCS_Sink    3
#define MAXCS_Case    2
#define MAXCS_Head    4

// initialize the Thermocouple
//Adafruit_MAX31855 thermocouple(MAXCLK, MAXCS, MAXDO);
Adafruit_MAX31855 MAXChiller = Adafruit_MAX31855(MAXCLK, MAXCS_Chiller, MAXDO);
Adafruit_MAX31855 MAXSink    = Adafruit_MAX31855(MAXCLK, MAXCS_Sink, MAXDO);
Adafruit_MAX31855 MAXCase    = Adafruit_MAX31855(MAXCLK, MAXCS_Case, MAXDO);
Adafruit_MAX31855 MAXHead    = Adafruit_MAX31855(MAXCLK, MAXCS_Head, MAXDO);


#define Heater         5
#define Peltier        6
#define Chiller        7
#define PelPlus        8
#define PelMinus       9
#define LowVoltage     10
#define HighVoltage    11
#define CoolingBoxLock 12
#define LockStatus1    18
#define LockStatus2    17
#define ChillerStatus  53


#define INTERLOCK_NORMAL             0
#define INTERLOCK_HIGHTEMP_LV        1
#define INTERLOCK_HIGHTEMP_PELTIER   2
#define INTERLOCK_HIGHTEMP_IDLE      3
#define INTERLOCK_HIGHDEW_PELCHILLER 4
#define INTERLOCK_HIGHDEW_IDLE       5
#define INTERLOCK_CHILLER_ALERT      6
#define INTERLOCK_NTC_DISONN         7

#define INTERLOCK_MAX_TEMP 45 // Cooling Head Max temperature
#define HEATER_MAX_TEMP    25 // Case Max temperature

#define CHILLER_STATE_NORMAL          0
#define CHILLER_STATE_ALERTED         1
#define CHILLER_STATE_BEING_RECOVERED 2


unsigned long OverFlowTime = -1;

// corresponds to ch [5-12] above, 
// heater (ch.5), peltier lines (ch.8&9) must be off by default
int swstate[] = {LOW,HIGH,HIGH,LOW,LOW,HIGH,HIGH,LOW};
int lockState1 = LOW;
int lockState2 = LOW;
int chillerState = LOW;

int chillerAlertState = 0;

double Temp_NTC;
double Temp_Chiller;
double Temp_Case;
double Temp_Sink;
double Temp_Head; 
double Temp_CarrierEnv;
double CarrierEnv_RH;
double CarrierEnv_DP;


int gInterlockStatus = 0;  // Interlock Status


unsigned long millis();
unsigned long Start, End;

unsigned long elapse(unsigned long Start, unsigned long End){
  if(Start > End){
    return ((End - Start) + OverFlowTime);
  }
  else{
    return (End - Start);
  }
}

String split(String STR){
  return  STR.substring(0, STR.indexOf(" "));
}

String splitInt(String STR){
  return STR.substring(STR.indexOf(" ") + 1 );
}

double readThermocoupleTemp( Adafruit_MAX31855& channel, double ref_value ) {
  double out = ref_value;

  const unsigned int maxTry = 10;
  const double compliance = 5.0;
  
  for( int k=0; k<maxTry; k++ ) {
    double tmp_Temp = channel.readCelsius();
    if( abs( tmp_Temp - ref_value ) > compliance ) {
      if( k == maxTry-1 ) {
        out = tmp_Temp;
      } else {
        delay( 100 );
      }
    } else {
      out = tmp_Temp;
      break;
    }
  }
  
  if( isnan( out ) ) {
    out = ref_value;
  }
  return out;
}


bool isCoverOpen() {
  return ( lockState1 == 0 || lockState2 == 0 );
}


//  #define HYT_ADDR 0x28
#define HYT_ADDR 0x35 

float Value_NTC = 0.0;

void doNTC(){
  float Volt_NTC, Resi_NTC;

  Value_NTC = analogRead(A4);
  
  Volt_NTC = Vin_NTC * (Value_NTC / 1023);
  Resi_NTC = Volt_NTC / (Vin_NTC - Volt_NTC)* Resi1_NTC;

  //Serial.print("NTCResistance:");
  //Serial.print(" ");
  //Serial.println(Resi_NTC);
  //Serial.println(Value_NTC);

  //Temp_NTC = 1/((1/B_NTC) * log(Resi_NTC / Resi0_NTC) + 1/(Base_Temp + Abs_Temp)) - Abs_Temp;

  float logResi = log(Resi_NTC);
  float absTemp_NTC = 1./(SteinHart_A + SteinHart_B * logResi + SteinHart_C * logResi * logResi * logResi);
  Temp_NTC = absTemp_NTC-Abs_Temp;
  //Temp_NTC = -9999;

  //Serial.print("NTCTemperature:");
  //Serial.print(" ");
  //Serial.println(Temp_NTC);
}


void doMAX31855() {
  Temp_Head    = readThermocoupleTemp( MAXHead, Temp_Head );
  Temp_Chiller = readThermocoupleTemp( MAXChiller, Temp_Chiller );
  Temp_Sink    = readThermocoupleTemp( MAXSink, Temp_Sink );
  Temp_Case    = readThermocoupleTemp( MAXCase, Temp_Case );
  
  //Serial.print("Chiller temp = ");
  Serial.print(Temp_Chiller);
  Serial.print(", ");
  //Serial.print("Case temp = ");
  Serial.print(Temp_Case);
  Serial.print(", ");
  //Serial.print("Sink temp = ");
  Serial.print(Temp_Sink);
  Serial.print(", ");
  //Serial.print("Head temp = ");
  Serial.print(Temp_Head);
  Serial.print(", ");

  double internal_chiller = MAXChiller.readInternal();
  double internal_sink    = MAXSink.readInternal();
  double internal_case    = MAXCase.readInternal();
  double internal_head    = MAXHead.readInternal();

  Serial.print(internal_chiller);
  Serial.print(", ");
  Serial.print(internal_case);
  Serial.print(", ");
  Serial.print(internal_sink);
  Serial.print(", ");
  Serial.println(internal_head);
  return;
}

#if USE_HYT271

void doHYT() {
  // All low-level things about HYT271
  

      
  Wire.beginTransmission(HYT_ADDR);   // Begin transmission with given device on I2C bus
  Wire.requestFrom(HYT_ADDR, 4);      // Request 4 bytes 
      
  // Read the bytes if they are available
  // The first two bytes are humidity the last two are temperature
  if(Wire.available() == 4) {    
                      
    int b1 = Wire.read();
    int b2 = Wire.read();
    int b3 = Wire.read();
    int b4 = Wire.read();
        
    Wire.endTransmission();           // End transmission and release I2C bus
        
    // combine humidity bytes and calculate humidity
    int rawCarrierEnv_RH = b1 << 8 | b2;
    // compound bitwise to get 14 bit measurement first two bits
    // are status/stall bit (see intro text)
    rawCarrierEnv_RH =  (rawCarrierEnv_RH &= 0x3FFF);
    CarrierEnv_RH = 100.0 / pow(2,14) * rawCarrierEnv_RH;
        
    // combine temperature bytes and calculate temperature
    b4 = (b4 >> 2); // Mask away 2 least significant bits see HYT 221 doc
    int rawTemperature = b3 << 6 | b4;
    Temp_CarrierEnv = 165.0 / pow(2,14) * rawTemperature - 40;


    double b = 18.678;
    double c = 257.14;
    double d = 234.5;
    double g = log(CarrierEnv_RH/100.) + b * Temp_CarrierEnv / (c + Temp_CarrierEnv);

    CarrierEnv_DP = c * g / (b - g );
        
  }
  else {
    //Serial.println("Not enough bytes available on wire.");
  }
  return;
}

#else

void doSHT(){
  sht.init();
  delay(2);
  sht.setAccuracy(SHTSensor::SHT_ACCURACY_HIGH);
  delay(2);
  sht.readSample();
  Temp_CarrierEnv = sht.getTemperature();
  CarrierEnv_RH = sht.getHumidity();
  
  double A = 5.64;
  double B = 0.861 * Temp_CarrierEnv - 46.4;

  CarrierEnv_DP = A * sqrt(CarrierEnv_RH) + B;
}

#endif



void doRelayOff(String str) {
  if(gInterlockStatus == 0){
    int channel = splitInt(str).toInt();
    if( channel < 1 || channel > 8 ) {
      Serial.println("ERROR [RelayOff]: Relay channel must be [1-8]!");
    }
    digitalWrite(channel + 4, LOW);
    swstate[splitInt(str).toInt() - 1] = LOW;
    Serial.print("Relay channel No.");
    Serial.print(splitInt(str).toInt());
    Serial.println(" has switched to LOW.");
  } else {
    Serial.println("rf rejected due to interlock!");
  }  
}

void doRelayOn( String str ) {
  if(gInterlockStatus == 0){
    int channel = splitInt(str).toInt();
    if( channel < 1 || channel > 8 ) {
      Serial.println("ERROR [RelayOff]: Relay channel must be [1-8]!");
    }
    digitalWrite(channel + 4, HIGH);
    swstate[splitInt(str).toInt() - 1] = HIGH;
    Serial.print("Relay channel No.");
    Serial.print(splitInt(str).toInt());
    Serial.println(" has switched to HIGH.");
  } else {
    Serial.println("rf rejected due to interlock!");
  }
}

void doRelayStates() {
  for( int k=0; k<8; k++) {
    Serial.print( swstate[k] );
    Serial.print(", ");
  }
  lockState1 = digitalRead( LockStatus1 );
  Serial.print( lockState1 );
  Serial.print(", ");
  lockState2 = digitalRead( LockStatus2 );
  Serial.print( lockState2 );
  Serial.print(", ");
  if( chillerState == 0 ) {
    chillerState = digitalRead( ChillerStatus );
  }
  Serial.print( chillerState );
  Serial.println();
}

void doPelPolInverted() {
  if(gInterlockStatus == 0){
    // Here, flip two Peltier pins simultaneously to HIGH
    //PORTB |= 0b00000011;   // If you use Arduino uno
          
    digitalWrite(PelPlus, HIGH);
    digitalWrite(PelMinus, HIGH); //If you use Arduino mega
           
    swstate[PelPlus - Heater] = HIGH;
    swstate[PelMinus - Heater] = HIGH; 
    Serial.println("Peltier polarity is inverted.");
  } else {
    Serial.println("PelPol rejected due to interlock!");
  }
}

void doPelPolNormal() {
  if(gInterlockStatus == 0){
    // Here, flip two Peltier pins simultaneously to LOW
    //PORTB &= ~ 0b00000011; //If you use Arduino uno

    digitalWrite(PelPlus, LOW);
    digitalWrite(PelMinus, LOW); // If you use Arduino mega
        
    swstate[PelPlus - Heater]  = LOW;
    swstate[PelMinus - Heater] = LOW; 
    Serial.println("Peltier polarity is normal.");
  } else {
    Serial.println("PelPol rejected due to interlock!");
  }  
}

void doUnlock() {
  digitalWrite(CoolingBoxLock, HIGH);
  delay(500);
  digitalWrite(CoolingBoxLock, LOW);
  Serial.println("Unlocked the box");
}

void doStatus() {
  //Serial.print("Current situation of the system is ");
  Serial.println(gInterlockStatus);
  //Serial.println(".");

}


void doReset() {
  // int swstate[] = {LOW,HIGH,HIGH,LOW,LOW,HIGH,HIGH,LOW};
  swstate[0] = LOW;
  swstate[1] = HIGH;
  swstate[2] = HIGH;
  swstate[3] = LOW;
  swstate[4] = LOW;
  swstate[5] = HIGH;
  swstate[6] = HIGH;
  swstate[7] = LOW;
  for( int k=0; k<8; k++ ) {
    digitalWrite( k+5, swstate[k] );
  }

  // Here, chillerState has a state transition of
  // 0 --> 1 --> 2 --> 0 --> ...
  // The latch 0 --> 1 is triggered by the chiller's alert signal.
  // The latch 1 --> 2 is set by the 1st reset command by the user.
  // The latch 2 --> 0 is set by the 2nd reset command by the user
  // where the interlock status recovers to Normal (0).
  if( chillerState == CHILLER_STATE_NORMAL ) {
    gInterlockStatus = INTERLOCK_NORMAL;
  } else if( chillerState == CHILLER_STATE_ALERTED ) {
    chillerState = 2;
  } else if( chillerState == CHILLER_STATE_BEING_RECOVERED ) {
    chillerState = 0;
    gInterlockStatus = INTERLOCK_NORMAL;
  }
   
  Serial.println("Interlock has been reset.");
}

/////////////////////////////////////////////////////////////////////////////////////////////
void doInterlockTest( String str){
  if(gInterlockStatus ==0){
    int mode = splitInt(str).toInt();
    if (mode<0 or mode>7){
      Serial.println("Error [InterlockTest]: mode must be [0-7]");
    }

    //INTERLOCK_NORMAL///////////////////////////////////////////
    if(mode == 0){
      Serial.println("Noting happened because INTERLOCK_NORMAL is the state with no abnormality.");
    }

    //INTERLOCK_HIGHTEMP_LV//////////////////////////////////////
    if(mode == 1){
      digitalWrite(LowVoltage,LOW);
      digitalWrite(HighVoltage,LOW);
      swstate[LowVoltage-5]=LOW;
      swstate[HighVoltage-5]=LOW;
      gInterlockStatus = INTERLOCK_HIGHTEMP_LV;
      Serial.println("INTERLOCK_HIGHTEMP_LV became active!");
    }

    //INTERLOCK_HIGHTEMP_PELTIER/////////////////////////////////
    if(mode == 2){
      digitalWrite(Peltier, LOW);
      // digitalWrite(Chiller, LOW); this is also a stupid bug
      swstate[Peltier-5]=LOW;
      gInterlockStatus = INTERLOCK_HIGHTEMP_PELTIER;
      Serial.println("INTERLOCK_HIGHTEMP_PELTIER became active!");
    }

    //INTERLOCK_HIGHTEMP_IDLE////////////////////////////////////
    if(mode == 3){
      digitalWrite(Chiller, LOW);
      swstate[Chiller-5]=LOW;
      gInterlockStatus = INTERLOCK_HIGHTEMP_IDLE;
      Serial.println("INTERLOCK_HIGHTEMP_IDLE became active!");
    }

    //INTERLOCK_HIGHDEW_PELCHILLER///////////////////////////////
    if(mode == 4){
      Serial.println("Nothing happened because INTERLOCK_HIDEW_PELCHILLER is not used now");
    }

    //INTERLOCK_HIGHDEW_IDLE/////////////////////////////////////
    if(mode == 5){
      digitalWrite(LowVoltage,LOW);
      digitalWrite(HighVoltage, LOW);
      swstate[LowVoltage-5]=LOW;
      swstate[HighVoltage-5]=LOW;
      gInterlockStatus = INTERLOCK_HIGHDEW_IDLE;
      Serial.println("INTERLOCK_HIGHDEW_IDLE became active!");
    }

    //INTERLOCK_CHILLER_ALERT////////////////////////////////////
    if(mode == 6){
      digitalWrite(Peltier, LOW);
      digitalWrite(Chiller, LOW);
      swstate[Peltier-5]=LOW;
      swstate[Chiller-5]=LOW;
      gInterlockStatus = INTERLOCK_CHILLER_ALERT;
      Serial.println("INTERLOCK_CHILLER_ALERT became active!");
    }

    //INTERLOCK_NTC_DISONN///////////////////////////////////////////
    if(mode == 7){
      digitalWrite(Peltier, LOW);
      digitalWrite(LowVoltage, LOW);
      digitalWrite(HighVoltage, LOW);
      gInterlockStatus = INTERLOCK_NTC_DISONN;
      Serial.println("INTERLOCK_NTC_DISONN became active!");
    }
  }
  else{Serial.println("[InterlockTest] rejected due to interlock which already exist!");}
}
////////////////////////////////////////////////////////////////////////////////////////////////


void setup() {
  Serial.begin( 9600 );

  Wire.begin();             // Join I2c Bus as master
                            // I2C Pins on the Arduino Uno are A4 for SDA and A5 for SCL

  pinMode(PelPlus,        OUTPUT);
  pinMode(PelMinus,       OUTPUT);
  pinMode(LowVoltage,     OUTPUT);
  pinMode(HighVoltage,    OUTPUT);
  pinMode(Peltier,        OUTPUT);
  pinMode(Chiller,        OUTPUT);
  pinMode(Heater,         OUTPUT);
  pinMode(CoolingBoxLock, OUTPUT);
  pinMode(ChillerStatus,  INPUT);

  for( int k=0; k<8; k++ ) {
    digitalWrite( k+5, swstate[k] );
  }

  sht.init();
  delay(2);
  sht.setAccuracy(SHTSensor::SHT_ACCURACY_HIGH);
  delay(2);
  
  while (!Serial) { delay(1); } // wait for Serial on Leonardo/Zero, etc
}


// Global Serial communication variables
String inputString = "";      // a String to hold incoming data
bool stringComplete = false;  // whether the string is complete


void loop() {
    
  // basic readout test, just print the current temp
  if ( stringComplete ) {
 
    if( String("Temps") == inputString ){
      doMAX31855();
    }

    else if( String("ModuleTemp") == inputString){
      doNTC();
      if( Temp_NTC > -60 ) {
	//Serial.println( analogRead( A4 ) );
	Serial.println(Temp_NTC);
      } else {
	Temp_NTC = -9999;
	Serial.println("IRREGULAR");
      }
    }

    else if( String("CarrierEnv") == inputString){
      
#if USE_HYT271
      doHYT();
#else
      doSHT();
#endif   
      Serial.print(Temp_CarrierEnv, 2);
      Serial.print(", ");
      Serial.print(CarrierEnv_RH, 2);
      Serial.print(", ");
      Serial.println(CarrierEnv_DP, 2);
    }

    else if( String("RelayOff") == split(inputString) ) {
      doRelayOff( inputString );
    }

    else if( String("RelayOn") == split(inputString) ){
      doRelayOn( inputString );
    }

    else if( String("RelayStates") == inputString ){
      doRelayStates();
    }

    else if(String("PelPolNormal") == inputString){
      doPelPolNormal();
    }

    else if(String("PelPolInverted") == inputString){
      doPelPolInverted();
    }

    else if(String("Unlock") == inputString ) {
      doUnlock();
    }

    else if( String("Status") == inputString || String("Stat") == inputString){
      doStatus();
    }

    else if( String("Reset") == inputString ) {
      doReset();
    }   

    else if( String("InterlockTest") == split(inputString) ){
      doInterlockTest( inputString );
    }

    else {
      Serial.println("Invalid command!: '" + inputString + "'");
    }
    
  }

  ///////////////case1 interlock by over-temperature /////////////
  Temp_Head    = readThermocoupleTemp( MAXHead, Temp_Head );
  Temp_Chiller = readThermocoupleTemp( MAXChiller, Temp_Chiller );
  Temp_Sink    = readThermocoupleTemp( MAXSink, Temp_Sink );
  Temp_Case    = readThermocoupleTemp( MAXCase, Temp_Case );    
    
  doNTC();
  if( Temp_NTC < -70 ) {
    Temp_NTC = -9999;
  }

  if( chillerState == 0 ) {
    chillerState = digitalRead( ChillerStatus );
  }
    
#if USE_HYT271
  doHYT();
#else
  doSHT();
#endif
    
  if( Temp_Case > HEATER_MAX_TEMP ) {
    digitalWrite(Heater, LOW);
    swstate[Heater-5] = LOW;
  }
  
  
  if( ( Temp_Head > INTERLOCK_MAX_TEMP || Temp_NTC > INTERLOCK_MAX_TEMP) &&
      gInterlockStatus == INTERLOCK_NORMAL){
      
    // High temperature interlock
    // Switch off module LV and HV
    digitalWrite(LowVoltage, LOW);
    digitalWrite(HighVoltage, LOW);
    swstate[LowVoltage-5] = LOW;
    swstate[HighVoltage-5] = LOW;
    
    gInterlockStatus = INTERLOCK_HIGHTEMP_LV;
    Start = millis();
  }
    
  if(gInterlockStatus == INTERLOCK_HIGHTEMP_LV){
    End = millis();
	
    if( Temp_Case > HEATER_MAX_TEMP ) {
      digitalWrite(Heater, LOW);
      swstate[Heater-5] = LOW;
    }
	
    gInterlockStatus = INTERLOCK_HIGHDEW_IDLE;
  }
    
    
//  // additional protection against high DP
//  if( false /*Temp_Head < 10.0 &&
//	      isCoverOpen() &&
//	      gInterlockStatus == INTERLOCK_NORMAL */) {
//      
//    // High dew point interlock
//    // Switch off Peltier and Chiller
//    digitalWrite( Peltier, LOW );
//    digitalWrite( Chiller, LOW );
//    swstate[Peltier-5] = LOW;
//    swstate[Chiller-5] = LOW;
//
//    // warm up by heater
//    digitalWrite(Heater, HIGH);
//    swstate[Heater-5] = HIGH;
//      
//    gInterlockStatus = INTERLOCK_HIGHDEW_PELCHILLER;
//       
//  }
    
  // chiller trigger alert
  if( chillerState == 1 ) {
    digitalWrite( Peltier, LOW );
    // digitalWrite( Chiller, LOW ); // commented out a stupid bug.
    swstate[Peltier-5] = LOW;
      
    gInterlockStatus = INTERLOCK_HIGHTEMP_PELTIER;

    Start = millis();
  }

  if(gInterlockStatus == INTERLOCK_HIGHTEMP_PELTIER){
    End = millis();
	
    // After two hours, turn off the chiller
    if(elapse(Start, End)/1000 > 7200){
      digitalWrite(Chiller, LOW);
      swstate[Chiller-5] = LOW;
      
      gInterlockStatus = INTERLOCK_HIGHTEMP_IDLE;
    }
  }
    
//#if 0
//  // protection against high DP
//  if( CarrierEnv_DP > Temp_Head - 2.0 &&
//      gInterlockStatus == INTERLOCK_NORMAL ) {
//      
//    // High dew point interlock
//    // Switch off Peltier and Chiller
//    digitalWrite( Peltier, LOW );
//    digitalWrite( Chiller, LOW );
//    swstate[Peltier-5] = LOW;
//    swstate[Chiller-5] = LOW;
//
//    // warm up by heater
//    if( Temp_Case < HEATER_MAX_TEMP ) {
//      digitalWrite(Heater, HIGH);
//      swstate[Heater-5] = HIGH;
//    } else {
//      digitalWrite(Heater, LOW);
//      swstate[Heater-5] = LOW;
//    }
//    //gInterlockStatus = INTERLOCK_HIGHDEW_PELCHILLER;
//       
//  }
//    
//  if ( gInterlockStatus == INTERLOCK_HIGHDEW_PELCHILLER ) {
//    if( Temp_Case > HEATER_MAX_TEMP ) {
//      digitalWrite(Heater, LOW);
//      swstate[Heater-5] = LOW;
//    }
//  }
//
//
//
//  if( gInterlockStatus == INTERLOCK_HIGHDEW_PELCHILLER && 
//      Temp_Head > CarrierEnv_DP + 10) {
//      
//    // Sufficiently high temperature than DP
//    // Switch off Module LV and HV
//    digitalWrite( LowVoltage, LOW );
//    digitalWrite( HighVoltage, LOW );
//    swstate[LowVoltage-5] = LOW;
//    swstate[HighVoltage-5] = LOW;
//	
//    if( Temp_Case > HEATER_MAX_TEMP ) {
//      digitalWrite(Heater, LOW);
//      swstate[Heater-5] = LOW;
//    }
//	
//    // gInterlockStatus = INTERLOCK_HIGHDEW_IDLE;
//  }
//    
//    
//  // additional protection against high DP
//  if( false /*Temp_Head < 10.0 &&
//	      isCoverOpen() &&
//	      gInterlockStatus == INTERLOCK_NORMAL*/ ) {
//      
//    // High dew point interlock
//    // Switch off Peltier and Chiller
//    digitalWrite( Peltier, LOW );
//    digitalWrite( Chiller, LOW );
//    swstate[Peltier-5] = LOW;
//    swstate[Chiller-5] = LOW;
//
//    // warm up by heater
//    digitalWrite(Heater, HIGH);
//    swstate[Heater-5] = HIGH;
//	
//    // gInterlockStatus = INTERLOCK_HIGHDEW_PELCHILLER;
//       
//  }
//#endif
    
//  // chiller trigger alert
//  if( chillerState == 1 ) {
//    digitalWrite( Peltier, LOW );
//    digitalWrite( Chiller, LOW );
//    swstate[Peltier-5] = LOW;
//    swstate[Chiller-5] = LOW;
//
//    gInterlockStatus = INTERLOCK_CHILLER_ALERT;
//  }
    
    
  // Module NTC disconnection while cover is closed
  if( Temp_NTC < -50 &&
      ( lockState1 == 1 || lockState2 == 1 ) &&
      gInterlockStatus == INTERLOCK_NORMAL) {
      
    digitalWrite( Peltier, LOW );
    digitalWrite( LowVoltage, LOW );
    digitalWrite( HighVoltage, LOW );
	
    gInterlockStatus = INTERLOCK_NTC_DISONN;
  }

  if( Temp_Case > HEATER_MAX_TEMP ) {
    digitalWrite(Heater, LOW);
    swstate[Heater-5] = LOW;
  }

  if( stringComplete || inputString == "" ) {
    reset();
  }
  
}



void reset() {
  inputString = "";
  stringComplete = false;
}



void serialEvent() {
  while (Serial.available()) {
    // get the new byte:
    char inChar = (char)Serial.read();
    // add it to the inputString:
    if (inChar == '\n') {
      stringComplete = true;
    } else {
      inputString += inChar;
    }
  }
}
