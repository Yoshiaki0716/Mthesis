#define Resi1_NTC 47000.0
#define Vin_NTC 5.0
#define B_NTC 3435.0
#define Resi0_NTC 10000
#define Abs_Temp 273.15
#define Base_Temp 25.0

#define SteinHart_A 0.8676453371787721e-3
#define SteinHart_B 2.541035850140508e-4
#define SteinHart_C 1.868520310774293e-7

double Temp_NTC;

void setup(){
  Serial.begin(9600);
}

void doNTC(){
  float Value_NTC, Volt_NTC, Resi_NTC;

  Value_NTC = analogRead(A0);
  
  Volt_NTC = Vin_NTC*1000*Value_NTC / 1023;
  Resi_NTC = Volt_NTC*Resi1_NTC / (Vin_NTC*1000 - Volt_NTC);

  //Serial.print("NTCResistance:");
  //Serial.print(" ");
  //Serial.println(Resi_NTC);

  //  Temp_NTC = 1/((1/B_NTC) * log(Resi_NTC / Resi0_NTC) + 1/(Base_Temp + Abs_Temp));

  float logResi = log(Resi_NTC);
  float absTemp_NTC = 1./(SteinHart_A + SteinHart_B * logResi + SteinHart_C * logResi * logResi * logResi);
  Temp_NTC = absTemp_NTC-Abs_Temp;


  //Serial.print("NTCTemperature:");
  //Serial.print(" ");
  //Serial.println(Temp_NTC);
}


void loop(){
  if( Serial.available() ){
    String str = Serial.readStringUntil('\n');


    if( String("ModuleTemp") == str){
      doNTC();
      if( Temp_NTC > -70 ) {
        Serial.println(Temp_NTC);
      } else {
        Serial.println("IRREGULAR");
      }
    } else {
      Serial.println("Invalid Command!");
    }
  }
  delay(100);
}
