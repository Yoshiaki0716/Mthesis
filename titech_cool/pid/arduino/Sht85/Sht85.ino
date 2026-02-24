#include <Wire.h>
#include "SHTSensor.h"

SHTSensor sht;
// To use a specific sensor instead of probing the bus use this command:
// SHTSensor sht(SHTSensor::SHT3X);

void setup() {
  // put your setup code here, to run once:

  Wire.begin();
  Serial.begin(9600);
  delay(1000); // let serial console settle

  if (sht.init()) {
      // Serial.print("Temperature\tRelativeHumidity\n");
  } else {
  }
  sht.setAccuracy(SHTSensor::SHT_ACCURACY_HIGH); // only supported by SHT3x

}

void loop() {
  // put your main code here, to run repeatedly:

  digitalWrite(13, HIGH);

  while (Serial.available()) {
    char c = Serial.read();
    if( c == 't' ) {
      if (sht.readSample()) {

        double T = sht.getTemperature();
        Serial.print(T, 7);
        Serial.print("\n");

      } else {
        Serial.print("Error in readSample()\n");
      }
    } else if( c == 'h' ) {
      if (sht.readSample()) {

        double H = sht.getHumidity();
        Serial.print(H, 7);
        Serial.print("\n");

      } else {
        Serial.print("Error in readSample()\n");
      }      
    } else if( c == 'a' ) {
      Serial.print("Active!\n");
    }
    
    delay(100);
    digitalWrite(13, LOW);
    delay(100);
  }


}
