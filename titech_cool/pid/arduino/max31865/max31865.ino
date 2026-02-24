#include <Wire.h>
#include <SPI.h>
#include <Adafruit_MAX31865.h>

// Use software SPI: CS, DI, DO, CLK
#define MAX_CS 13
#define MAX_DI 12
#define MAX_DO 11
#define MAX_CLK 10
Adafruit_MAX31865 max = Adafruit_MAX31865(MAX_CS, MAX_DI, MAX_DO, MAX_CLK);
float RREF = 1000.;

int inPin =7;
int outPin = 13;

void setup() {
  // put your setup code here, to run once:
  Wire.begin();
  Serial.begin(9600);
  max.begin(MAX31865_2WIRE);
  while (Serial.available () == 0) {}
}

void loop() {
  while (Serial.available()) {
    char c = Serial.read();
    int rtd = max.readRTD();
    //Serial.print("RTD value: ");
    //Serial.println(rtd);
    float ratio = float(rtd) / 32768.0;
    float res_raw = RREF*ratio;
    float res_calib = res_raw + ( 0.33433 - 0.00622344 * res_raw );
    //Serial.print("Ratio = "); Serial.println(ratio,8);
    //Serial.print("Resistance = ");
    Serial.println(res_calib,8);
  }
  // put your main code here, to run repeatedly:
}
