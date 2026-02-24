/*

*/
int i;
char key;
const int N=100;
String s[N];
int inPin1=11;
int inPin2=12;
int stat=0;
int t,u;

void setup() {
  Serial.begin( 9600 );     // シリアル通信を初期化する。通信速度は9600bps
  //Serial.println( "Hello Arduino!" );   // 最初に1回だけメッセージを表示する
  pinMode(13, OUTPUT);    
  pinMode(inPin1,OUTPUT);
  pinMode(inPin2,OUTPUT);
  pinMode(10, OUTPUT);
  digitalWrite(13,LOW);
  
}

void loop() {
  String ss;
  String str;
  
  // 受信データがあった時だけ、処理を行う
  while( Serial.available() ) {       // 受信データがあるか？
    key = Serial.read();            // 1文字だけ読み込む
    s[i]=String(key);
    str.concat(s[i]);
    if(i>N) break;
    i++;
  }
  str.trim();//remove last "\n"

  if(i>0){
    //Serial.write( key );            // 1文字送信。受信データをそのまま送り返す。
    //Serial.println(str.charAt(0));
    //Serial.println( str );// 1文字送信。受信データをそのまま送り返す。
    //Serial.println(str.compareTo("OFF"));
    /*
    Serial.println(str.substring(0)); 
    Serial.println(OFF.substring(0)); 
    Serial.println(str.length()); 
    Serial.println(OFF.length()); 
    //Serial.write( key );            // 1文字送信。受信データをそのまま送り返す.
    */
  }
  i=0;
  
  
  if (str.compareTo("ON")==0){
      digitalWrite(13,HIGH);
      digitalWrite(inPin1, HIGH);   // turn the LED on (HIGH is the voltage level)
      delay(100);
      digitalWrite(inPin2, HIGH);
      delay(100);
      //Serial.println("ArduinoON");
      stat=1;
  }
  if (str.compareTo("OFF")==0){
      digitalWrite(13,LOW);
      digitalWrite(inPin1, LOW);   // turn the LED on (HIGH is the voltage level)
      delay(100);
      digitalWrite(inPin2, LOW);
      delay(100);
      //Serial.println("ArduinoOFF");
      stat=0;
  }
  if (str.compareTo("UNLOCK")==0){
      digitalWrite(13,LOW);
      digitalWrite(10, HIGH);   // turn the LED on (HIGH is the voltage level)
      delay(50);
      digitalWrite(10, LOW);
      //Serial.println("UNLOCKED");
      stat=0;
  }
  //Pin status query
  if (str.compareTo("QUERY")==0){
      u=digitalRead(inPin1);
      t=digitalRead(inPin2);
      Serial.println(2*u+t,BIN);
  }
  Serial.flush();
  delay(100);
  
  //check the status of pins
  /*
  if(stat==0){
    if(digitalRead(inPin1)==HIGH || digitalRead(inPin2)==HIGH){
      Serial.println("==Arduino Error==");
      Serial.println("==Caution! Relay Switch doesn't get work! ==");
      Serial.println("Please check the Relay Box");    
    
    } 
  }
  if(stat==1){
    if(digitalRead(inPin1)==LOW || digitalRead(inPin2)==LOW){
      Serial.println("==Arduino Error==");
      Serial.println("==Caution! Relay Switch doesn't get work! ==");
      Serial.println("Please check the Relay Box");    
    } 
  }
  */
}
