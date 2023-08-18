/****************************************************************************************************
 *  Laser Barrier Script - Anne-Kristin Lenz - Version 4.0 - 25/06/2020
 ****************************************************************************************************
 * This Arduino sketch is used for a laser barrier triggering measurements with chronos 1.4 highspeed
 * cmaeras.
 * 3 laser transmitters (KY-008) and 3 laser receivers, 2 pairs of those are orthogonal (beams 
 * build cross) to each other and parallel to the third one, are used to detect an object (plastic 
 * bead, water drop) passing through. If the object passes through all three laser beams the
 * measurement will be triggerd. Additionallly the lapse of time between the disruption of the first 
 * and second and between the second and thrid is being displayed on the serial monitor. Only
 * drops that fall through upper and lower laser in certain time will tigger the highspeed cameras 
 * (0-40 ms)
***************************************************************************************************/

/****************************************************************************************************  
 * Declarations
 * 
 * Define all pins of arduino used (Contol output of all receivers and lasers, control input of all
 * lasers with one pin each) and declare variables
 * 
 * Set the threshold time: time between interruption of laser 2 and laser 3 has to be smaller than
 * threshold to trigger the measurement
 * 
 * OLED display declarations and inclusion of libraries
***************************************************************************************************/
// OLED display library inclusions (Install libraries (Sketch> Include Library > Manage Librariers > 
// SSD 1306 by Adafruit  & GFX by Adafruit)) and definitions of screen resolution
#include <Wire.h>
#include <Adafruit_GFX.h>
#include <Adafruit_SSD1306.h>

#define SCREEN_WIDTH 128 // OLED display width, in pixels
#define SCREEN_HEIGHT 64 // OLED display height, in pixels

// Declaration for an SSD1306 display connected to I2C (SDA, SCL pins)
Adafruit_SSD1306 display(SCREEN_WIDTH, SCREEN_HEIGHT, &Wire, -1);

// Potentiometer (analog), control LED, trigger cable to PSV/ High speed camera and reset button 
# define Potentiometer 0
# define Reset 2
# define Control 3
# define Trigger 4


// Receiver Output pins
# define Receiver_1_Out 5
# define Receiver_2_Out 6
# define Receiver_3_Out 7

// Laser Output pins
# define Laser_1_Out 8
# define Laser_2_Out 9
# define Laser_3_Out 10

// Laser Input pins
# define Laser_1_In 11
# define Laser_2_In 12
# define Laser_3_In 13

// Variables for saving time and setting threshold time
unsigned long thresholdTime;  // in us -> Potentiometer for getting threshold
unsigned long time1;
unsigned long time2;
unsigned long time3;
unsigned long dt_1;
unsigned long dt_2;

/**************************************************************************************************** 
 * setup() - function
 *  
 * This function initializes all pins, checks the functioning of all components and gives the signal 
 * (output on serial monitor) that the measurement can be started.
****************************************************************************************************/
void setup() {
  // Start the display
  Serial.begin(9600);
  // Test if display is connected. If not working: output on Serial.monitor and endless for loop.
  // SSD1306_SWITCHCAPVCC = generate display voltage from 3.3V internally
  if(!display.begin(SSD1306_SWITCHCAPVCC, 0x3C)) { // Address 0x3D for 128x64
    Serial.println(F("SSD1306 allocation failed"));
    for(;;);
  }
  // Output on the serial monitor to test, if arduino is connected to computer
  display.clearDisplay();
  display.setTextSize(1);
  display.setTextColor(WHITE);
  display.setCursor(0, 10);
  display.println(F("LASER BARRIER"));
  display.println();
  display.println(F("Starting setup..."));
  display.display();
  delay(2000);

  // initialize all pins
  pinMode(Trigger, OUTPUT);
  pinMode(Control, OUTPUT);
  pinMode(Reset, INPUT);
  pinMode(Receiver_1_Out, INPUT);
  pinMode(Receiver_2_Out, INPUT);
  pinMode(Receiver_3_Out, INPUT);
  pinMode(Laser_1_Out, INPUT);
  pinMode(Laser_2_Out, INPUT);
  pinMode(Laser_3_Out, INPUT);
  pinMode(Laser_1_In, OUTPUT);
  pinMode(Laser_2_In, OUTPUT);
  pinMode(Laser_3_In, OUTPUT);

 // setting up components for measurement
  startSetup();
}

/**************************************************************************************************** 
 * loop() - function
 *  
 * This function checks whether the laser beam is interrupted. If both the laser beams of the upper 
 * laser and the lower 2 lasers are interrupted after each other and the interruption of laser 2 and 
 * laser 3 is within defined threshold time, the measurement is triggered. The time elapsed between 
 * the disruption of the beams is displayed in the serial monitor.
****************************************************************************************************/
void loop() {
  // if laser barrier 1 is interrupted, save time and shut down laser 1
  if(digitalRead(Receiver_1_Out) == LOW){
    time1 = micros();
    digitalWrite(Laser_1_In, LOW);
    
    // wait until laser barrier 2 is interrupted and save time
    while (digitalRead(Receiver_2_Out) == HIGH) {
      }
    time2 = micros();
    digitalWrite(Laser_2_In, LOW);
    
    // wait until laser barrier 3 is interrupted and save time
    while (digitalRead(Receiver_3_Out) == HIGH) {
      }
    time3 = micros();    
    digitalWrite(Laser_3_In, LOW);
    
    // if the time difference between interruption of laser 2 and 3 is shorter than 1 ms 
    // and the time difference between interruption of laser 1 and 2 is shorter than the chosen 
    // thresholdTime the measurement will be triggered (plateau for 1 s) and the time between 
    // the interruptions be displayed
    dt_1 = time2 - time1;
    dt_2 = time3 - time2;
    
    if (dt_2 <= 1000) {
      if (dt_1 < thresholdTime) {
        digitalWrite(Trigger, HIGH);
        delay(1000);
        digitalWrite(Trigger, LOW);
        display.clearDisplay();
        display.setCursor(0, 0);
        display.println(F("Successful trial!"));
        display.println();
        display.print(F("dt_1 (us): "));
        display.println(dt_1);
        display.print(F("dt_2 (us): "));
        display.println(dt_2);
        display.display(); 
       } // if the time difference between interruption of laser 1 and 2 is longer than the chosen 
         // threshold time an error message will be displayed
       else {
        display.clearDisplay();
        display.setCursor(0,0);
        display.println(F("TIME OUT"));
        display.print(F("dt_1 >  thresholdTime, dt_1 (ms):"));
        display.println(dt_1/1000);
        display.display();
       }
    } // if the time difference between interruption of laser 2 and 3 is longer than the theshold
      // time (thresholdTime), an error message will be displayed
    else {
      display.clearDisplay();
      display.setCursor(0, 0);
      display.println(F("TIME OUT"));
      display.print(F("dt_2 > 1 ms, dt_2 (ms):"));
      display.println(dt_2/1000);
      display.display();
    }
    // shut down system in 10 s
    shutdownSetup(0);
    
  } // if laser barrier 2 or 3 get interrupted before laser barrier 1, shutdown of system
  else if(digitalRead(Receiver_2_Out) == LOW || digitalRead(Receiver_3_Out) == LOW) {
    // display error message
    display.clearDisplay();
    display.setCursor(0, 0);
    display.println(F("ERROR: Laser 2 or 3  interrupted before   Laser 1"));
    display.display();
    // shutdown setup instantly
    shutdownSetup(0);
  }
}

/**************************************************************************************************** 
 * startSetup() - function
 * 
 * This function is running at beginning of measurement to turn on lasers and check whether all 
 * lasers and receivers are working. If the setup is ready for the measurement to begin, the line 
 * "Ready for measurement" is displayed in the serial monitor and controll LED is switched off. If 
 * not all lasers and receivers are working, the line "Error in Setup" will be displayed and the first 
 * non-working component detected will be named. The lasers will be switched off after 10s.
****************************************************************************************************/
void startSetup() {
  // Set trigger to LOW
  digitalWrite(Trigger, LOW);
  // Turn on Control LED to show setup starts
  digitalWrite(Control, HIGH);

  // Read Potentiometer and show reading in display.
  thresholdTime = (analogRead(Potentiometer)/ 100 ) * 5000;
  delay(100);
  display.clearDisplay();
  display.setCursor(0, 10);
  display.print(F("Threshold time set to  (ms): "));
  display.println(thresholdTime/1000); 
  display.display();
  
  // Turn on lasers
  digitalWrite(Laser_1_In, HIGH);
  digitalWrite(Laser_2_In, HIGH);
  digitalWrite(Laser_3_In, HIGH);
  delay(1000);
  display.clearDisplay();
  display.setCursor(0, 10);
  
  // Check if all lasers are on and aimed at receivers. If one ore multiple components are not
  // working correctly, first non-working component is detected and an error displayed via the serial
  // monitor. The lasers are shutdown after 10 s.

  int laser[] = {digitalRead(Laser_1_Out),digitalRead(Laser_2_Out),digitalRead(Laser_3_Out)};
  int receiver[] = {digitalRead(Receiver_1_Out), digitalRead(Receiver_2_Out), digitalRead(Receiver_3_Out)};
  int i = 0;
  
  for (i = 0; i <= 2; i++){
    if(laser[i] != HIGH) {
      display.print(F("ERROR: Laser "));
      display.println(i+1); 
      display.display();    
      shutdownSetup(10);
    }
    
    if(receiver[i] != HIGH) {
      display.print(F("ERROR: Receiver "));
      display.println(i+1);
      display.display();
      shutdownSetup(10);
    }
  }
  
  // If all components are working the setup is ready to use
  display.clearDisplay();
  display.setCursor(0, 10);
  display.println(F("Ready for measurement"));
  display.display();
  // Switch off control LED to show measurement can start
  digitalWrite(Control, LOW);
}

/**************************************************************************************************** 
 * shutdownSetup() - function
 * 
 * This function is used to shut down the setup (turing off lasers) in the time defined in the 
 * variable delayed (in s), control LED will be switched on to show measurement is finished or some-
 * thing failed
****************************************************************************************************/
void shutdownSetup(int delayed){
  // The time until lasers are switched off will be displayed
  display.println();
  display.print(F("Laser shutdown in "));
  display.println(delayed);
  display.println(F("Press reset button to start next trial."));
  display.display();
  // Delay in shuting down lasers
  delay(delayed*1000);
  // Turn off lasers
  digitalWrite(Laser_1_In, LOW);
  digitalWrite(Laser_2_In, LOW);
  digitalWrite(Laser_3_In, LOW);
  // wait until reset button is pressed and start setup again
  while (digitalRead(Reset) == HIGH) {
      }
  // clear display before start
  display.clearDisplay();
  display.display();
  startSetup();
}
