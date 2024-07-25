// Motor pins   , driving according to Joshua's code

// Motor pins
const int M1 = 9;
const int E1 = 8;
const int M2 = 11;
const int E2 = 10;

// Number of steps per revolution
const int step_rev = 50; 
const int step45 = 6;  // 45 degrees steps
const int step25 = step_rev /50; // 5 degrees steps
const int stepDelay = 10;

// Switch pins
const int switchPinA0 = A0;
const int switchPinA1 = A1;
const int reset_pin = A2;
const int controlPinA3 = A3;

int curr_pos = 0;  // Track the current step position

void setup() {
  // Set the pin modes
  pinMode(M1, OUTPUT);
  pinMode(E1, OUTPUT);
  pinMode(M2, OUTPUT);
  pinMode(E2, OUTPUT);

  Serial.begin(9600);

  // Initialize the switch pins as input
  pinMode(switchPinA0, INPUT);
  pinMode(switchPinA1, INPUT);
  pinMode(reset_pin, INPUT);
  pinMode(controlPinA3, INPUT);


}

void loop() {
  // Read the switch states
  bool switchStateA0 = digitalRead(switchPinA0);
  bool switchStateA1 = digitalRead(switchPinA1);
  bool resetButtonState = digitalRead(reset_pin);
  int control_vol = analogRead(controlPinA3);
  float vol_converted = control_vol * (5.0 / 1023.0); // Convert analog reading to voltage

  Serial.print("Switch State A0: ");
  Serial.print(switchStateA0);
  Serial.print(" Switch State A1: ");
  Serial.print(switchStateA1);
  Serial.print(" Reset Button State: ");
  Serial.print(resetButtonState);
  Serial.print(" Control Voltage: ");
  Serial.println(vol_converted);
  Serial.print("Current Position: ");
  Serial.println(curr_pos);

  // Primary logic for switch control
  if (resetButtonState == LOW) {  
    if (switchStateA0 == HIGH && switchStateA1 == LOW && curr_pos == 0) {
      rotateMotor(step45, true);  // 45 degrees clockwise
      curr_pos = 1;  // Update the current position
    } else if (switchStateA1 == HIGH && switchStateA0 == LOW && curr_pos == 1) {
      rotateMotor(step45, false);  // 45 degrees counterclockwise
      curr_pos = 0;  // Update the current position
    } else if (switchStateA0 == LOW && switchStateA1 == LOW) {
      if (vol_converted > 2.0 && curr_pos == 0) {
        rotateMotor(step45, true);  // 45 degrees clockwise
        curr_pos = 1;  // Update the current position
      } else if (vol_converted <= 2.0 && curr_pos == 1) {
        rotateMotor(step45, false);  // 45 degrees counterclockwise
        curr_pos = 0;  // Update the current position
      }
    }
  }

  // Rotate the motor 2.5 degrees when the reset button is pressed
  if (resetButtonState == HIGH) {
    rotateMotor(step25, true);  // Rotate 2.5 degrees
  }

  delay(100);
}

void rotateMotor(int steps, bool clockwise) {
  for (int i = 0; i < steps; i++) {
    stepMotor(clockwise);
  }
}

void stepMotor(bool clockwise) {
  if (clockwise) {
    step1();
    delay(stepDelay); // Adjust the delay to control the speed
    step2();
    delay(stepDelay);
    step3();
    delay(stepDelay);
    step4();
    delay(stepDelay);
  } else {
    step4();
    delay(stepDelay);
    step3();
    delay(stepDelay);
    step2();
    delay(stepDelay);
    step1();
    delay(stepDelay);
  }
}

void step1() {
  digitalWrite(E1, HIGH);
  digitalWrite(M1, HIGH);
  digitalWrite(E2, LOW);
  digitalWrite(M2, LOW);
  
}

void step2() {
  digitalWrite(E1, LOW);
  digitalWrite(M1, LOW);
  digitalWrite(E2, HIGH);
  digitalWrite(M2, HIGH);
}

void step3() {
  digitalWrite(E1, HIGH);
  digitalWrite(M1, LOW);
  digitalWrite(E2, LOW);
  digitalWrite(M2, LOW);
}

void step4() {
  digitalWrite(E1, LOW);
  digitalWrite(M1, HIGH);
  digitalWrite(E2, HIGH);
  digitalWrite(M2, LOW);
}

