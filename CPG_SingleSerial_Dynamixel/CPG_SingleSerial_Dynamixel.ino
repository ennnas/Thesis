#include <BasicLinearAlgebra.h>     // Libreria per l'algebra matriciale
#include <Math.h>                   
#include <SoftwareSerial.h>         // Libreria per l'utilizzo di porte seriali software
#include <Dynamixel_Serial.h>       // Libreria per il controllo dei servo

#define dt 0.05 // incremento del metodo di eulero
#define L_size 8 // Nt: dimensione della matrice Laplaciana
#define START_POSITION 2048 // Posizione iniziale di 180°
#define OFFSET_SPALLA 285 // Oscillazione angolo spalla di 25°
#define OFFSET_GOMITO 341 // Oscillazione angolo gomito di 30°
#define N_CAMP 40 // Numero di campioni da scartare 
#define GAIT 1 // Tipo di gait da eseguire 1 Camminata Lenta , 2 Trotto , 3 Galoppo

// ######################  SINGOLA SERIALE ###############
#define ID_L1_SPALLA 1
#define ID_L1_GOMITO 2
#define ID_R1_SPALLA 11
#define ID_R1_GOMITO 12
#define ID_L2_SPALLA 21
#define ID_L2_GOMITO 22
#define ID_R2_SPALLA 31
#define ID_R2_GOMITO 32

// software serial pins
SoftwareSerial mySerial(10,11);

// velocità di rotazione del servo [0,4095]
int16_t speed = 0xFFF;

// velocita di comunicazione baudrate
const long unsigned int motor_baudrate = 57600;
const long unsigned int serial_baudrate = 1000000;

// Matrici per il calcolo del CPG
Matrice<1, L_size> x_p, y_p, f, ones;
Matrice<L_size, L_size> L;

// software serial pins
#define SOFT_RX_PIN_L1 6
#define SOFT_TX_PIN_L1 7
#define SOFT_RX_PIN_L2 8
#define SOFT_TX_PIN_L2 9
#define SOFT_RX_PIN_R1 10
#define SOFT_TX_PIN_R1 11
#define SOFT_RX_PIN_R2 12
#define SOFT_TX_PIN_R2 13

// Spalla -> angoli[0] , Gomito -> angoli[1]
float angoli[4][2]; 

// Matrice Laplaciana per slow walk gait
float Ls[L_size][L_size] = {
  {3, 0, 0, -1, 0, 1, 1, 0,},
  {0, 3, 1, 0, -1, 0, 0, 1,},
  {0, 1, 3, 0, 1, 0, 0, -1,},
  { -1, 0, 0, 3, 0, 1, 1, 0,},
  {0, -1, 1, 0, 3, 0, 0, 1,},
  {1, 0, 0, 1, 0, 3, -1, 0,},
  {1, 0, 0, 1, 0, -1, 3, 0,},
  {0, 1, -1, 0, 1, 0, 0, 3,}
};

// Matrice Laplaciana per il trotto
float Lt[L_size][L_size] = {
  {3,0,1,0,1,0,-1,0,},
  {0,3,0,1,0,1,0,-1,},
  {1,0,3,0,-1,0,1,0,},
  {0,1,0,3,0,-1,0,1,},
  {1,0,-1,0,3,0,1,0,},
  {0,1,0,-1,0,3,0,1,},
  {-1,0,1,0,1,0,3,0,},
  {0,-1,0,1,0,1,0,3,}
};

// Matrice Laplaciana per il galoppo
float Lg[L_size][L_size] = {
  {3.00,0.00,-0.81,-0.59,1.00,-0.00,0.81,0.59,},
  {0.00,3.00,0.59,-0.81,0.00,1.00,-0.59,0.81,},
  {-0.81,0.59,3.00,0.00,0.81,-0.59,1.00,-0.00,},
  {-0.59,-0.81,0.00,3.00,0.59,0.81,0.00,1.00,},
  {1.00,0.00,0.81,0.59,3.00,0.00,-0.81,-0.59,},
  {-0.00,1.00,-0.59,0.81,0.00,3.00,0.59,-0.81,},
  {0.81,-0.59,1.00,0.00,-0.81,0.59,3.00,0.00,},
  {0.59,0.81,-0.00,1.00,-0.59,-0.81,0.00,3.00,}

};

float eq1(float x, float y1, float y2) {
  /* Matlab code
    function z = eq1(x,y)
    z = -x(1)+1.7*y(1)-y(2);
  */
  float z = -x + 1.7 * y1 - y2;
  return z;
}

float eq2(float x, float y1, float y2) {
  /* Matlab code
    function z = eq2(x,y)
    z = -x(2)+1.7*y(2)+y(1);
  */
  float z = -x + 1.7 * y2 + y1;
  return z;
}

Matrice<1, L_size> func(Matrice<1, L_size> f, Matrice<1, L_size> x) {
  /* Matlab code
    function z = func(f,x,L)
    z = f.' - 3*L*(x.');
  */
  float k = 3.0;
  Matrice<L_size,1> z;
  z = ~f - (L * (~x)) * k;
  return ~z; // Matrice<1,8>
}

// Converto gli angoli generati dal CPG in angoli per il Dynamixel
void raw_angles(int i,float x) {
  angoli[i][0] = START_POSITION - (x / 2.3) * OFFSET_SPALLA;
  angoli[i][1] = START_POSITION + (x / 2.3) * OFFSET_GOMITO;
}

// Imposta il laplaciano per le andature
void setLaplacian(int val){
  switch (val){
    case 1:
      L = Ls;
      break;
    case 2:
      L = Lt;
      break;
    case 3:
      L = Lg;
      break;
  }
}
void setup() 
{
  Serial.begin(serial_baudrate);
  // SETUP DEL CPG
  setLaplacian(GAIT);
  x_p.Init();
  y_p.Init();
  f.Fill(0);
  ones.Fill(1);
  mySerial.begin(motor_baudrate);

  Dynamixel.begin(mySerial);                 // Set Ardiuno Serial speed to factory default speed of 57600
  for(int i=0; i<40;i+=10){
    Dynamixel.setMode(i+1,SERVO,0x000,speed);              // Turn mode to SERVO, must be WHEEL if using wheel mode
    Dynamixel.setMode(i+2,SERVO,0x000,speed);              // Turn mode to SERVO, must be WHEEL if using wheel mode
  }
  
  // SETUP DELLA ZAMPA, MESSA IN POSIZIONE
  for(int i=0; i<40;i+=10){
    Dynamixel.setMaxTorque(i+1, 0x2FF);
    Dynamixel.servo(i+1,START_POSITION,0x100);
    Dynamixel.setBaudRate(i+1,motor_baudrate);           // Set Dynamixel to new serial speed 
    Dynamixel.setMaxTorque(i+2, 0x2FF);
    Dynamixel.servo(i+2,START_POSITION,0x100);
    Dynamixel.setBaudRate(i+2,motor_baudrate);           // Set Dynamixel to new serial speed 
  }

  Serial.println("Setup Fatto");
  delay(3000);
}

void loop() {
    int start_time = millis();
    // Avanzo di N campioni prima di dare l'angolo al motore
    for(int i=0; i<N_CAMP; i++){
      // calcolo la dinamica del sistema -> f(x,t)
          for (int j = 0; j < L_size; j += 2) {
                f(0, j) = eq1(x_p(0, j), y_p(0, j), y_p(0, j + 1));
                f(0, j + 1) = eq2(x_p(0, j + 1), y_p(0, j), y_p(0, j + 1));
          }
          
          // matrici temporane per il salvataggio dello stato successivo
          Matrice<1, L_size> x_temp, y_temp;
        
          x_temp = x_p + (func(f, x_p)) * dt;
          y_temp = ((x_p + ones).Absolute() - (x_p - ones).Absolute()) * 0.5;
          x_p = x_temp;
          y_p = y_temp;
    }
    //int finish_time = millis() - start_time;
    // converto l'angolo generato in raw
    for(int j=0;j<4;j++){
          raw_angles(j,x_p(0,j*2));
    }
    //Serial << angoli[0][0] << "\n";
    for(int i=0; i<40;i+=10){
          Dynamixel.servo(i+1,angoli[i/10][0],0x100);
          Dynamixel.servo(i+2,angoli[i/10][1],0x100);
          //Serial << angoli[i][0] << " " << angoli[i][1] << "\n";
    }
    int finish_time = millis() - start_time;
    Serial << finish_time << "\n";
}

