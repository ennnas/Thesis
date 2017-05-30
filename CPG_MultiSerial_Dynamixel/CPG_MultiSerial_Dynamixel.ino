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
#define N_Legs 4 // Numero di zampe del robot
// ###################### SERIALI MULTIPLE ###############

// Definizione delle interfaccie Dynamixel
DynamixelClass L1;
DynamixelClass R1;
DynamixelClass L2;
DynamixelClass R2;
DynamixelClass Legs[N_Legs] = {L1,R1,L2,R2};

// software serial pins
#define SOFT_TX_PIN_L1 8
#define SOFT_TX_PIN_R1 9
#define SOFT_TX_PIN_L2 10
#define SOFT_TX_PIN_R2 11

// Definizione delle seriali software (RX,TX)
SoftwareSerial mySerial[N_Legs] ={SoftwareSerial(0,SOFT_TX_PIN_L1),
                                  SoftwareSerial(0,SOFT_TX_PIN_R1),
                                  SoftwareSerial(0,SOFT_TX_PIN_L2),
                                  SoftwareSerial(0,SOFT_TX_PIN_R2)};

// id dei motori (devono essere sequenziali)
#define ID_SPALLA 1
#define ID_GOMITO 2

// velocità di rotazione del servo [0,4095]
int16_t speed = 0xFFF;

// velocita di comunicazione baudrate
const long unsigned int motor_baudrate = 1000000;

// Matrici per il calcolo del CPG
Matrice<1, L_size> x_p, y_p, f, ones;
Matrice<L_size, L_size> L;


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
  Serial.begin(motor_baudrate);
  // SETUP DEL CPG
  setLaplacian(GAIT);
  x_p.Init();
  y_p.Init();
  f.Fill(0);
  ones.Fill(1);
  //  Setup delle seriali
  for(int i=0;i<N_Legs;i++){
    mySerial[i].begin(motor_baudrate);
    Legs[i].begin(mySerial[i]);
  }

  for(int i=0; i<4;i++){
    Legs[i].setMode(ID_SPALLA,SERVO,0x000,speed);              // Turn mode to SERVO, must be WHEEL if using wheel mode
    Legs[i].setMode(ID_GOMITO,SERVO,0x000,speed);              // Turn mode to SERVO, must be WHEEL if using wheel mode
  }
  
  // SETUP DELLA ZAMPA, MESSA IN POSIZIONE
  for(int i=0; i<4;i++){
    Legs[i].setMaxTorque(ID_SPALLA, 0x2FF);
    Legs[i].servo(ID_SPALLA,START_POSITION,0x100);
    Legs[i].setBaudRate(ID_SPALLA,motor_baudrate);           // Set Dynamixel to new serial speed 
    Legs[i].setMaxTorque(ID_GOMITO, 0x2FF);
    Legs[i].servo(ID_GOMITO,START_POSITION,0x100);
    Legs[i].setBaudRate(ID_GOMITO,motor_baudrate);           // Set Dynamixel to new serial speed 
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
    Serial << angoli[0][0] << "\n";
    for(int i=0; i<4;i++){
          Legs[i].servo(ID_SPALLA,angoli[i/10][0],0x100);
          Legs[i].servo(ID_GOMITO,angoli[i/10][1],0x100);
          //Serial << angoli[i][0] << " " << angoli[i][1] << "\n";
    }
    //int finish_time = millis() - start_time;
    //Serial << finish_time << "\n";
}
