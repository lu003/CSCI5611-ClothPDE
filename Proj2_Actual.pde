//PDEs and Integration
//CSCI 5611 Swinging Rope [Exercise]
//Stephen J. Guy <sjguy@umn.edu>

//NOTE: The simulation starts paused, press "space" to unpause

//TODO:
//  1. The rope moves very slowly now, this is because the timestep is 1/20 of realtime
//      a. Make the timestep realtime (20 times faster than the inital code), what happens?
//      b. Call the small 1/20th timestep update, 20 times each frame (in a for loop) -- why is it different?
//  2. When the rope hanging down fully the spacing between the links is not equal, even though they
//     where initalized with an even spacing between each node. What is this?
//      - If this is a bug, fix the corisponding code
//      - If this why is not a bug, explain why this is the expected behavior
//  3. By default, the rope starts vertically. Change initScene() so it starts at an angle. The rope should
//     then swing backc and forth.
//  4. Try changing the mass and the k value. How do they interact wich each other?
//  5. Set the kv value very low, does the rope bounce a lot? What about a high kv (not too high)?
//     Why doesn't the rope stop swinging at high values of kv?
//  6. Add friction so that the rope eventually stops. An easy friction model is a scaled force 
//     in the opposite direction of a nodes current velocity. 

//Challenge:
//  - Set the top of the rope to be wherever the user’s mouse is, and allow the user to drag the rope around the scene.
//  - Keep the top of the rope fixed, but allow the user to click and drag one of the balls of the rope to move it around.
//  - Place a medium-sized, static 2D ball in the scene, have the nodes on the rope experience a “bounce” force if they collide with this ball.


//Create Window
String windowTitle = "Swinging Rope";
void setup() {
  size(800, 800, P3D);
  camera = new Camera();
  surface.setTitle(windowTitle);
  initScene();
}

class Camera
{
  Camera()
  {
    position      = new PVector( -500, 0, 500 ); // initial position
    theta         = -3.14/4; // rotation around Y axis. Starts with forward direction as ( 0, 0, -1 )
    phi           = 0; // rotation around X axis. Starts with up direction as ( 0, 1, 0 )
    moveSpeed     = 50;
    turnSpeed     = 1.57; // radians/sec
    boostSpeed    = 10;  // extra speed boost for when you press shift
    
    // dont need to change these
    shiftPressed = false;
    negativeMovement = new PVector( 0, 0, 0 );
    positiveMovement = new PVector( 0, 0, 0 );
    negativeTurn     = new PVector( 0, 0 ); // .x for theta, .y for phi
    positiveTurn     = new PVector( 0, 0 );
    fovy             = PI / 4;
    aspectRatio      = width / (float) height;
    nearPlane        = 0.1;
    farPlane         = 10000;
  }
  
  void Update(float dt)
  {
    theta += turnSpeed * ( negativeTurn.x + positiveTurn.x)*dt;
    
    // cap the rotation about the X axis to be less than 90 degrees to avoid gimble lock
    float maxAngleInRadians = 85 * PI / 180;
    phi = min( maxAngleInRadians, max( -maxAngleInRadians, phi + turnSpeed * ( negativeTurn.y + positiveTurn.y ) * dt ) );
    
    // re-orienting the angles to match the wikipedia formulas: https://en.wikipedia.org/wiki/Spherical_coordinate_system
    // except that their theta and phi are named opposite
    float t = theta + PI / 2;
    float p = phi + PI / 2;
    PVector forwardDir = new PVector( sin( p ) * cos( t ),   cos( p ),   -sin( p ) * sin ( t ) );
    PVector upDir      = new PVector( sin( phi ) * cos( t ), cos( phi ), -sin( t ) * sin( phi ) );
    PVector rightDir   = new PVector( cos( theta ), 0, -sin( theta ) );
    PVector velocity   = new PVector( negativeMovement.x + positiveMovement.x, negativeMovement.y + positiveMovement.y, negativeMovement.z + positiveMovement.z );
    position.add( PVector.mult( forwardDir, moveSpeed * velocity.z * dt ) );
    position.add( PVector.mult( upDir,      moveSpeed * velocity.y * dt ) );
    position.add( PVector.mult( rightDir,   moveSpeed * velocity.x * dt ) );
    
    aspectRatio = width / (float) height;
    perspective( fovy, aspectRatio, nearPlane, farPlane );
    camera( position.x, position.y, position.z,
            position.x + forwardDir.x, position.y + forwardDir.y, position.z + forwardDir.z,
            upDir.x, upDir.y, upDir.z );
  }
  
  // only need to change if you want difrent keys for the controls
  void HandleKeyPressed()
  {
    if ( key == 'w' || key == 'W' ) positiveMovement.z = 1;
    if ( key == 's' || key == 'S' ) negativeMovement.z = -1;
    if ( key == 'a' || key == 'A' ) negativeMovement.x = -1;
    if ( key == 'd' || key == 'D' ) positiveMovement.x = 1;
    if ( key == 'q' || key == 'Q' ) positiveMovement.y = 1;
    if ( key == 'e' || key == 'E' ) negativeMovement.y = -1;
    if ( key == 'i'  ) theta += 0.01;
    if ( key == 'k'  ) theta += -0.01;
    
    
    
    if ( keyCode == LEFT )  negativeTurn.x = 1;
    if ( keyCode == RIGHT ) positiveTurn.x = -0.5;
    if ( keyCode == UP )    positiveTurn.y = 0.5;
    if ( keyCode == DOWN )  negativeTurn.y = -1;
    
    if ( keyCode == SHIFT ) shiftPressed = true; 
    if (shiftPressed){
      positiveMovement.mult(boostSpeed);
      negativeMovement.mult(boostSpeed);
    }
    
  }
  
  // only need to change if you want difrent keys for the controls
  void HandleKeyReleased()
  {
    if ( key == 'w' || key == 'W' ) positiveMovement.z = 0;
    if ( key == 'q' || key == 'Q' ) positiveMovement.y = 0;
    if ( key == 'd' || key == 'D' ) positiveMovement.x = 0;
    if ( key == 'a' || key == 'A' ) negativeMovement.x = 0;
    if ( key == 's' || key == 'S' ) negativeMovement.z = 0;
    if ( key == 'e' || key == 'E' ) negativeMovement.y = 0;
    
    if ( keyCode == LEFT  ) negativeTurn.x = 0;
    if ( keyCode == RIGHT ) positiveTurn.x = 0;
    if ( keyCode == UP    ) positiveTurn.y = 0;
    if ( keyCode == DOWN  ) negativeTurn.y = 0;
    
    if ( keyCode == SHIFT ){
      shiftPressed = false;
      positiveMovement.mult(1.0/boostSpeed);
      negativeMovement.mult(1.0/boostSpeed);
    }
  }
  
  // only necessary to change if you want different start position, orientation, or speeds
  PVector position;
  float theta;
  float phi;
  float moveSpeed;
  float turnSpeed;
  float boostSpeed;
  
  // probably don't need / want to change any of the below variables
  float fovy;
  float aspectRatio;
  float nearPlane;
  float farPlane;  
  PVector negativeMovement;
  PVector positiveMovement;
  PVector negativeTurn;
  PVector positiveTurn;
  boolean shiftPressed;
};




// ----------- Example using Camera class -------------------- //
Camera camera;



void keyReleased()
{
  camera.HandleKeyReleased();
}

//Simulation Parameters
float floor = 500;
Vec3 gravity = new Vec3(0,400,0);
float radius = 2;
Vec3 stringTop = new Vec3(0,0,0);
float restLen = 10;
float mass = 1.0; //TRY-IT: How does changing mass affect resting length of the rope?
float k = 200; //TRY-IT: How does changing k affect resting length of the rope?
float kv = 30; //TRY-IT: How big can you make kv?
float kf = 0.8; 
Vec3 Obstacle=new Vec3(100,100,100);
float obstacleR=30;
//Initial positions and velocities of masses
static int maxNodes = 100;
//Vec3 pos[] = new Vec3[maxNodes];
//Vec3 vel[] = new Vec3[maxNodes];
//Vec3 acc[] = new Vec3[maxNodes];

Vec3 pos[][]= new Vec3[maxNodes][maxNodes];
Vec3 vel[][]= new Vec3[maxNodes][maxNodes];
Vec3 acc[][]= new Vec3[maxNodes][maxNodes];

int numNodes = 10;
int numThreads=15;

void initScene(){
  for(int j=0;j<numThreads;j++)
{  for (int i = 0; i < numNodes; i++){
    pos[j][i] = new Vec3(0,0,0);
    pos[j][i].x = stringTop.x+restLen*j;
    pos[j][i].y = stringTop.y + restLen*i; //Make each node a little lower
    pos[j][i].z = stringTop.z+20*i;
    vel[j][i] = new Vec3(0,0,0);
  } }
}

void update(float dt){

  //Reset accelerations each timestep (momenum only applies to velocity)
  for (int j = 0; j < numThreads; j++){
  for (int i = 0; i < numNodes; i++){
    acc[j][i] = new Vec3(0,0,0);
    acc[j][i].add(gravity);
  }}
  
  for (int j = 0; j < numThreads; j++){
  
  //Compute (damped) Hooke's law for each spring
  for (int i = 0; i < numNodes-1; i++){
    Vec3 diff = pos[j][i+1].minus(pos[j][i]);
    
    float stringF = -k*(diff.length() - restLen);
    //Vec3 frictionForce = vel[j][i].times(-kf);
    //println(stringF,diff.length(),restLen);
    
    Vec3 stringDir = diff.normalized();
    float projVbot = dot(vel[j][i], stringDir);
    float projVtop = dot(vel[j][i+1], stringDir);
    float dampF = -kv*(projVtop - projVbot);
    
    Vec3 force = stringDir.times(stringF+dampF);
    acc[j][i].add(force.times(-1.0/mass));
    //acc[j][i].add(frictionForce);
    acc[j][i+1].add(force.times(1.0/mass));
    
  }
  if(j<numThreads-1){
  for (int i = 0; i < numNodes; i++){
    Vec3 diff = pos[j+1][i].minus(pos[j][i]);
    
    float stringF = -k*(diff.length() - restLen);
    //Vec3 frictionForce = vel[j][i].times(-kf);
    //println(stringF,diff.length(),restLen);
    
    Vec3 stringDir = diff.normalized();
    float projVbot = dot(vel[j][i], stringDir);
    float projVtop = dot(vel[j+1][i], stringDir);
    float dampF = -kv*(projVtop - projVbot);
    
    Vec3 force = stringDir.times(stringF+dampF);
    acc[j][i].add(force.times(-1.0/mass));
    //acc[j][i].add(frictionForce);
    //println(i,j,acc[j+1][i]);
    acc[j+1][i].add(force.times(1.0/mass));
    
  }}

  //Eulerian integration
  for (int i = 1; i < numNodes; i++){
    vel[j][i].add(acc[j][i].times(dt));
    pos[j][i].add(vel[j][i].times(dt));
  }
  
  //Collision detection and response
  for (int i = 0; i < numNodes; i++){
    if (pos[j][i].y+radius > floor){
      vel[j][i].y *= -.9;
      pos[j][i].y = floor - radius;
    }
  }
  

    for (int i = 0; i < numNodes; i++){
      if (pos[j][i].distanceTo(Obstacle) < (obstacleR+radius)){
      Vec3 normal = (pos[j][i].minus(Obstacle)).normalized();
      pos[j][i] = Obstacle.plus(normal.times(obstacleR+radius).times(1.01));
      Vec3 velNormal = normal.times(dot(vel[j][i],normal));
      vel[j][i].subtract(velNormal.times(1 + 0.7));
    }
  } 
  
 
 
  
}

 //camera(mouseX, mouseY, 220.0, // eyeX, eyeY, eyeZ
 //        0.0, 0.0, 0.0, // centerX, centerY, centerZ
 //        0.0, 1.0, 0.0); // upX, upY, upZ
}

//Draw the scene: one sphere per mass, one line connecting each pair
boolean paused = true;
void draw() {
  background(255,255,255);
  noLights();

  camera.Update(1.0/frameRate);
  if (!paused) {for(int i=1;i<110 ;i++)update(1/(120*frameRate));}
  //if (!paused) {update(1/frameRate);}
  fill(0,0,0);
  for (int j = 0; j < numThreads; j++){
  for (int i = 0; i < numNodes-1; i++){
    pushMatrix();
    if(i<numNodes-1) {line(pos[j][i].x,pos[j][i].y, pos[j][i].z, pos[j][i+1].x,pos[j][i+1].y, pos[j][i+1].z);}
    //if(j<numThreads-1) {line(pos[j][i].x,pos[j][i].y,pos[j][i+1].x,pos[j][i+1].y);}
    translate(pos[j][i+1].x,pos[j][i+1].y,pos[j][i+1].z);
    sphere(radius);
    popMatrix();
  } }
  fill(50,100,70);
  pushMatrix();
    translate(Obstacle.x,Obstacle.y,Obstacle.z);
    sphere(obstacleR);
    popMatrix();
  
  if (paused)
    surface.setTitle(windowTitle + " [PAUSED]");
  else
    surface.setTitle(windowTitle + " "+ nf(frameRate,0,2) + "FPS");
}

void keyPressed(){
  if (key == ' ')
    paused = !paused;
   camera.HandleKeyPressed();
   
   if (key == 'r') {initScene();}
}


///////////////////
// Vec3D Library
///////////////////

public class Vec3 {
  public float x, y, z;
  
  public Vec3(float x, float y, float z){
    this.x = x;
    this.y = y;
    this.z = z;
  }
  
  public String toString(){
    return "(" + x+ ", " + y +", " + z +")";
  }
  
  public float length(){
    return sqrt(x*x+y*y+z*z);
  }
  
  public float lengthSqr(){
    return x*x+y*y+z*z;
  }
  
  public Vec3 plus(Vec3 rhs){
    return new Vec3(x+rhs.x, y+rhs.y, z+rhs.z);
  }
  
  public void add(Vec3 rhs){
    x += rhs.x;
    y += rhs.y;
    z += rhs.z;
  }
  
  public Vec3 minus(Vec3 rhs){
    return new Vec3(x-rhs.x, y-rhs.y, z-rhs.z);
  }
  
  public void subtract(Vec3 rhs){
    x -= rhs.x;
    y -= rhs.y;
    z -= rhs.z;
  }
  
  public Vec3 times(float rhs){
    return new Vec3(x*rhs, y*rhs, z*rhs);
  }
  
  public void mul(float rhs){
    x *= rhs;
    y *= rhs;
    z *= rhs;
  }
  
  public void normalize(){
    float magnitude = sqrt(x*x + y*y + z*z);
    x /= magnitude;
    y /= magnitude;
    z /= magnitude;
  }
  
  public Vec3 normalized(){
    float magnitude = sqrt(x*x + y*y + z*z);
    return new Vec3(x/magnitude, y/magnitude, z/magnitude);
  }
  
  public void clampToLength(float maxL){
    float magnitude = sqrt(x*x + y*y + z*z);
    if (magnitude > maxL){
      x *= maxL/magnitude;
      y *= maxL/magnitude;
      z *= maxL/magnitude;
    }
  }
  
  public void setToLength(float newL){
    float magnitude = sqrt(x*x + y*y + z*z);
    x *= newL/magnitude;
    y *= newL/magnitude;
    z *= newL/magnitude;
  }
  
  public float distanceTo(Vec3 rhs){
    float dx = rhs.x - x;
    float dy = rhs.y - y;
    float dz = rhs.z - z;
    return sqrt(dx*dx + dy*dy + dz*dz);
  }
  
}

Vec3 interpolate(Vec3 a, Vec3 b, float t){
  return a.plus((b.minus(a)).times(t));
}

float interpolate(float a, float b, float t){
  return a + ((b-a)*t);
}

float dot(Vec3 a, Vec3 b){
  return a.x*b.x + a.y*b.y+ a.z*b.z;
}

Vec3 projAB(Vec3 a, Vec3 b){
  return b.times(a.x*b.x + a.y*b.y +  a.z*b.z);
}
