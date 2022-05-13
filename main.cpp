#include <iostream>
#include <time.h>
#include "CImg.h"
#include "DEIntegrator.h" 
#include <math.h>

#define pi 3.1415926535897932384

#define min(a,b) (!(b<a)?a:b)
#define max(a,b) ((b<a)?a:b)
#define clamp(x,low,high) min(max(x, low), high)
#define dot(x,y) (x[0]*y[0] + x[1]*y[1] + x[2]*y[2])
#define cross(x,y) {x[1]*y[2] - x[2]*y[1], x[2]*y[0] - x[0]*y[2], x[0]*y[1] - x[1]*y[0]}
//#define cross(z,x,y) z[0]=x[1]*y[2] - x[2]*y[1]; z[1]=x[2]*y[0] - x[0]*y[2]; z[2]= x[0]*y[1] - x[1]*y[0]

//Remez Rational function for the integral of exp(-sqrt(1+x^2)) from 0 to 5 
//(0.00006052784551827664 + (0.09635940796022202215 + 0.030794431055283765117*x)*x) numerator
//(0.26505636418297421598 + (0.07373876951945245253 + 0.058730315948628762954*x)*x) denominator

//Remez Rational function for the integral of exp(-sqrt(1+x^2)), then divided by exp(x) from 0 to 5 
//0.000310541706096794180 + (0.024107390518790483710 - 0.0047822317672900930816*x)*x
//0.08848134212047325803 + (-0.02551218768843382409 + 0.10403190689073187223*x)*x

using namespace cimg_library;

const double r_e = 6371000; //Earth's radius in metres
const double atmos = 30000; //Additional radius of atmosphere, constant density being assumed.


// Applied Rodrigues' formula for rotation around an axis.
//inline Vector3d rotate(Vector3d axis,double angle,Vector3d vec){return cos(angle)*vec +sin(angle)*axis.cross(vec) + (1-cos(angle))*axis.dot(vec)*axis;}

const int width = 500;
const int height = 500;
float xFOV = 0.5; //Field of View in x and y directions
float yFOV = xFOV*height/width;
CImg<unsigned char> image = CImg<int8_t>(width,height,1,3,0);
CImgDisplay disp(image,"Sunset");
bool update = 1;

double thetaS=0;
double phiS=0;
double phiV = 0;
double thetaV = 0;

double temp[3] = {};

double expint = 0.601907230197; //integral of exp(-sqrt(1+x^2)) from 0 to infinity.

const double blue = 435.8e-9; //Wavelengths of R, G and B light
const double green = 546.1e-9;
const double red = 700e-9;

const double blue4 = pow(blue,4); //Some constants used in Rayleigh scattering later
const double green4 = pow(green,4);
const double red4 = pow(red,4);

const double rb = 6.656458924474147; // red4/blue4, used in a pow later
const double rg = 2.699625078466617; // red4/green4, used in a pow later

inline double distance(double* r,double* d, double h){ //Distance that a line in direction d travels through the atmosphere from r.
  temp[0] = dot(r,d);
  return sqrt(h*h - dot(r,r) + temp[0]*temp[0]) - temp[0];
}

inline double distance1(double x, double r1, double r2, double r3, double sun1, double sun2, double sun3, double h){
  //Same as above but written out in terms of components
  return sqrt(h*h - (r1*r1*x*x + r2*r2*x*x + (r_e+r3*x)*(r_e+r3*x)) + (r1*x*sun1 + r2*x*sun2 + (r_e+r3*x)*sun3)*(r1*x*sun1 + r2*x*sun2 + (r_e+r3*x)*sun3)) - (r1*x*sun1 + r2*x*sun2 + (r_e+r3*x)*sun3);
}

double powrb(double x) //remez polynomial for pow(x,rb) for x on [0.1,0.9]
{
    double u = 4.9299156782812124;
    u = u * x + -6.6277602664640747;
    u = u * x + 3.2418485928370049;
    u = u * x + -6.4070680132140313e-1;
    return u * x + 4.078590289906773e-2;
}

double powrg(double x) //remez polynomial for pow(x,rg) for x on [0.1,0.9]
{
    double u = 6.7347141422882724e-1;
    u = u * x + 3.8538863754414188e-1;
    u = u * x + -6.019373471105763e-2;
    return u * x + 3.8490669167774629e-3;
}

inline double absorb(double d,double lambda4){ //Can later replace with a more sophisticated model
  const double sigma = green4*1e-5; // the constant of proportionality representing the probability of a photon being scattered per meter.
  // currently, a reduction of 10^-5 per metre for green light
  return exp(-d*sigma/lambda4);
}

class Function1 //This is the function to integrate, it must be univariate so the other arguments are passed in as properties
{
public:
	double r1; double r2; double r3;
	double sun1; double sun2; double sun3;
	double operator()(double x) const
	{
		return absorb(x+distance1(x,r1,r2,r3,sun1,sun2,sun3,r_e+atmos),red4);
	}
};
//integral of (absorption from remaining atmosphere for scattered wave to travel through)
//            *(scattering probability at that angle and wavelength)
//            *(absorption from atmsophere travelled through to get to that particle)
//Not actually used, but the value this gives is what the integral function should give
void integrate(double high,double precision, double* d, double* sun,double* integral){
  Function1 f1;
  f1.r1 = d[0]; f1.r1 = d[1]; f1.r1 = d[2];
  f1.sun1 = sun[0]; f1.sun2 = sun[1]; f1.sun3 = sun[2];
  double total = DEIntegrator<Function1>::Integrate(f1, 0, high, precision);
  double temp = (1+dot(d,sun)*dot(d,sun));
  integral[0] = temp*total/red4; integral[1] = temp*powrg(total)/green4; integral[2] = temp*powrb(total)/blue4;
//return total*(1+dot(d,sun)*dot(d,sun))/lambda4;
} 

// A very crude approximation to the above integral, it goes as far out as 2x the intended value in some cases but still looks decent
// It takes the value of the integrand at the centre of the domain, so the midpoint rule with n=1
void integrate2(double high, double* d, double* sun, double* integral){
  double r[3] = {high*d[0]/2,high*d[1]/2,r_e+high*d[2]/2};
  double total = absorb(high/2 + distance(r,sun,r_e+atmos),red4);
  double temp = high*(1+dot(d,sun)*dot(d,sun));
  integral[0] = temp*total/red4; integral[1] = temp*powrg(total)/green4; integral[2] = temp*powrb(total)/blue4;
}

int main(){
  double start[3] = {0,0,r_e};
  int b=0; //frame counter
  double integral[3] = {};
  double sunC[3]; //The sun's colour
  disp.hide_mouse();
  while (!disp.is_closed()){
  b+=1;
  if (disp.is_keyARROWUP()!=disp.is_keyARROWDOWN() or disp.is_keyARROWRIGHT()!=disp.is_keyARROWLEFT() or disp.mouse_x()!=width/2 or disp.mouse_y()!=height/2){update = 1;}
  phiS+=(disp.is_keyARROWUP()-disp.is_keyARROWDOWN())*0.02;
  thetaS+=(disp.is_keyARROWRIGHT()-disp.is_keyARROWLEFT())*0.02;
  thetaV = fmod(thetaV + (disp.mouse_x() - width/2)/100.0, 2*pi);
  phiV = clamp(phiV + (disp.mouse_y() - height/2)/100.0,0,pi);
  disp.set_mouse(width/2,height/2);
  if (update){
  
  disp.set_mouse(width/2,height/2);
  double sun[3] = {sin(phiS)*cos(thetaS),sin(phiS)*sin(thetaS),cos(phiS)};

  sunC[0]=absorb(distance(start,sun,r_e+atmos),red4);
  sunC[1]=absorb(distance(start,sun,r_e+atmos),green4);
  sunC[2]=absorb(distance(start,sun,r_e+atmos),blue4);

  double view[3] = {sin(phiV)*cos(thetaV),sin(phiV)*sin(thetaV),cos(phiV)}; //centre of view
  double yaxis[3] = {sin(thetaV),-cos(thetaV),0}; // axis of rotation for vertical movement
  double xaxis[3] = cross(view,yaxis); // axis of rotation for horizontal movement


    for (int x=0; x<width; x++){
      for (int y=0; y<height; y++){
	//For each pixel, cast a line in that direction, and sum the contribution from each particle along that line
        //So an integral of (absorption from remaining atmosphere for scattered wave to travel through)
        //                 *(scattering probability at that angle and wavelength)
        //                 *(absorption from atmsophere travelled through to get to that particle)
	temp[0] = cos(x*xFOV/width -xFOV/2) ; temp[1] = (-sin(x*xFOV/width-xFOV/2));
        double d[3] = {view[0]*temp[0] + yaxis[0]*temp[1],view[1]*temp[0] + yaxis[1]*temp[1],view[2]*temp[0] + yaxis[2]*temp[1]};

	temp[0] = cos(y*yFOV/height -yFOV/2) ; temp[1] = (-sin(y*yFOV/height-yFOV/2));
	d[0] = d[0]*temp[0] + xaxis[0]*temp[1]; d[1] = d[1]*temp[0] + xaxis[1]*temp[1]; d[2] = d[2]*temp[0] + xaxis[2]*temp[1];

	if (d[2]>0){
        //integrate(distance(start,d,r_e+atmos),1e20,d,sun,integral);
        integrate2(distance(start,d,r_e+atmos),d,sun,integral);

        image(x,height-y-1,0)=clamp(integral[0]/5e27,0,255);
        image(x,height-y-1,1)=clamp(integral[1]/5e27,0,255);
        image(x,height-y-1,2)=clamp(integral[2]/5e27,0,255);
	if(dot(d,sun)>0.9995){
          image(x,height-y-1,0)=min(255*sunC[0]+image(x,height-y-1,0),255);
          image(x,height-y-1,1)=min(255*sunC[1]+image(x,height-y-1,1),255);
          image(x,height-y-1,2)=min(255*sunC[2]+image(x,height-y-1,2),255);
          }
	}
	else{
          image(x,height-y-1,0)=30;
          image(x,height-y-1,1)=140;
          image(x,height-y-1,2)=0;
        }
      }
    }
  image.display(disp);
  std::cout << disp.frames_per_second() << "\n";
  update = 0;
  }
  cimg::sleep(1);
  }
return 0;
}