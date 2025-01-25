//
// BFC + HSMAC method 
// 
// Programmed by Dr Tomasz BEDNARZ of JCU (School of Engineering)
// Email: tomasz.bednarz@jcu.edu.au
//
// November 2006, version 20061101.3
//
// All rights reserved (C)
//
// If you have this code, it means that it was given to you by
// the author.  If you use this code or any part of it for your
// own purposes, please mention the name of author and add his
// name to any papers or reports published with this ideas here.
// Also if you prefer F90, the code written in easy way to
// convert it to F90 - just copy and paste and change the 
// looping functions, etc...
// However, I recommend you to use C++ version due to many
// changes coming in the near future... 
//
// to compile the program you need to have:
// GNU C/C++
// Blitz++ library
//
// The code is solving NS equation on orthogonal meshes.
//
// REMEBER:
// The code is still under development process - please ask the
// me every 2nd month for the new and updated version.
//
// Things under developent:
// - corrections for non-orthogonality and highly skew grids
// - XML configuration files for universality
// - relative iteration progress notifications
// - particle transport 
// - postprocessing functions
//
//////////////////////////////////////////////////////////////////

#define VersionString "20061030" // ONLY 8 chars stored as a version number!

#define SQR(a) ((a)*(a))


Array<double,2> x;
Array<double,2> y;

Array<double,2> p;
Array<double,2> uo;
Array<double,2> vo;
Array<double,2> un;
Array<double,2> vn;
Array<double,2> tn;
Array<double,2> to;

Array<double,2> dppp;

Array<double,2> jacp;
Array<double,2> jaceu;
Array<double,2> jacnv;

Array<double,2> psi;

Array<double,2> ap11;
Array<double,2> ap12;
Array<double,2> ap21;
Array<double,2> ap22;

Array<double,2> ap2d1; // for 2nd derivatives 
Array<double,2> ap2d2; // for 2nd derivatives 
Array<double,2> an2d1; // for 2nd derivatives 
Array<double,2> an2d2; // for 2nd derivatives 
Array<double,2> ae2d1; // for 2nd derivatives 
Array<double,2> ae2d2; // for 2nd derivatives 

Array<double,2> ae11;
Array<double,2> ae12;
Array<double,2> ae21;
Array<double,2> ae22;

Array<double,2> an11;
Array<double,2> an12;
Array<double,2> an21;
Array<double,2> an22;


#define sgn(x) ((x)>0 ? 1 : ((x)==0 ? 0:(-1)))

int NX, NY;

#define LKSI 0.032f  //0.3f
#define LETA 0.032f  //0.015f

int iiter;
int niter;

int iter;

double pfactor = 1.7f; // pressure relaxation factor
double dksi; //,dxi;
double deta; //,det;
double dt;

double AA;
double BB;
double CC;
double DD;
double T_0;

void boundary_velocity();
void boundary_pressure();
void boundary_temperature();
void init_pressure_velocity_correction();

double dmax;
int inner;

void init_tables(int nx, int ny);
