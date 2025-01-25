//
// GRID GENERATOR 
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
// This code is solving Poisson eq using Thompson method
// for grid generation.
// Currently am working on a new code for this purpose for
// orthogonal meshes... Please ask me for new complately 
// changed version at the end of January.
//
//////////////////////////////////////////////////////////////////

// v0.2 :  attraction to line
// v0.1 :  computing elliptic grids without control points

#include <blitz/array.h>   

using namespace blitz;


Array<double,2> x;
Array<double,2> y;

Array<double,1> p;
Array<double,1> q;

Array<double,1> ai;   // ai coeff for eta-lines attraction
Array<double,1> ci;   // ci coeff for eta-lines attraction
Array<double,1> ksii; // ksii - eta to which line is attracted

Array<double,1> aj;   // aj coeff for eta-lines attraction
Array<double,1> cj;   // cj coeff for eta-lines attraction
Array<double,1> etaj; // etaj - eta to which line is attracted

unsigned int NXCL;
unsigned int NYCL;

double rel;

double xksi;
double yksi;
double xeta;
double yeta;
double alp;
double bet;
double gam;                                                             
double xksieta;
double yksieta;

double dksi;
double deta;

double errx;
double erry;
double errmax = 1.0f;

unsigned int NX;
unsigned int NY;

unsigned int iter;

double error;

// my research
void boundary_conditions()
{
   firstIndex  i;
   secondIndex j;

// south boundary
   x(Range(0,NX-1), Range(0)) = 0.2f+(2.0f-0.2f)*(i/(double)(NX-1));
   y(Range(0,(NX-1)/2), Range(0)) = -0.02-0.08*(i/(double)((NX-1)/2));
   y(Range((NX-1)/2+1,NX-1), Range(0)) = -0.1;
// north boundary
   x(Range(0,NX-1), Range(NY-1)) = 0.2f+(2.0f-0.2f)*(i/(double)(NX-1));
   y(Range(0,NX-1), Range(NY-1)) = 0.0f;
// west boundary
   x(Range(0), Range(0,NY-1)) = 0.2f;
   y(Range(0), Range(0,NY-1)) = - (-0.02*(j/(double)(NY-1)) + 0.02);  
// east boundary
   x(Range(NX-1), Range(0,NY-1)) = 2.0f;
   y(Range(NX-1), Range(0,NY-1)) = - (-0.1*(j/(double)(NY-1)) + 0.1);  
/*
// controls
   NYCL = 1;
   aj.resize(shape(NYCL));
   cj.resize(shape(NYCL));
   etaj.resize(shape(NYCL));  
   aj(0) = -10.0f;
   cj(0) = 0.01f;
   etaj(0) = 0.0f;
*/
   NXCL = 0;
/*
   ai.resize(shape(NXCL));
   ci.resize(shape(NXCL));
   ksii.resize(shape(NXCL));  
   ai(0) = -15.0f;
   ci(0) = 0.01f;
   ksii(0) = 1.0f;
*/
}


// cube ?
void boundary_conditions2()
{
   firstIndex  i;
   secondIndex j;

// south boundary
   x(Range(0,NX-1), Range(0)) = 0.2f+(2.0f-0.2f)*(i/(double)(NX-1));
   y(Range(0,NX-1), Range(0)) = -0.02-0.08*(i/(double)(NX-1));

//   y(Range(0,N-1), 0) = -0.1;
// north boundary
   x(Range(0,NX-1), Range(NY-1)) = 0.2f+(2.0f-0.2f)*(i/(double)(NX-1));
   y(Range(0,NX-1), Range(NY-1)) = 0.0f;
//   y(Range::all(), M-1) = 1;
// west boundary
   x(Range(0), Range(0,NY-1)) = 0.2f;
   y(Range(0), Range(0,NY-1)) = - (-0.02*(j/(double)(NY-1)) + 0.02);  

// east boundary
   x(Range(NX-1), Range(0,NY-1)) = 2.0f;
   y(Range(NX-1), Range(0,NY-1)) = - (-0.1*(j/(double)(NY-1)) + 0.1);  
}


// zaokraglone ?
void boundary_conditions3()
{
   firstIndex  i;
   secondIndex j;

// south boundary
   x(Range(0,NX-1), Range(0)) = (2.0f)*(i/(double)(NX-1));
   y(Range(0,NX-1), Range(0)) = -1.0f*(sin(3.14*i/(double)(NX-1)))-0;

// north boundary
   x(Range(0,NX-1), Range(NY-1)) = (2.0f)*(i/(double)(NX-1));
   y(Range(0,NX-1), Range(NY-1)) = 2.0f;

// west boundary
   x(Range(0), Range(0,NY-1)) = 0.0f;
   y(Range(0), Range(0,NY-1)) = (2.0*(j/(double)(NY-1)));  

// east boundary
   x(Range(NX-1), Range(0,NY-1)) = 2.0f;
   y(Range(NX-1), Range(0,NY-1)) = (2.0*(j/(double)(NY-1)));  


}

// ============================================================================ 
//  solve poisson equation with SSOR method
// ============================================================================
void solve_xy()
{ 
   double ap, ae, aw, an, as, left, cj, xold, yold;

   iter = 0;
   while (errmax > error)
   {
      errmax = 1.0;
      for (int j=1; j<(NY-1); j++)
      {
         for (int i=1; i<(NX-1); i++)
         {
            xksi = ( x(i+1,j) - x(i-1,j) ) / (2*dksi);
            yksi = ( y(i+1,j) - y(i-1,j) ) / (2*dksi);
            xeta = ( x(i,j+1) - x(i,j-1) ) / (2*deta);
            yeta = ( y(i,j+1) - y(i,j-1) ) / (2*deta);
         
            alp = xeta*xeta + yeta*yeta;
            bet = xksi*xeta + yksi*yeta;
            gam = xksi*xksi + yksi*yksi;

            xksieta = ( x(i+1,j+1) - x(i-1,j+1) - x(i+1,j-1) + x(i-1,j-1) ) / ( 4*dksi*deta );
            yksieta = ( y(i+1,j+1) - y(i-1,j+1) - y(i+1,j-1) + y(i-1,j-1) ) / ( 4*dksi*deta );

            cj = xksi*yeta - xeta*yksi;

            ap = 2*( (alp/(dksi*dksi)) + (gam/(deta*deta)) ); 
            ae = alp / (dksi*dksi) + cj*cj*(p(i) / (2*dksi));
            aw = alp / (dksi*dksi) - cj*cj*(p(i) / (2*dksi));
            an = gam / (deta*deta) + cj*cj*(q(j) / (2*deta));
            as = gam / (deta*deta) - cj*cj*(q(j) / (2*deta)); 

            left = ae*x(i+1,j) + aw*x(i-1,j) + an*x(i,j+1) + as*x(i,j-1) - 2*bet*xksieta;
            xold = x(i,j);
            x(i,j) = rel*( left )/ap + (1-rel)*x(i,j);
            errx = abs( xold - x(i,j) );

            left = ae*y(i+1,j) + aw*y(i-1,j) + an*y(i,j+1) + as*y(i,j-1) - 2*bet*yksieta;
            yold = y(i,j);
            y(i,j) = rel*( left )/ap + (1-rel)*y(i,j);
            erry = abs( yold - y(i,j) );

            errmax = max(errx,erry);         
          }
      }
      iter++;
      printf ("iter: %-5d, errmax: %-5.6e\n",iter,errmax);
   }
}

double attraction_line(double a, double c, double xi)
{
   return (double)( a*xi/abs(xi)*exp(-c*abs(xi)) );
}

double attraction_point(double b, double d, double ksid, double etad)
{
   return (double)( b*ksid/abs(ksid)*exp(-d*sqrt(ksid*ksid+etad*etad)) );
}

void calculate_controls()
{
   double difeta, difksi;
   int i,i1;
   int j,j1;

   for (i=0; i<NX; i++)
   {
      p(i) = 0.0f;

      if (NXCL > 0)
      for (i1=0; i1<NXCL; i1++)
      {
         difksi = ksii(i1)-dksi*i;
         p(i) = p(i) - attraction_line(ai(i1),ci(i1),difksi);       
      }
   }

   for (j=0; j<NY; j++)
   {
      q(j) = 0.0f;

      if (NYCL > 0)
      for (j1=0; j1<NYCL; j1++)
      {
         difeta = etaj(j1)-deta*j;
         q(j) = q(j) - attraction_line(aj(j1),cj(j1),difeta);       
      }
   }

/*
   int MM = 1,i,j,i1,j1;  
   double dd,bb;

   for (i=0; i<NX; i++)
   {
      for (i1=0; i1<MM; i1++)
      {
//         dd = 0 - (double)i;
//         if (i != 0) p(i) = p(i) - attraction_line(15, 0.1, dd);


         dd = 0.5 - i/(double)(NX-1);
//         if (i != 0) 
//         p(i) = p(i) - attraction_line(1.0, 0.1, dd);
      }   
   }

   for (j=0; j<NY; j++)
   {
      for (j1=0; j1<MM; j1++)
      {

         dd = 0.25 - j/(double)(NX-1);
         bb = 0.25 - j/(double)(NY-1);
//         if (i != 0) 
         p(j) = p(j) - attraction_point(-0.25, 0.01, dd,bb);
      }   
   }

//p(15)=30;

   for (j=0; j<NY; j++)
   {
      for (j1=0; j1<MM; j1++)
      {
//         dd = 0 - (double)j;
//         if (j != 0) q(j) = q(j) - ss(-20, 0.01, dd);
      }   
   }


   for (j=0; j<NY; j++)
   {
      for (j1=0; j1<MM; j1++)
      {

         dd = 0.25 - j/(double)(NX-1);
         bb = 0.25 - j/(double)(NY-1);
//         if (i != 0) 
         q(j) = q(j) - attraction_point(-0.25, 0.01, dd,bb);
      }   
   }


//q(1)=-50;

*/
}


void write_dat()
{
   FILE *fp;
   fp=fopen("grid.dat","wt");

   fprintf(fp,"VARIABLES = x, y\n");
   fprintf(fp,"ZONE I=%d, J=%d\n",NX,NY);

   for (int j=0; j<NY; j++) {
      for (int i=0; i<NX; i++) {
         fprintf (fp,"%f  %f\n",x(i,j),y(i,j));
      }
   }
   fclose(fp);
}


int main()
{
   // this to be read from file
   NX = 40;
   NY = 40;
   dksi = 2.0f / (double)(NX-1);
   deta = 0.1f / (double)(NY-1);  
   rel = 1.6;
   error = 1e-16;
   iter = 0;

   allocateArrays(shape(NX,NY),x,y);
   p.resize(shape(NX));
   q.resize(shape(NY));
   p = 0;
   q = 0;

   boundary_conditions3();
   calculate_controls();

   solve_xy();

   write_dat();

   return 0;
}

