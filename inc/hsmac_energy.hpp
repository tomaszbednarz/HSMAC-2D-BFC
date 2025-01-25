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

#define energy_adv 3
// 1 - central
// 2 - UTOPIA
// 3 - UTOPIA 2

void calculate_energy()
{
   int i,j;

   double xksi, yksi, xeta, yeta;
   double tksi, teta;
   double dif, o1;
   double alpha, beta, gamma, sigma, tau;
   double tksiksi, tetaeta, tksieta;
   double xksiksi, yksiksi, xetaeta, yetaeta, xksieta, yksieta;

   for (i=1; i<=NX; i++)
   {
      for (j=1; j<=NY; j++)
      {

         tksi = ( to(i+1,j)-to(i-1,j) ) / (2.0f*dksi);

         teta = ( to(i,j+1)-to(i,j-1) ) / (2.0f*deta);

         tksiksi = ( to(i+1,j)-2.0f*to(i,j)+to(i-1,j) ) / (dksi*dksi);

         tetaeta = ( to(i,j+1)-2.0f*to(i,j)+to(i,j-1) ) / (deta*deta);

         tksieta = ( (to(i+1,j+1)-(to(i+1,j-1)))/(2.0f*deta) 
                   - (to(i-1,j+1)-(to(i-1,j-1)))/(2.0f*deta) ) / (2.0f*dksi);

         double ut = ( uo(i,j)+uo(i-1,j) ) / 2.0f;
         double vt = ( vo(i,j)+vo(i,j-1) ) / 2.0f;

         // calculate FTX, FTY

#if (energy_adv == 1)
         #include "hsmac_energy_central.hpp"
#endif

#if (energy_adv == 2)
         #include "hsmac_energy_utopia.hpp"
#endif

#if (energy_adv == 3)
         #include "hsmac_energy_utopia2.hpp"
#endif

         double DIFT = DD * (
                            ( SQR(ap11(i,j))+SQR(ap21(i,j)) )*tksiksi + 
                            ( SQR(ap12(i,j))+SQR(ap22(i,j)) )*tetaeta +
                            ( (ap11(i,j)*ap12(i,j)) + (ap21(i,j)*ap22(i,j)) )*tksieta*2.0f 
                            + tksi*ap2d1(i,j) + teta*ap2d2(i,j) );

         tn(i,j) = to(i,j) + dt * ( - FTX - FTY + DIFT );

      }
   }

   boundary_temperature();
}

