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

#define momentum_adv 3

void calculate_velocity()
{
   int i,j;

   double xksi, yksi, xeta, yeta;
   double pksi, peta;
   double alpha, beta, gamma, sigma, tau;
   double uksiksi, uksieta, uetaeta;
   double xksiksi, xksieta, xetaeta;
   double yksiksi, yksieta, yetaeta;
   double vksiksi, vksieta, vetaeta;

   // =======================================
   // u component
   // =======================================
   for (i=1; i<NX; i++)
   {
      for (j=1; j<=NY; j++)
      {

#include "hsmac_momentum_pressure_x.hpp"

         double uksi = ( uo(i+1,j)-uo(i-1,j) ) / (2.0f*dksi);
         double ueta = ( uo(i,j+1)-uo(i,j-1) ) / (2.0f*deta);

#if (momentum_adv == 2)
	#include "hsmac_momentum_utopia_x.hpp"
#endif

#if (momentum_adv == 3)
	#include "hsmac_momentum_utopia_x_3.hpp"
#endif

         uksiksi = ( ( (uo(i+1, j )-uo( i , j ))/dksi ) 
                   - ( (uo( i , j )-uo(i-1, j ))/dksi ) ) / dksi; 

         uetaeta = ( ( (uo( i ,j+1)-uo( i , j ))/(deta) ) 
                   - ( (uo( i , j )-uo( i ,j-1))/(deta) ) ) / deta;                   

         uksieta = ( ( (uo(i+1,j+1)-uo(i-1,j+1))/(2.0f*dksi) ) 
                   - ( (uo(i+1,j-1)-uo(i-1,j-1))/(2.0f*dksi) ) ) / (2.0f*deta); 


         double DIFFU = BB * (
                             ( SQR(ae11(i,j))+SQR(ae21(i,j)) )*uksiksi + 
                             ( SQR(ae12(i,j))+SQR(ae22(i,j)) )*uetaeta +
                             ( (ae11(i,j)*ae12(i,j)) + (ae21(i,j)*ae22(i,j)) )*uksieta*2.0f +
                               uksi*ae2d1(i,j) + ueta*ae2d2(i,j) );


         un(i,j) = uo(i,j) + dt * ( PRESSU + DIFFU - FUX - FUY );
      }
   }

   // =======================================
   // v component
   // =======================================
   for (i=1; i<=NX; i++)
   {
      for (j=1; j<NY; j++)
      {
	#include "hsmac_momentum_pressure_y.hpp"

         double vksi = ( vo(i+1,j)-vo(i-1,j) ) / (2.0f*dksi);
         double veta = ( vo(i,j+1)-vo(i,j-1) ) / (2.0f*deta);

#if (momentum_adv == 2)
	#include "hsmac_momentum_utopia_y.hpp"
#endif

#if (momentum_adv == 3)
	#include "hsmac_momentum_utopia_y_3.hpp"
#endif

         vksiksi = ( ( (vo(i+1, j )-vo( i , j ))/dksi ) 
                   - ( (vo( i , j )-vo(i-1, j ))/dksi ) ) / dksi; 

         vetaeta = ( ( (vo( i ,j+1)-vo( i , j ))/(deta) ) 
                   - ( (vo( i , j )-vo( i ,j-1))/(deta) ) ) / deta;                   

         vksieta = ( ( (vo(i+1,j+1)-vo(i+1,j-1))/(2*deta) ) 
                   - ( (vo(i-1,j+1)-vo(i-1,j-1))/(2*deta) ) ) / (2*dksi); 


         double DIFFV = BB * (
                             ( SQR(an11(i,j))+SQR(an21(i,j)) )*vksiksi + 
                             ( SQR(an12(i,j))+SQR(an22(i,j)) )*vetaeta +
                             ( (an11(i,j)*an12(i,j)) + (an21(i,j)*an22(i,j)) )*vksieta*2.0f +
                               vksi*an2d1(i,j) + veta*an2d2(i,j) );
  
       
         double CONV = CC * ( (to(i,j)+(to(i,j+1)))*0.5 - T_0);

         vn(i,j) = vo(i,j) + dt * ( PRESSV + DIFFV - FVX - FVY + CONV );
          
      }
   }
   
   boundary_velocity();
}

