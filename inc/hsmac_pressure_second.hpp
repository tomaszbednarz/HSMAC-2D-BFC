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


void pressure_velocity_correction()
{
   double xksi,yksi;
   double xeta,yeta;

   int i,j;

   inner = 0; 

   do {

      dmax = 1e-80;

      for (i=1; i<=NX; i++)
      {
         for (j=1; j<=NY; j++) 
         {

            double dp = ap11(i,j) * (1.0f/dksi)        * (un(i,j)-un(i-1,j)) +
                        ap12(i,j) * (1.0f/(4.0f*deta)) * (un(i,j+1)+un(i-1,j+1)-un(i,j-1)-un(i-1,j-1)) +
                        ap21(i,j) * (1.0f/(4.0f*dksi)) * (vn(i+1,j)+vn(i+1,j-1)-vn(i-1,j)-vn(i-1,j-1)) +
                        ap22(i,j) * (1.0f/deta)        * (vn(i,j)-vn(i,j-1));

            double delp = - pfactor * dp / dppp(i,j);

            if ( abs(dp) > abs(dmax) ) dmax = dp;

         // pressure correction
            p(i,j) += delp;

         // velocity correction
            un(i  ,j  ) += AA * dt * ae11(i  ,j  ) * (1/dksi) * delp;
            un(i-1,j  ) -= AA * dt * ae11(i-1,j  ) * (1/dksi) * delp;
            vn(i  ,j  ) += AA * dt * an22(i  ,j  ) * (1/deta) * delp;
            vn(i  ,j-1) -= AA * dt * an22(i  ,j-1) * (1/deta) * delp;
         }
      }
      inner ++;

      boundary_velocity();
      boundary_pressure();

   } while ( (abs(dmax) > 1.0e-3) && (inner < iiter) );

}


void init_pressure_velocity_correction()
{
   double xksi,yksi;
   double xeta,yeta;

   int i,j;

   cout << "*** calculating pressure velocity correction constants (dppp, jacobians, etc.)" << endl;
   cout << "*** version.second" << endl;

   inner = 0; 

   for (i=0; i<=NX+1; i++)
   {
      for (j=0; j<=NY+1; j++) 
      {
         // point P

         xksi = ( ((x(i,j)+x(i,j-1))/2.0f) - ((x(i-1,j)+x(i-1,j-1))/2.0f) ) / dksi;
         yksi = ( ((y(i,j)+y(i,j-1))/2.0f) - ((y(i-1,j)+y(i-1,j-1))/2.0f) ) / dksi;
         xeta = ( ((x(i,j)+x(i-1,j))/2.0f) - ((x(i,j-1)+x(i-1,j-1))/2.0f) ) / deta;
         yeta = ( ((y(i,j)+y(i-1,j))/2.0f) - ((y(i,j-1)+y(i-1,j-1))/2.0f) ) / deta;

         jacp(i,j) = xksi*yeta - xeta*yksi;

         ap11(i,j) =  yeta/jacp(i,j);
         ap12(i,j) = -yksi/jacp(i,j);
         ap21(i,j) = -xeta/jacp(i,j);
         ap22(i,j) =  xksi/jacp(i,j);  
      }
   }

   for (i=0; i<=NX+1; i++)
   {
      for (j=0; j<=NY+1; j++) 
      {
         // point e 

         xksi = ( ((x(i+1,j)+x(i+1,j-1))/2.0f) - ((x(i-1,j)+x(i-1,j-1))/2.0f) ) / (2.0f*dksi);
         yksi = ( ((y(i+1,j)+y(i+1,j-1))/2.0f) - ((y(i-1,j)+y(i-1,j-1))/2.0f) ) / (2.0f*dksi);
         xeta = ( x(i,j)-x(i,j-1) ) / deta;
         yeta = ( y(i,j)-y(i,j-1) ) / deta;

         jaceu(i,j) = xksi*yeta - xeta*yksi;

         ae11(i,j) =  yeta/jaceu(i,j);
         ae12(i,j) = -yksi/jaceu(i,j);
         ae21(i,j) = -xeta/jaceu(i,j);
         ae22(i,j) =  xksi/jaceu(i,j);  

         // point n 

         xksi = ( x(i,j)-x(i-1,j) ) / dksi;
         yksi = ( y(i,j)-y(i-1,j) ) / dksi;
         xeta = ( ((x(i,j+1)+x(i-1,j+1))/2.0f) - ((x(i,j-1)+x(i-1,j-1))/2.0f) ) / (2.0f*deta);
         yeta = ( ((y(i,j+1)+y(i-1,j+1))/2.0f) - ((y(i,j-1)+y(i-1,j-1))/2.0f) ) / (2.0f*deta);

         jacnv(i,j) = xksi*yeta - xeta*yksi;

         an11(i,j) =  yeta/jacnv(i,j);
         an12(i,j) = -yksi/jacnv(i,j);
         an21(i,j) = -xeta/jacnv(i,j);
         an22(i,j) =  xksi/jacnv(i,j);  

      }
   }

   for (i=1; i<=NX; i++)
   {
      for (j=1; j<=NY; j++) 
      {
         // for pressure correction         

         double dt1 = ap11(i,j) * (1.0f/(dksi*dksi))       * dt * (ae11(i,j)+ae11(i-1,j)); 
         double dt2 = ap12(i,j) * (1.0f/(16.0f*deta*deta)) * dt * (ae12(i,j+1)+ae12(i,j-1)+ae12(i-1,j+1)+ae12(i-1,j-1)); 
         double dt3 = ap21(i,j) * (1.0f/(16.0f*dksi*dksi)) * dt * (an21(i+1,j)+an21(i+1,j-1)+an21(i-1,j)+an21(i-1,j-1)); 
         double dt4 = ap22(i,j) * (1.0f/(deta*deta))       * dt * (an22(i,j)+an22(i,j-1)); 

         dppp(i,j) = dt1 + dt4 + dt2 + dt3;

         // for second derivatives
         double da11dksi = (ae11(i,j)-ae11(i-1,j))/(dksi);
         double da11deta = (an11(i,j)-an11(i,j-1))/(deta);
         double da21dksi = (ae21(i,j)-ae21(i-1,j))/(dksi);
         double da21deta = (an21(i,j)-an21(i,j-1))/(deta);
         
         ap2d1(i,j) = ap11(i,j)*da11dksi + ap12(i,j)*da11deta + ap21(i,j)*da21dksi + ap22(i,j)*da21deta;

         double da12dksi = (ae12(i,j)-ae12(i-1,j))/(dksi);
         double da12deta = (an12(i,j)-an12(i,j-1))/(deta);
         double da22dksi = (ae22(i,j)-ae22(i-1,j))/(dksi);
         double da22deta = (an22(i,j)-an22(i,j-1))/(deta);
         
         ap2d2(i,j) = ap11(i,j)*da12dksi + ap12(i,j)*da12deta + ap21(i,j)*da22dksi + ap22(i,j)*da22deta;
      }
   }

   // for 2nd derivaties at e
   for (i=1; i<NX; i++)
   {
      for (j=1; j<=NY; j++)
      {
         double da11dksi = (ap11(i+1,j)-ap11(i,j))/(dksi);
         double da11deta = (((an11(i+1,j)+an11(i,j))/2.0f)-((an11(i+1,j-1)+an11(i,j-1))/2.0f))/deta;
         double da21dksi = (ap21(i+1,j)-ap21(i,j))/(dksi);
         double da21deta = (((an21(i+1,j)+an21(i,j))/2.0f)-((an21(i+1,j-1)+an21(i,j-1))/2.0f))/deta;

         ae2d1(i,j) = ap11(i,j)*da11dksi + ap12(i,j)*da11deta + ap21(i,j)*da21dksi + ap22(i,j)*da21deta;

         double da12dksi = (ap12(i+1,j)-ap12(i,j))/(dksi);
         double da12deta = (((an12(i+1,j)+an12(i,j))/2.0f)-((an12(i+1,j-1)+an12(i,j-1))/2.0f))/deta;
         double da22dksi = (ap22(i+1,j)-ap22(i,j))/(dksi);
         double da22deta = (((an22(i+1,j)+an22(i,j))/2.0f)-((an22(i+1,j-1)+an22(i,j-1))/2.0f))/deta;

         ae2d2(i,j) = ap11(i,j)*da12dksi + ap12(i,j)*da12deta + ap21(i,j)*da22dksi + ap22(i,j)*da22deta;
      }
   }

   for (i=1; i<=NX; i++)
   {
      for (j=1; j<NY; j++)
      {
         double da11dksi = (((ae11(i,j+1)+ae11(i,j))/2.0f)-((ae11(i-1,j+1)+ae11(i-1,j))/2.0f))/dksi;
         double da11deta = (ap11(i,j+1)-ap11(i,j))/(deta);
         double da21dksi = (((ae21(i,j+1)+ae21(i,j))/2.0f)-((ae21(i-1,j+1)+ae21(i-1,j))/2.0f))/dksi;
         double da21deta = (ap21(i,j+1)-ap21(i,j))/(deta);
         
         an2d1(i,j) = ap11(i,j)*da11dksi + ap12(i,j)*da11deta + ap21(i,j)*da21dksi + ap22(i,j)*da21deta;

         double da12dksi = (((ae12(i,j+1)+ae12(i,j))/2.0f)-((ae12(i-1,j+1)+ae12(i-1,j))/2.0f))/dksi;
         double da12deta = (ap12(i,j+1)-ap12(i,j))/(deta);
         double da22dksi = (((ae22(i,j+1)+ae22(i,j))/2.0f)-((ae22(i-1,j+1)+ae22(i-1,j))/2.0f))/dksi;
         double da22deta = (ap22(i,j+1)-ap22(i,j))/(deta);
         
         an2d2(i,j) = ap11(i,j)*da12dksi + ap12(i,j)*da12deta + ap21(i,j)*da22dksi + ap22(i,j)*da22deta;
      }
   }
}


