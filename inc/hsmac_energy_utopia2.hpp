

         // UTOPIA, method 2 

         double unew = ap11(i,j)*ut + ap21(i,j)*vt;
         double vnew = ap12(i,j)*ut + ap22(i,j)*vt;

         if ((i==1) || (i==NX)) {
            o1 = ( (to(i+1,j)-to(i-1,j)) - 1.0f*sgn(unew)*(to(i+1,j)-2*to(i,j)+to(i-1,j)) ) / (2*dksi);
         } else {
            if (unew>0.0f) dif=1.0f; else if (unew<0.0f) dif=-1.0f; else dif=0.0f;
            o1= (  (dif-1)/12 * to(i+2,j)
                  -(dif-2)/ 3 * to(i+1,j)
                  +(dif  )/ 2 * to(i  ,j)                     
                  -(dif+2)/ 3 * to(i-1,j)
                  +(dif+1)/12 * to(i-2,j) ) / dksi;
         }
         double FTX = unew*o1; 

         if ((j==1) || (j==NY)) {
            o1 = ( (to(i,j+1)-to(i,j-1)) - 1.0*sgn(vnew)*(to(i,j+1)-2*to(i,j)+to(i,j-1)) ) / (2*deta);
         } else {
            if (vnew>0.0f) dif=1.0f; else if (vnew<0.0f) dif=-1.0f; else dif=0.0f;
            o1= (  (dif-1)/12 * to(i,j+2)
                  -(dif-2)/ 3 * to(i,j+1)
                  +(dif  )/ 2 * to(i,j  )                     
                  -(dif+2)/ 3 * to(i,j-1)
                  +(dif+1)/12 * to(i,j-2) ) / deta;
         }
         double FTY = vnew*o1;


