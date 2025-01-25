

         // UTOPIA 

         if ((i==1) || (i==NX)) {
            o1 = ( (to(i+1,j)-to(i-1,j)) - 1.0f*sgn(ut)*(to(i+1,j)-2*to(i,j)+to(i-1,j)) ) / (2*dksi);
         } else {
            if (ut>0.0f) dif=1.0f; else if (ut<0.0f) dif=-1.0f; else dif=0.0f;
            o1= (  (dif-1)/12 * to(i+2,j)
                  -(dif-2)/ 3 * to(i+1,j)
                  +(dif  )/ 2 * to(i  ,j)                     
                  -(dif+2)/ 3 * to(i-1,j)
                  +(dif+1)/12 * to(i-2,j) ) / dksi;
         }
         double FTX = ut / J_t * (yeta * o1 - yksi * teta);

         if ((j==1) || (j==NY)) {
            o1 = ( (to(i,j+1)-to(i,j-1)) - 1.0*sgn(vt)*(to(i,j+1)-2*to(i,j)+to(i,j-1)) ) / (2*deta);
         } else {
            if (vt>0.0f) dif=1.0f; else if (vt<0.0f) dif=-1.0f; else dif=0.0f;
            o1= (  (dif-1)/12 * to(i,j+2)
                  -(dif-2)/ 3 * to(i,j+1)
                  +(dif  )/ 2 * to(i,j  )                     
                  -(dif+2)/ 3 * to(i,j-1)
                  +(dif+1)/12 * to(i,j-2) ) / deta;
         }
         double FTY = vt / J_t * (xksi * o1 - xeta * tksi);
