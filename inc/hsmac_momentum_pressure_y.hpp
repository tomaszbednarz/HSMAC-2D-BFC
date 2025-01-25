

         peta = ( p(i,j+1)-p(i,j) ) / deta;
         pksi = ( (p(i+1,j+1)+p(i+1,j))-(p(i-1,j+1)+p(i-1,j)) ) / (4.0f*dksi);

         double PRESSV = -AA * (an21(i,j)*pksi + an22(i,j)*peta);
