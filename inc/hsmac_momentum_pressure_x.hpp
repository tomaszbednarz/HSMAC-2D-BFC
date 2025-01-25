

         pksi = ( p(i+1,j)-p(i,j) ) / dksi;
         peta = ( (p(i+1,j+1)+p(i,j+1))-(p(i+1,j-1)+p(i,j-1)) ) / (4*deta);

         double PRESSU = -AA * (ae11(i,j)*pksi + ae12(i,j)*peta);
