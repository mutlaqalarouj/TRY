[max_y, max_x, max_z]=size(P);
w = 1;
len = length(NNCCL);
C = sparse(1,ne);
D = sparse(1,ne);
HOSTNUM3D = reshape(HOSTNUM,DY,DX,DZ);

for ip = 1:DZ
        
        HOSTNUM1(:,:,ip) = HOSTNUM3D(:,:,ip)';
    end

for w = 1:len
    
             i = NNCCL(w,1);
             j = NNCCL(w,2);
             k = NNCCL(w,3);
   
    if (NNCG(w) - HOSTNUM1(i,j,k)) == 1 
        
%         i = NNCCL(w,1); j = NNCCL(w,2); k = NNCCL(w,3); 
        n=(k-1)*(max_y*max_x)+(j-1)*max_y+i;
        b_Wx = -L_ek(i,j,k)*L_ek(i,j-1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(i,j,k)*0.5*(x(j)-x(j-1))+L_ek(i,j-1,k)*0.5*(x(j+1)-x(j)));
        b_P =  b_Wx ;
        B(n,n) = -b_P + B(n,n);
        C(n) = b_P * Pn(NNCCG(w,1),NNCCG(w,2),NNCCG(w,3)) + C(n);
        a_Wx = -sigma(i,j,k)*sigma(i,j-1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(sigma(i,j,k)*0.5*(x(j)-x(j-1))+sigma(i,j-1,k)*0.5*(x(j+1)-x(j)));
        a_P =  a_Wx ;
        A(n,n) = -a_P + A(n,n);
        D(n) = a_P * SPn(NNCCG(w,1),NNCCG(w,2),NNCCG(w,3)) + D(n);
           
    elseif (NNCG(w) - HOSTNUM1(i,j,k)) == 139
        
%         i = NNCCL(w,1); j = NNCCL(w,2); k = NNCCL(w,3); 
        n=(k-1)*(max_y*max_x)+(j-1)*max_y+i;
        b_Wy = -L_ek(i,j,k)*L_ek(i-1,j,k)*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(i,j,k)*0.5*(y(i)-y(i-1))+L_ek(i-1,j,k)*0.5*(y(i+1)-y(i)));
        b_P =  b_Wy ;
        B(n,n) = -b_P + B(n,n);
        C(n) = b_P * Pn(NNCCG(w,1),NNCCG(w,2),NNCCG(w,3)) + C(n);
        a_Wy = -sigma(i,j,k)*sigma(i-1,j,k)*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(i,j,k)*0.5*(y(i)-y(i-1))+sigma(i-1,j,k)*0.5*(y(i+1)-y(i)));
        a_P =  a_Wy ;
        A(n,n) = -a_P + A(n,n);
        D(n) = a_P * SPn(NNCCG(w,1),NNCCG(w,2),NNCCG(w,3)) + D(n);
        
                
    elseif (NNCG(w) - HOSTNUM1(i,j,k)) == 6672 
        
%         i = NNCCL(w,1); j = NNCCL(w,2); k = NNCCL(w,3); 
        n=(k-1)*(max_y*max_x)+(j-1)*max_y+i;
        b_Wz = -L_ek(i,j,k)*L_ek(i,j,k-1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(L_ek(i,j,k)*0.5*(z(k)-z(k-1))+L_ek(i,j,k-1)*0.5*(z(k+1)-z(k)));
        b_P =  b_Wz ;
        B(n,n) = -b_P + B(n,n);
        C(n) = b_P * Pn(NNCCG(w,1),NNCCG(w,2),NNCCG(w,3)) + C(n);
        a_Wz = -sigma(i,j,k)*sigma(i,j,k-1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(sigma(i,j,k)*0.5*(z(k)-z(k-1))+sigma(i,j,k-1)*0.5*(z(k+1)-z(k)));
        a_P =  a_Wz ;
        A(n,n) = -a_P + A(n,n); 
        D(n) = a_P * SPn(NNCCG(w,1),NNCCG(w,2),NNCCG(w,3)) + D(n);
        
    elseif (HOSTNUM1(i,j,k)-NNCG(w)) == 139
        
%         i = NNCCL(w,1); j = NNCCL(w,2); k = NNCCL(w,3); 
        n=(k-1)*(max_y*max_x)+(j-1)*max_y+i;
        b_Ey = -L_ek(i,j,k)*L_ek(i+1,j,k)*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(i+1,j,k)*0.5*(y(i+1)-y(i))+L_ek(i,j,k)*0.5*(y(i+2)-y(i+1)));
        b_P =  b_Ey ;
        B(n,n) = -b_P + B(n,n);
        C(n) = b_P * Pn(NNCCG(w,1),NNCCG(w,2),NNCCG(w,3)) + C(n);
        a_Ey = -sigma(i,j,k)*sigma(i+1,j,k)*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(i+1,j,k)*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1)));
        a_P =  a_Ey ;
        A(n,n) = -a_P + A(n,n);
        D(n) = a_P * SPn(NNCCG(w,1),NNCCG(w,2),NNCCG(w,3)) + D(n);
        
    elseif (HOSTNUM1(i,j,k)-NNCG(w)) == 6672  
        
%         i = NNCCL(w,1); j = NNCCL(w,2); k = NNCCL(w,3); 
        n=(k-1)*(max_y*max_x)+(j-1)*max_y+i;
        b_Ez = -L_ek(i,j,k)*L_ek(i,j,k+1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(L_ek(i,j,k+1)*0.5*(z(k+1)-z(k))+L_ek(i,j,k)*0.5*(z(k+2)-z(k+1)));
        b_P =  b_Ez ;
        B(n,n) = -b_P + B(n,n);
        C(n) = b_P * Pn(NNCCG(w,1),NNCCG(w,2),NNCCG(w,3)) + C(n);
        a_Ez = -sigma(i,j,k)*sigma(i,j,k+1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(sigma(i,j,k+1)*0.5*(z(k+1)-z(k))+sigma(i,j,k)*0.5*(z(k+2)-z(k+1)));
        a_P =  a_Ez ;
        A(n,n) = -a_P + A(n,n);
        D(n) = a_P * SPn(NNCCG(w,1),NNCCG(w,2),NNCCG(w,3)) + D(n);
        
    elseif (HOSTNUM1(i,j,k)-NNCG(w)) == 1    
        
%         i = NNCCL(w,1); j = NNCCL(w,2); k = NNCCL(w,3); 
        n=(k-1)*(max_y*max_x)+(j-1)*max_y+i;
        b_Ex = -L_ek(i,j,k)*L_ek(i,j+1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(i,j+1,k)*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1)));
        b_P =  b_Ex ;
        B(n,n) = -b_P + B(n,n); 
        C(n) = b_P * Pn(NNCCG(w,1),NNCCG(w,2),NNCCG(w,3)) + C(n);
        a_Ex = -sigma(i,j,k)*sigma(i,j+1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(sigma(i,j+1,k)*0.5*(x(j+1)-x(j))+sigma(i,j,k)*0.5*(x(j+2)-x(j+1)));
        a_P =  a_Ex ;
        A(n,n) = -a_P + A(n,n); 
        D(n) = a_P * SPn(NNCCG(w,1),NNCCG(w,2),NNCCG(w,3)) + D(n);
    end
    
end








