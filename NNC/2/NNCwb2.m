[max_y, max_x, max_z]=size(sigma); 
% ne = max_y*max_x*max_z;
% B = sparse(ne,ne);
w = 1;
len = length(NNCI1b)+1;


while w < len
    
    
%     check = find(NNC.NNC2 == NNC.NNC2(w));
    check = find(NNCI2b == NNCI2b(w));
       
    w = check(end);
    wsize = size(check);
    
    NNCOV;
    
    i = NNCI2b_dex(w,1);
    j = NNCI2b_dex(w,2);
    k = NNCI2b_dex(w,3);  
    
    a = zeros(5,1);
    b = zeros(5,1);
    c = zeros(5,1);
    
        
    if wsize(1) > 1
        
        for wcheck = 1:wsize(1)
                    
            a(wcheck) = NNCI1b_dex(w,1);
            b(wcheck) = NNCI1b_dex(w,2);
            c(wcheck) = NNCI1b_dex(w,3); 
            w = w+1;            
        end
        
    else
            a(1) = NNCI1b_dex(w,1);
            b(1) = NNCI1b_dex(w,2);
            c(1) = NNCI1b_dex(w,3); 
            
            w = w+1;
    end
    
    n=(k-1)*(max_y*max_x)+(j-1)*max_y+i;
    nWx1 = (c(1)-1)*(max_y*max_x)+(b(1)-1)*max_y+a(1); 
    nWx2 = (c(2)-1)*(max_y*max_x)+(b(2)-1)*max_y+a(2);
    nWx2(nWx2<0)=n-1;
    nWx3 = (c(3)-1)*(max_y*max_x)+(b(3)-1)*max_y+a(3);
    nWx3(nWx3<0)=n-1;
    nWx4 = (c(4)-1)*(max_y*max_x)+(b(4)-1)*max_y+a(4);
    nWx4(nWx4<0)=n-1;
    nWx5 = (c(5)-1)*(max_y*max_x)+(b(5)-1)*max_y+a(5);
    nWx5(nWx5<0)=n-1;
                
                %----------------defines corners of resistance matrix ---------
                 
                
     if (j==1 && i==1) && k==1
%                     B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) = zeros;
                    b_Ex = -L_ek(i,j,k)*L_ek(i,j+1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(i,j+1,k)*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1)));
                    b_Ey = -L_ek(i,j,k)*L_ek(i+1,j,k)*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(i+1,j,k)*0.5*(y(i+1)-y(i))+L_ek(i,j,k)*0.5*(y(i+2)-y(i+1)));
                    b_Ez = -L_ek(i,j,k)*L_ek(i,j,k+1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(L_ek(i,j,k+1)*0.5*(z(k+1)-z(k))+L_ek(i,j,k)*0.5*(z(k+2)-z(k+1)));
                    b_P  =  -(b_Ex + (b_Ey+b_Ey2+b_Ey3) + b_Ez);
                    b_Ex(isnan(b_Ex)) = 0;
                    b_Ey(isnan(b_Ey)) = 0;
                    b_Ez(isnan(b_Ez)) = 0;
                    b_P(isnan(b_P)) = 0;
                B(n,[n (n+1) (n+max_y) (n+(max_y*max_x))]) = [b_P b_Ey b_Ex b_Ez];
                    
                elseif (j==max_x && i==max_y) && k==max_z
%                   B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) = zeros;
                    b_Wx1 = (-L_ek(i,j,k)*L_ek(a(1),b(1),c(1))*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(a(1),b(1),c(1))*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1))))/1;
                    if c(2) == 0
                        b_Wx2 = 0;
                    else
                        b_Wx2 = (-L_ek(i,j,k)*L_ek(a(2),b(2),c(2))*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(a(2),b(2),c(2))*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1))))/1;
                    end
                    if c(3) == 0 
                        b_Wx3 = 0;
                    else
                    b_Wx3 = (-L_ek(i,j,k)*L_ek(a(3),b(3),c(3))*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(a(3),b(3),c(3))*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1))))/1;
                    end
                    if c(4) == 0 
                        b_Wx4 = 0;
                    else
                    b_Wx4 = (-L_ek(i,j,k)*L_ek(a(4),b(4),c(4))*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(a(4),b(4),c(4))*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1))))/1;
                    end
                    if c(5) == 0 
                        b_Wx5 = 0;
                    else
                    b_Wx5 = (-L_ek(i,j,k)*L_ek(a(5),b(5),c(5))*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(a(5),b(5),c(5))*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1))))/1;
                    end
                    b_Wy = -wsize(1)*L_ek(i,j,k)*L_ek(i-1,j,k)*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(i,j,k)*0.5*(y(i)-y(i-1))+L_ek(i-1,j,k)*0.5*(y(i+1)-y(i)));
                    b_Wz = -wsize(1)*L_ek(i,j,k)*L_ek(i,j,k-1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(L_ek(i,j,k)*0.5*(z(k)-z(k-1))+L_ek(i,j,k-1)*0.5*(z(k+1)-z(k)));
                    b_P  =  -wsize(1)*(b_Wx1  +b_Wx2 +b_Wx3 +b_Wx4 +b_Wx5 +b_Wy + b_Wz);
                    b_P(isnan(b_P)) = 0;
                    b_Wz(isnan(b_Wz)) = 0;
                    b_Wx1(isnan(b_Wy1)) = 0;
                    b_Wx2(isnan(b_Wy2)) = 0;
                    b_Wx3(isnan(b_Wy3)) = 0;
                    b_Wx4(isnan(b_Wy4)) = 0;
                    b_Wx5(isnan(b_Wy5)) = 0;
                    b_Wy(isnan(b_Wx)) = 0;
                B(n,[n-(max_y*max_x) nWx1 nWx2 nWx3 nWx4 nWx5 (n-1) n]) = [b_Wz b_Wx1 b_Wx2 b_Wx3 b_Wx4 b_Wx5 b_Wy b_P];     
                    
                elseif (j==1 && i==1) && k==max_z
%                   B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) = zeros; 
                    b_Ex = -L_ek(i,j,k)*L_ek(i,j+1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(i,j+1,k)*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1)));
                    b_Ey = -L_ek(i,j,k)*L_ek(i+1,j,k)*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(i+1,j,k)*0.5*(y(i+1)-y(i))+L_ek(i,j,k)*0.5*(y(i+2)-y(i+1)));
                    b_Wz = -L_ek(i,j,k)*L_ek(i,j,k-1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(L_ek(i,j,k)*0.5*(z(k)-z(k-1))+L_ek(i,j,k-1)*0.5*(z(k+1)-z(k)));
                    b_P  =  -(b_Ex + b_Ey + b_Wz);
                    b_Ex(isnan(b_Ex)) = 0;
                    b_Ey(isnan(b_Ey)) = 0;
                    b_P(isnan(b_P)) = 0;
                    b_Wz(isnan(b_Wz)) = 0;
                B(n,[(n-(max_y*max_x)) n (n+1) (n+max_y)]) = [b_Wz b_P b_Ey b_Ex]; 
                
                elseif (j==max_x && i==1) && k==max_z
%                     B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) = zeros;
                    b_Wx = -L_ek(i,j,k)*L_ek(i,j-1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(i,j,k)*0.5*(x(j)-x(j-1))+L_ek(i,j-1,k)*0.5*(x(j+1)-x(j)));
                    b_Ey = -L_ek(i,j,k)*L_ek(i+1,j,k)*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(i+1,j,k)*0.5*(y(i+1)-y(i))+L_ek(i,j,k)*0.5*(y(i+2)-y(i+1)));
                    b_Wz = -L_ek(i,j,k)*L_ek(i,j,k-1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(L_ek(i,j,k)*0.5*(z(k)-z(k-1))+L_ek(i,j,k-1)*0.5*(z(k+1)-z(k)));
                    b_P  =  -(b_Wx + b_Ey + b_Wz);
                    b_Ey(isnan(b_Ey)) = 0;
                    b_P(isnan(b_P)) = 0;
                    b_Wz(isnan(b_Wz)) = 0;
                    b_Wx(isnan(b_Wx)) = 0;
                B(n,[(n-(max_y*max_x)) (n-max_y) n (n+1)]) = [b_Wz b_Wx b_P b_Ey];
                    
                elseif (j==max_x && i==1) && k==1
%                     B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) = zeros;
                    b_Wx = -L_ek(i,j,k)*L_ek(i,j-1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(i,j,k)*0.5*(x(j)-x(j-1))+L_ek(i,j-1,k)*0.5*(x(j+1)-x(j)));
                    b_Ey = -L_ek(i,j,k)*L_ek(i+1,j,k)*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(i+1,j,k)*0.5*(y(i+1)-y(i))+L_ek(i,j,k)*0.5*(y(i+2)-y(i+1)));
                    b_Ez = -L_ek(i,j,k)*L_ek(i,j,k+1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(L_ek(i,j,k+1)*0.5*(z(k+1)-z(k))+L_ek(i,j,k)*0.5*(z(k+2)-z(k+1)));
                    b_P  =  -(b_Wx + b_Ey + b_Ez);
                    b_Ey(isnan(b_Ey)) = 0;
                    b_Ez(isnan(b_Ez)) = 0;
                    b_P(isnan(b_P)) = 0;
                    b_Wx(isnan(b_Wx)) = 0;
                B(n,[(n-max_y) n (n+max_y) (n+(max_y*max_x))]) = [b_Wx b_P b_Ey b_Ez];
                    
                elseif (j==max_x && i==max_y) && k==1
%                     B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) = zeros;
                    b_Wx1 = (-L_ek(i,j,k)*L_ek(a(1),b(1),c(1))*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(a(1),b(1),c(1))*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1))))/1;
                    if c(2) == 0
                        b_Wx2 = 0;
                    else
                        b_Wx2 = (-L_ek(i,j,k)*L_ek(a(2),b(2),c(2))*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(a(2),b(2),c(2))*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1))))/1;
                    end
                    if c(3) == 0 
                        b_Wx3 = 0;
                    else
                    b_Wx3 = (-L_ek(i,j,k)*L_ek(a(3),b(3),c(3))*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(a(3),b(3),c(3))*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1))))/1;
                    end
                    if c(4) == 0 
                        b_Wx4 = 0;
                    else
                    b_Wx4 = (-L_ek(i,j,k)*L_ek(a(4),b(4),c(4))*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(a(4),b(4),c(4))*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1))))/1;
                    end
                    if c(5) == 0 
                        b_Wx5 = 0;
                    else
                    b_Wx5 = (-L_ek(i,j,k)*L_ek(a(5),b(5),c(5))*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(a(5),b(5),c(5))*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1))))/1;
                    end
                    b_Wy = -wsize(1)*L_ek(i,j,k)*L_ek(i-1,j,k)*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(i,j,k)*0.5*(y(i)-y(i-1))+L_ek(i-1,j,k)*0.5*(y(i+1)-y(i)));
                    b_Ez = -wsize(1)*L_ek(i,j,k)*L_ek(i,j,k+1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(L_ek(i,j,k+1)*0.5*(z(k+1)-z(k))+L_ek(i,j,k)*0.5*(z(k+2)-z(k+1)));
                    b_P  =  -wsize(1)*(b_Wx1  +b_Wx2 +b_Wx3 +b_Wx4 +b_Wx5 +b_Wy + b_Ez);
                    b_Ez(isnan(b_Ez)) = 0;
                    b_P(isnan(b_P)) = 0;
                    b_Wx1(isnan(b_Wy1)) = 0;
                    b_Wx2(isnan(b_Wy2)) = 0;
                    b_Wx3(isnan(b_Wy3)) = 0;
                    b_Wx4(isnan(b_Wy4)) = 0;
                    b_Wx5(isnan(b_Wy5)) = 0;
                    b_Wy(isnan(b_Wx)) = 0;
                B(n,[nWx1 nWx2 nWx3 nWx4 nWx5 (n-1) n (n+(max_y*max_x))]) = [b_Wx1 b_Wx2 b_Wx3 b_Wx4 b_Wx5 b_Wy b_P b_Ez];
                    
                elseif (j==1 && i==max_y) && k==max_z
%                     B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) = zeros;
                    b_Ex = -wsize(1)*L_ek(i,j,k)*L_ek(i,j+1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(i,j+1,k)*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1)));
                    b_Wy = -wsize(1)*L_ek(i,j,k)*L_ek(i-1,j,k)*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(i,j,k)*0.5*(y(i)-y(i-1))+L_ek(i-1,j,k)*0.5*(y(i+1)-y(i)));
                    b_Wz = -wsize(1)*L_ek(i,j,k)*L_ek(i,j,k-1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(L_ek(i,j,k)*0.5*(z(k)-z(k-1))+L_ek(i,j,k-1)*0.5*(z(k+1)-z(k)));
                    b_P  =  -wsize(1)*(b_Ex +b_Wy + b_Wz);
                    b_Ex(isnan(b_Ex)) = 0;
                    b_P(isnan(b_P)) = 0;
                    b_Wz(isnan(b_Wz)) = 0;
                    b_Wx1(isnan(b_Wy1)) = 0;
                    b_Wx2(isnan(b_Wy2)) = 0;
                    b_Wx3(isnan(b_Wy3)) = 0;
                    b_Wx4(isnan(b_Wy4)) = 0;
                    b_Wx5(isnan(b_Wy5)) = 0;
                    B(n,[(n-(max_y*max_x)) (n-1) n (n+max_y)]) = [b_Wz b_Wy b_P b_Ex];
                    
                elseif (j==1 && i==max_y) && k==1
%                     B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) = zeros;
                    b_Ex = -wsize(1)*L_ek(i,j,k)*L_ek(i,j+1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(i,j+1,k)*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1)));
                    b_Wy = -wsize(1)*L_ek(i,j,k)*L_ek(i-1,j,k)*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(i,j,k)*0.5*(y(i)-y(i-1))+L_ek(i-1,j,k)*0.5*(y(i+1)-y(i)));
                    b_Ez = -wsize(1)*L_ek(i,j,k)*L_ek(i,j,k+1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(L_ek(i,j,k+1)*0.5*(z(k+1)-z(k))+L_ek(i,j,k)*0.5*(z(k+2)-z(k+1)));
                    b_P  =  -wsize(1)*(b_Ex +b_Wy + b_Ez);
                    b_Ex(isnan(b_Ex)) = 0;
                    b_Ez(isnan(b_Ez)) = 0;
                    b_P(isnan(b_P)) = 0;
                    
                B(n,[(n-1) n (n+max_y) (n+(max_y*max_x))]) = [b_Wy1  b_Wy b_P b_Ex b_Ez];
                    
                    %----------------defines edges of resistance matrix -----------
                elseif (i==1 && k==1)
%                     B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) = zeros;
                    b_Ex = -L_ek(i,j,k)*L_ek(i,j+1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(i,j+1,k)*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1)));
                    b_Wx = -L_ek(i,j,k)*L_ek(i,j-1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(i,j,k)*0.5*(x(j)-x(j-1))+L_ek(i,j-1,k)*0.5*(x(j+1)-x(j)));
                    b_Ey = -L_ek(i,j,k)*L_ek(i+1,j,k)*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(i+1,j,k)*0.5*(y(i+1)-y(i))+L_ek(i,j,k)*0.5*(y(i+2)-y(i+1)));
                    b_Ez = -L_ek(i,j,k)*L_ek(i,j,k+1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(L_ek(i,j,k+1)*0.5*(z(k+1)-z(k))+L_ek(i,j,k)*0.5*(z(k+2)-z(k+1)));
                    b_P  =  -(b_Ex + b_Wx + b_Ey + b_Ez);
                    b_Ex(isnan(b_Ex)) = 0;
                    b_Ey(isnan(b_Ey)) = 0;
                    b_Ez(isnan(b_Ez)) = 0;
                    b_P(isnan(b_P)) = 0;
                B(n,[(n-max_y) n (n+1) (n+max_y) (n+(max_y*max_x))]) = [b_Wx b_P b_Ey b_Ex b_Ez];
                
                elseif (i==1 && j==1)
%                     B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) = zeros;
                    b_Ex = -L_ek(i,j,k)*L_ek(i,j+1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(i,j+1,k)*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1)));
                    b_Ey = -L_ek(i,j,k)*L_ek(i+1,j,k)*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(i+1,j,k)*0.5*(y(i+1)-y(i))+L_ek(i,j,k)*0.5*(y(i+2)-y(i+1)));
                    b_Ez = -L_ek(i,j,k)*L_ek(i,j,k+1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(L_ek(i,j,k+1)*0.5*(z(k+1)-z(k))+L_ek(i,j,k)*0.5*(z(k+2)-z(k+1)));
                    b_Wz = -L_ek(i,j,k)*L_ek(i,j,k-1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(L_ek(i,j,k)*0.5*(z(k)-z(k-1))+L_ek(i,j,k-1)*0.5*(z(k+1)-z(k)));
                    b_P  =  -(b_Ex + b_Ey + b_Ez + b_Wz);
                    b_Ex(isnan(b_Ex)) = 0;
                    b_Ey(isnan(b_Ey)) = 0;
                    b_Ez(isnan(b_Ez)) = 0;
                    b_P(isnan(b_P)) = 0;
                    b_Wz(isnan(b_Wz)) = 0;
                B(n,[(n-(max_y*max_x)) n (n+1) (n+max_y) (n+(max_y*max_x))]) = [b_Wz b_P b_Ey b_Ex b_Ez];
                    
                elseif (i==1 && k==max_z)
%                     B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) = zeros;
                    b_Ex = -L_ek(i,j,k)*L_ek(i,j+1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(i,j+1,k)*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1)));
                    b_Wx = -L_ek(i,j,k)*L_ek(i,j-1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(i,j,k)*0.5*(x(j)-x(j-1))+L_ek(i,j-1,k)*0.5*(x(j+1)-x(j)));
                    b_Ey = -L_ek(i,j,k)*L_ek(i+1,j,k)*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(i+1,j,k)*0.5*(y(i+1)-y(i))+L_ek(i,j,k)*0.5*(y(i+2)-y(i+1)));
                    b_Wz = -L_ek(i,j,k)*L_ek(i,j,k-1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(L_ek(i,j,k)*0.5*(z(k)-z(k-1))+L_ek(i,j,k-1)*0.5*(z(k+1)-z(k)));
                    b_P  =  -(b_Ex + b_Wx + b_Ey + b_Wz);
                    b_Ex(isnan(b_Ex)) = 0;
                    b_Ey(isnan(b_Ey)) = 0;
                    b_P(isnan(b_P)) = 0;
                    b_Wz(isnan(b_Wz)) = 0;
                    b_Wx(isnan(b_Wx)) = 0;
                B(n,[(n-(max_y*max_x)) (n-max_y) n (n+1) (n+max_y)]) = [b_Wz b_Wx b_P b_Ey b_Ex];
                    
                elseif (i==1 && j==max_x)
%                     B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) = zeros;
                    b_Wx = -L_ek(i,j,k)*L_ek(i,j-1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(i,j,k)*0.5*(x(j)-x(j-1))+L_ek(i,j-1,k)*0.5*(x(j+1)-x(j)));
                    b_Ey = -L_ek(i,j,k)*L_ek(i+1,j,k)*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(i+1,j,k)*0.5*(y(i+1)-y(i))+L_ek(i,j,k)*0.5*(y(i+2)-y(i+1)));
                    b_Ez = -L_ek(i,j,k)*L_ek(i,j,k+1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(L_ek(i,j,k+1)*0.5*(z(k+1)-z(k))+L_ek(i,j,k)*0.5*(z(k+2)-z(k+1)));
                    b_Wz = -L_ek(i,j,k)*L_ek(i,j,k-1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(L_ek(i,j,k)*0.5*(z(k)-z(k-1))+L_ek(i,j,k-1)*0.5*(z(k+1)-z(k)));
                    b_P  =  -(b_Wx + b_Ey + b_Ez + b_Wz);
                    b_Ey(isnan(b_Ey)) = 0;
                    b_Ez(isnan(b_Ez)) = 0;
                    b_P(isnan(b_P)) = 0;
                    b_Wz(isnan(b_Wz)) = 0;
                    b_Wx(isnan(b_Wx)) = 0;
                B(n,[(n-(max_y*max_x)) (n-max_y) n (n+1) (n+(max_y*max_x))]) = [b_Wz b_Wx1 b_Wx2 b_Wx3 b_Wx4 b_Wx5 b_P b_Ey b_Ez];
                   
                elseif (i==max_y && k==1)
%                     B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) = zeros;
                    b_Ex = -wsize(1)*L_ek(i,j,k)*L_ek(i,j+1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(i,j+1,k)*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1)));
                    b_Wx1 = (-L_ek(i,j,k)*L_ek(a(1),b(1),c(1))*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(a(1),b(1),c(1))*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1))))/1;
                    if c(2) == 0
                        b_Wx2 = 0;
                    else
                        b_Wx2 = (-L_ek(i,j,k)*L_ek(a(2),b(2),c(2))*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(a(2),b(2),c(2))*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1))))/1;
                    end
                    if c(3) == 0 
                        b_Wx3 = 0;
                    else
                    b_Wx3 = (-L_ek(i,j,k)*L_ek(a(3),b(3),c(3))*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(a(3),b(3),c(3))*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1))))/1;
                    end
                    if c(4) == 0 
                        b_Wx4 = 0;
                    else
                    b_Wx4 = (-L_ek(i,j,k)*L_ek(a(4),b(4),c(4))*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(a(4),b(4),c(4))*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1))))/1;
                    end
                    if c(5) == 0 
                        b_Wx5 = 0;
                    else
                    b_Wx5 = (-L_ek(i,j,k)*L_ek(a(5),b(5),c(5))*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(a(5),b(5),c(5))*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1))))/1;
                    end
                    b_Wy = -wsize(1)*L_ek(i,j,k)*L_ek(i-1,j,k)*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(i,j,k)*0.5*(y(i)-y(i-1))+L_ek(i-1,j,k)*0.5*(y(i+1)-y(i)));
                    b_Ez = -wsize(1)*L_ek(i,j,k)*L_ek(i,j,k+1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(L_ek(i,j,k+1)*0.5*(z(k+1)-z(k))+L_ek(i,j,k)*0.5*(z(k+2)-z(k+1)));
                    b_P  =  -wsize(1)*(b_Ex + b_Wx1  +b_Wx2 +b_Wx3 +b_Wx4 +b_Wx5 +b_Wy + b_Ez);
                    b_Ex(isnan(b_Ex)) = 0;
                    b_Ez(isnan(b_Ez)) = 0;
                    b_P(isnan(b_P)) = 0;
                    b_Wx1(isnan(b_Wy1)) = 0;
                    b_Wx2(isnan(b_Wy2)) = 0;
                    b_Wx3(isnan(b_Wy3)) = 0;
                    b_Wx4(isnan(b_Wy4)) = 0;
                    b_Wx5(isnan(b_Wy5)) = 0;
                    b_Wy(isnan(b_Wx)) = 0;
                B(n,[nWx1 nWx2 nWx3 nWx4 nWx5 (n-1) n (n+max_y) (n+(max_y*max_x))]) = [b_Wx1 b_Wx2 b_Wx3 b_Wx4 b_Wx5 b_Wy b_P b_Ex b_Ez];
                    
                elseif (i==max_y && j==max_x)
%                     B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) = zeros;
                    b_Wx1 = (-L_ek(i,j,k)*L_ek(a(1),b(1),c(1))*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(a(1),b(1),c(1))*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1))))/1;
                    if c(2) == 0
                        b_Wx2 = 0;
                    else
                        b_Wx2 = (-L_ek(i,j,k)*L_ek(a(2),b(2),c(2))*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(a(2),b(2),c(2))*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1))))/1;
                    end
                    if c(3) == 0 
                        b_Wx3 = 0;
                    else
                    b_Wx3 = (-L_ek(i,j,k)*L_ek(a(3),b(3),c(3))*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(a(3),b(3),c(3))*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1))))/1;
                    end
                    if c(4) == 0 
                        b_Wx4 = 0;
                    else
                    b_Wx4 = (-L_ek(i,j,k)*L_ek(a(4),b(4),c(4))*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(a(4),b(4),c(4))*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1))))/1;
                    end
                    if c(5) == 0 
                        b_Wx5 = 0;
                    else
                    b_Wx5 = (-L_ek(i,j,k)*L_ek(a(5),b(5),c(5))*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(a(5),b(5),c(5))*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1))))/1;
                    end
                    b_Wy = -wsize(1)*L_ek(i,j,k)*L_ek(i-1,j,k)*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(i,j,k)*0.5*(y(i)-y(i-1))+L_ek(i-1,j,k)*0.5*(y(i+1)-y(i)));
                    b_Ez = -wsize(1)*L_ek(i,j,k)*L_ek(i,j,k+1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(L_ek(i,j,k+1)*0.5*(z(k+1)-z(k))+L_ek(i,j,k)*0.5*(z(k+2)-z(k+1)));
                    b_Wz = -wsize(1)*L_ek(i,j,k)*L_ek(i,j,k-1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(L_ek(i,j,k)*0.5*(z(k)-z(k-1))+L_ek(i,j,k-1)*0.5*(z(k+1)-z(k)));
                    b_P  =  -wsize(1)*(b_Wx1  +b_Wx2 +b_Wx3 +b_Wx4 +b_Wx5 +b_Wy + b_Ez + b_Wz);
                    b_Ez(isnan(b_Ez)) = 0;
                    b_P(isnan(b_P)) = 0;
                    b_Wz(isnan(b_Wz)) = 0;
                    b_Wx1(isnan(b_Wy1)) = 0;
                    b_Wx2(isnan(b_Wy2)) = 0;
                    b_Wx3(isnan(b_Wy3)) = 0;
                    b_Wx4(isnan(b_Wy4)) = 0;
                    b_Wx5(isnan(b_Wy5)) = 0;
                    b_Wy(isnan(b_Wx)) = 0;
                B(n,[(n-(max_y*max_x)) nWx1 nWx2 nWx3 nWx4 nWx5 (n-1) n (n+(max_y*max_x))]) = [b_Wz b_Wx1 b_Wx2 b_Wx3 b_Wx4 b_Wx5 b_Wy b_P b_Ez];
                    
                elseif (i==max_y && j==1)
%                     B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) = zeros;
                    b_Ex = -wsize(1)*L_ek(i,j,k)*L_ek(i,j+1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(i,j+1,k)*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1)));
                    b_Wy1 = (-L_ek(i,j,k)*L_ek(a(1),b(1),c(1))*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(a(1),b(1),c(1))*0.5*(y(i+1)-y(i))+L_ek(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    b_Wy = -wsize(1)*L_ek(i,j,k)*L_ek(i-1,j,k)*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(i,j,k)*0.5*(y(i)-y(i-1))+L_ek(i-1,j,k)*0.5*(y(i+1)-y(i)));
                    b_Ez = -wsize(1)*L_ek(i,j,k)*L_ek(i,j,k+1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(L_ek(i,j,k+1)*0.5*(z(k+1)-z(k))+L_ek(i,j,k)*0.5*(z(k+2)-z(k+1)));
                    b_Wz = -wsize(1)*L_ek(i,j,k)*L_ek(i,j,k-1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(L_ek(i,j,k)*0.5*(z(k)-z(k-1))+L_ek(i,j,k-1)*0.5*(z(k+1)-z(k)));
                    b_P  =  -wsize(1)*(b_Ex +b_Wy + b_Ez + b_Wz);
                    b_Ex(isnan(b_Ex)) = 0;
                    b_Ez(isnan(b_Ez)) = 0;
                    b_P(isnan(b_P)) = 0;
                    b_Wz(isnan(b_Wz)) = 0;
                    b_Wx1(isnan(b_Wy1)) = 0;
                    b_Wx2(isnan(b_Wy2)) = 0;
                    b_Wx3(isnan(b_Wy3)) = 0;
                    b_Wx4(isnan(b_Wy4)) = 0;
                    b_Wx5(isnan(b_Wy5)) = 0;
                    b_Wy(isnan(b_Wx)) = 0;
                B(n,[(n-(max_y*max_x)) (n-1) n (n+max_y) (n+(max_y*max_x))]) = [b_Wz b_Wy b_P b_Ex b_Ez];
                    
                elseif (i==max_y && k==max_z)
%                     B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) = zeros;
                    b_Ex = -wsize(1)*L_ek(i,j,k)*L_ek(i,j+1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(i,j+1,k)*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1)));
                    b_Wx1 = (-L_ek(i,j,k)*L_ek(a(1),b(1),c(1))*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(a(1),b(1),c(1))*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1))))/1;
                    if c(2) == 0
                        b_Wx2 = 0;
                    else
                        b_Wx2 = (-L_ek(i,j,k)*L_ek(a(2),b(2),c(2))*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(a(2),b(2),c(2))*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1))))/1;
                    end
                    if c(3) == 0 
                        b_Wx3 = 0;
                    else
                    b_Wx3 = (-L_ek(i,j,k)*L_ek(a(3),b(3),c(3))*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(a(3),b(3),c(3))*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1))))/1;
                    end
                    if c(4) == 0 
                        b_Wx4 = 0;
                    else
                    b_Wx4 = (-L_ek(i,j,k)*L_ek(a(4),b(4),c(4))*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(a(4),b(4),c(4))*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1))))/1;
                    end
                    if c(5) == 0 
                        b_Wx5 = 0;
                    else
                    b_Wx5 = (-L_ek(i,j,k)*L_ek(a(5),b(5),c(5))*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(a(5),b(5),c(5))*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1))))/1;
                    end
                    b_Wy = -wsize(1)*L_ek(i,j,k)*L_ek(i-1,j,k)*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(i,j,k)*0.5*(y(i)-y(i-1))+L_ek(i-1,j,k)*0.5*(y(i+1)-y(i)));
                    b_Wz = -L_ek(i,j,k)*L_ek(i,j,k-1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(L_ek(i,j,k)*0.5*(z(k)-z(k-1))+L_ek(i,j,k-1)*0.5*(z(k+1)-z(k)));
                    b_P  =  -wsize(1)*(b_Wx1  +b_Wx2 +b_Wx3 +b_Wx4 +b_Wx5 + b_Wx +b_Wy+ b_Wz);
                    b_Ex(isnan(b_Ex)) = 0;
                    b_P(isnan(b_P)) = 0;
                    b_Wz(isnan(b_Wz)) = 0;
                    b_Wx1(isnan(b_Wy1)) = 0;
                    b_Wx2(isnan(b_Wy2)) = 0;
                    b_Wx3(isnan(b_Wy3)) = 0;
                    b_Wx4(isnan(b_Wy4)) = 0;
                    b_Wx5(isnan(b_Wy5)) = 0;
                    b_Wy(isnan(b_Wx)) = 0;
                B(n,[(n-(max_y*max_x)) nWx1 nWx2 nWx3 nWx4 nWx5 (n-1) n (n+max_y)]) = [b_Wz b_Wx1 b_Wx2 b_Wx3 b_Wx4 b_Wx5 b_Wy b_Wy5 b_P b_Ex];
                    
                    
                elseif (j==max_x && k==1)
%                     B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) = zeros;
                    b_Wx1 = (-L_ek(i,j,k)*L_ek(a(1),b(1),c(1))*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(a(1),b(1),c(1))*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1))))/1;
                    if c(2) == 0
                        b_Wx2 = 0;
                    else
                        b_Wx2 = (-L_ek(i,j,k)*L_ek(a(2),b(2),c(2))*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(a(2),b(2),c(2))*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1))))/1;
                    end
                    if c(3) == 0 
                        b_Wx3 = 0;
                    else
                    b_Wx3 = (-L_ek(i,j,k)*L_ek(a(3),b(3),c(3))*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(a(3),b(3),c(3))*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1))))/1;
                    end
                    if c(4) == 0 
                        b_Wx4 = 0;
                    else
                    b_Wx4 = (-L_ek(i,j,k)*L_ek(a(4),b(4),c(4))*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(a(4),b(4),c(4))*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1))))/1;
                    end
                    if c(5) == 0 
                        b_Wx5 = 0;
                    else
                    b_Wx5 = (-L_ek(i,j,k)*L_ek(a(5),b(5),c(5))*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(a(5),b(5),c(5))*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1))))/1;
                    end
                    b_Ey = -wsize(1)*L_ek(i,j,k)*L_ek(i+1,j,k)*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(i+1,j,k)*0.5*(y(i+1)-y(i))+L_ek(i,j,k)*0.5*(y(i+2)-y(i+1)));
                    b_Wy = -wsize(1)*L_ek(i,j,k)*L_ek(i-1,j,k)*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(i,j,k)*0.5*(y(i)-y(i-1))+L_ek(i-1,j,k)*0.5*(y(i+1)-y(i)));
                    b_Ez = -wsize(1)*L_ek(i,j,k)*L_ek(i,j,k+1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(L_ek(i,j,k+1)*0.5*(z(k+1)-z(k))+L_ek(i,j,k)*0.5*(z(k+2)-z(k+1)));
                    b_P  =  -wsize(1)*(b_Wx1  +b_Wx2 +b_Wx3 +b_Wx4 +b_Wx5 + b_Ey +b_Wy + b_Ez);
                    b_Ey(isnan(b_Ey)) = 0;
                    b_Ez(isnan(b_Ez)) = 0;
                    b_P(isnan(b_P)) = 0;
                    b_Wx1(isnan(b_Wy1)) = 0;
                    b_Wx2(isnan(b_Wy2)) = 0;
                    b_Wx3(isnan(b_Wy3)) = 0;
                    b_Wx4(isnan(b_Wy4)) = 0;
                    b_Wx5(isnan(b_Wy5)) = 0;
                    b_Wy(isnan(b_Wx)) = 0;
                B(n,[nWx1 nWx2 nWx3 nWx4 nWx5 (n-1) n (n+1) (n+(max_y*max_x))]) = [b_Wx1 b_Wx2 b_Wx3 b_Wx4 b_Wx5 b_Wy b_P b_Ey b_Ez];
                    
                elseif (j==max_x && k==max_z)
%                     B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) = zeros;
                    b_Wx1 = (-L_ek(i,j,k)*L_ek(a(1),b(1),c(1))*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(a(1),b(1),c(1))*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1))))/1;
                    if c(2) == 0
                        b_Wx2 = 0;
                    else
                        b_Wx2 = (-L_ek(i,j,k)*L_ek(a(2),b(2),c(2))*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(a(2),b(2),c(2))*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1))))/1;
                    end
                    if c(3) == 0 
                        b_Wx3 = 0;
                    else
                    b_Wx3 = (-L_ek(i,j,k)*L_ek(a(3),b(3),c(3))*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(a(3),b(3),c(3))*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1))))/1;
                    end
                    if c(4) == 0 
                        b_Wx4 = 0;
                    else
                    b_Wx4 = (-L_ek(i,j,k)*L_ek(a(4),b(4),c(4))*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(a(4),b(4),c(4))*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1))))/1;
                    end
                    if c(5) == 0 
                        b_Wx5 = 0;
                    else
                    b_Wx5 = (-L_ek(i,j,k)*L_ek(a(5),b(5),c(5))*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(a(5),b(5),c(5))*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1))))/1;
                    end
                    b_Ey = -wsize(1)*L_ek(i,j,k)*L_ek(i+1,j,k)*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(i+1,j,k)*0.5*(y(i+1)-y(i))+L_ek(i,j,k)*0.5*(y(i+2)-y(i+1)));
                    b_Wy = -wsize(1)*L_ek(i,j,k)*L_ek(i-1,j,k)*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(i,j,k)*0.5*(y(i)-y(i-1))+L_ek(i-1,j,k)*0.5*(y(i+1)-y(i)));
                    b_Wz = -wsize(1)*L_ek(i,j,k)*L_ek(i,j,k-1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(L_ek(i,j,k)*0.5*(z(k)-z(k-1))+L_ek(i,j,k-1)*0.5*(z(k+1)-z(k)));
                    b_P  =  -wsize(1)*(b_Wx1  +b_Wx2 +b_Wx3 +b_Wx4 +b_Wx5 + b_Ey +b_Wy + b_Wz);
                    b_Ey(isnan(b_Ey)) = 0;
                    b_P(isnan(b_P)) = 0;
                    b_Wz(isnan(b_Wz)) = 0;
                    b_Wx1(isnan(b_Wy1)) = 0;
                    b_Wx2(isnan(b_Wy2)) = 0;
                    b_Wx3(isnan(b_Wy3)) = 0;
                    b_Wx4(isnan(b_Wy4)) = 0;
                    b_Wx5(isnan(b_Wy5)) = 0;
                    b_Wy(isnan(b_Wx)) = 0;
                B(n,[(n-(max_y*max_x)) nWx1 nWx2 nWx3 nWx4 nWx5 (n-1) n (n+1)]) = [b_Wz b_Wx1 b_Wx2 b_Wx3 b_Wx4 b_Wx5 b_Wy b_P b_Ey];
                    
                elseif (j==1 && k==max_z)
%                     B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) = zeros;
                    b_Ex = -wsize(1)*L_ek(i,j,k)*L_ek(i,j+1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(i,j+1,k)*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1)));
                    b_Ey = -wsize(1)*L_ek(i,j,k)*L_ek(i+1,j,k)*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(i+1,j,k)*0.5*(y(i+1)-y(i))+L_ek(i,j,k)*0.5*(y(i+2)-y(i+1)));
                    b_Wy = -wsize(1)*L_ek(i,j,k)*L_ek(i-1,j,k)*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(i,j,k)*0.5*(y(i)-y(i-1))+L_ek(i-1,j,k)*0.5*(y(i+1)-y(i)));
                    b_Wz = -wsize(1)*L_ek(i,j,k)*L_ek(i,j,k-1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(L_ek(i,j,k)*0.5*(z(k)-z(k-1))+L_ek(i,j,k-1)*0.5*(z(k+1)-z(k)));
                    b_P  =  -wsize(1)*(b_Ex + b_Ey +b_Wy + b_Wz);
                    b_Ex(isnan(b_Ex)) = 0;
                    b_Ey(isnan(b_Ey)) = 0;
                    b_P(isnan(b_P)) = 0;
                    b_Wz(isnan(b_Wz)) = 0;
                    b_Wy1(isnan(b_Wy1)) = 0;
                    b_Wy2(isnan(b_Wy2)) = 0;
                    b_Wy3(isnan(b_Wy3)) = 0;
                    b_Wy4(isnan(b_Wy4)) = 0;
                    b_Wy5(isnan(b_Wy5)) = 0;
                B(n,[(n-(max_y*max_x)) (n-1) n (n+1) (n+max_y)]) = [b_Wz b_Wy b_P b_Ey b_Ex];
                    
                elseif (j==1 && k==1)
%                     B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) = zeros;
                    b_Ex = -wsize(1)*L_ek(i,j,k)*L_ek(i,j+1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(i,j+1,k)*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1)));
                    b_Ey = -wsize(1)*L_ek(i,j,k)*L_ek(i+1,j,k)*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(i+1,j,k)*0.5*(y(i+1)-y(i))+L_ek(i,j,k)*0.5*(y(i+2)-y(i+1)));
                    b_Wy = -wsize(1)*L_ek(i,j,k)*L_ek(i-1,j,k)*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(i,j,k)*0.5*(y(i)-y(i-1))+L_ek(i-1,j,k)*0.5*(y(i+1)-y(i)));
                    b_Ez = -wsize(1)*L_ek(i,j,k)*L_ek(i,j,k+1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(L_ek(i,j,k+1)*0.5*(z(k+1)-z(k))+L_ek(i,j,k)*0.5*(z(k+2)-z(k+1)));
                    b_P  =  -wsize(1)*(b_Ex + b_Ey +b_Wy + b_Ez);
                    b_Ex(isnan(b_Ex)) = 0;
                    b_Ey(isnan(b_Ey)) = 0;
                    b_Ez(isnan(b_Ez)) = 0;
                    b_P(isnan(b_P)) = 0;
                    b_Wx1(isnan(b_Wy1)) = 0;
                    b_Wx2(isnan(b_Wy2)) = 0;
                    b_Wx3(isnan(b_Wy3)) = 0;
                    b_Wx4(isnan(b_Wy4)) = 0;
                    b_Wx5(isnan(b_Wy5)) = 0;
                    b_Wy(isnan(b_Wx)) = 0;
                B(n,[(n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) = [b_Wy b_P b_Ey b_Ex b_Ez];
                    
                    %----------------this part defines the surfaces
                    %defining the surface boundaries---------------
                elseif i==1
%                     B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) = zeros;
                    b_Ex = -L_ek(i,j,k)*L_ek(i,j+1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(i,j+1,k)*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1)));
                    b_Wx = -L_ek(i,j,k)*L_ek(i,j-1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(i,j,k)*0.5*(x(j)-x(j-1))+L_ek(i,j-1,k)*0.5*(x(j+1)-x(j)));
                    b_Ey = -L_ek(i,j,k)*L_ek(i+1,j,k)*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(i+1,j,k)*0.5*(y(i+1)-y(i))+L_ek(i,j,k)*0.5*(y(i+2)-y(i+1)));
                    b_Ez = -L_ek(i,j,k)*L_ek(i,j,k+1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(L_ek(i,j,k+1)*0.5*(z(k+1)-z(k))+L_ek(i,j,k)*0.5*(z(k+2)-z(k+1)));
                    b_Wz = -L_ek(i,j,k)*L_ek(i,j,k-1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(L_ek(i,j,k)*0.5*(z(k)-z(k-1))+L_ek(i,j,k-1)*0.5*(z(k+1)-z(k)));
                    b_P  =  -(b_Ex + b_Wx + b_Ey + b_Ez + b_Wz);
                    b_Ex(isnan(b_Ex)) = 0;
                    b_Ey(isnan(b_Ey)) = 0;
                    b_Ez(isnan(b_Ez)) = 0;
                    b_P(isnan(b_P)) = 0;
                    b_Wz(isnan(b_Wz)) = 0;
                    b_Wx(isnan(b_Wx)) = 0;
                B(n,[(n-(max_y*max_x)) (n-max_y) n (n+1) (n+max_y) (n+(max_y*max_x))]) = [b_Wz b_Wx b_P b_Ey b_Ex b_Ez];
                    
                elseif i==max_y
%                     B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) = zeros;
                    b_Ex = -wsize(1)*L_ek(i,j,k)*L_ek(i,j+1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(i,j+1,k)*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1)));
                    b_Wx1 = (-L_ek(i,j,k)*L_ek(a(1),b(1),c(1))*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(a(1),b(1),c(1))*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1))))/1;
                    if c(2) == 0
                        b_Wx2 = 0;
                    else
                        b_Wx2 = (-L_ek(i,j,k)*L_ek(a(2),b(2),c(2))*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(a(2),b(2),c(2))*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1))))/1;
                    end
                    if c(3) == 0 
                        b_Wx3 = 0;
                    else
                    b_Wx3 = (-L_ek(i,j,k)*L_ek(a(3),b(3),c(3))*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(a(3),b(3),c(3))*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1))))/1;
                    end
                    if c(4) == 0 
                        b_Wx4 = 0;
                    else
                    b_Wx4 = (-L_ek(i,j,k)*L_ek(a(4),b(4),c(4))*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(a(4),b(4),c(4))*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1))))/1;
                    end
                    if c(5) == 0 
                        b_Wx5 = 0;
                    else
                    b_Wx5 = (-L_ek(i,j,k)*L_ek(a(5),b(5),c(5))*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(a(5),b(5),c(5))*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1))))/1;
                    end
                    b_Wy = -wsize(1)*L_ek(i,j,k)*L_ek(i-1,j,k)*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(i,j,k)*0.5*(y(i)-y(i-1))+L_ek(i-1,j,k)*0.5*(y(i+1)-y(i)));
                    b_Ez = -wsize(1)*L_ek(i,j,k)*L_ek(i,j,k+1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(L_ek(i,j,k+1)*0.5*(z(k+1)-z(k))+L_ek(i,j,k)*0.5*(z(k+2)-z(k+1)));
                    b_Wz = -wsize(1)*L_ek(i,j,k)*L_ek(i,j,k-1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(L_ek(i,j,k)*0.5*(z(k)-z(k-1))+L_ek(i,j,k-1)*0.5*(z(k+1)-z(k)));
                    b_P  =  -wsize(1)*(b_Ex + b_Wx1  +b_Wx2 +b_Wx3 +b_Wx4 +b_Wx5 +b_Wy + b_Ez + b_Wz);
                    b_Ex(isnan(b_Ex)) = 0;
                    b_Ez(isnan(b_Ez)) = 0;
                    b_P(isnan(b_P)) = 0;
                    b_Wz(isnan(b_Wz)) = 0;
                    b_Wx1(isnan(b_Wy1)) = 0;
                    b_Wx2(isnan(b_Wy2)) = 0;
                    b_Wx3(isnan(b_Wy3)) = 0;
                    b_Wx4(isnan(b_Wy4)) = 0;
                    b_Wx5(isnan(b_Wy5)) = 0;
                    b_Wy(isnan(b_Wx)) = 0;
                B(n,[(n-(max_y*max_x)) nWx1 nWx2 nWx3 nWx4 nWx5 (n-1) n (n+max_y) (n+(max_y*max_x))]) = [b_Wz b_Wx1 b_Wx2 b_Wx3 b_Wx4 b_Wx5 b_Wy b_P b_Ex b_Ez];
                    
                elseif j==1
%                     B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) = zeros;
                    b_Ex = -wsize(1)*L_ek(i,j,k)*L_ek(i,j+1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(i,j+1,k)*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1)));
                    b_Ey = -wsize(1)*L_ek(i,j,k)*L_ek(i+1,j,k)*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(i+1,j,k)*0.5*(y(i+1)-y(i))+L_ek(i,j,k)*0.5*(y(i+2)-y(i+1)));
                    b_Wy = -wsize(1)*L_ek(i,j,k)*L_ek(i-1,j,k)*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(i,j,k)*0.5*(y(i)-y(i-1))+L_ek(i-1,j,k)*0.5*(y(i+1)-y(i)));
                    b_Ez = -wsize(1)*L_ek(i,j,k)*L_ek(i,j,k+1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(L_ek(i,j,k+1)*0.5*(z(k+1)-z(k))+L_ek(i,j,k)*0.5*(z(k+2)-z(k+1)));
                    b_Wz = -wsize(1)*L_ek(i,j,k)*L_ek(i,j,k-1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(L_ek(i,j,k)*0.5*(z(k)-z(k-1))+L_ek(i,j,k-1)*0.5*(z(k+1)-z(k)));
                    b_P  =  -wsize(1)*(b_Ex + b_Ey +b_Wy + b_Ez + b_Wz);
                    b_Ex(isnan(b_Ex)) = 0;
                    b_Ey(isnan(b_Ey)) = 0;
                    b_Ez(isnan(b_Ez)) = 0;
                    b_P(isnan(b_P)) = 0;
                    b_Wz(isnan(b_Wz)) = 0;
                    b_Wx1(isnan(b_Wy1)) = 0;
                    b_Wx2(isnan(b_Wy2)) = 0;
                    b_Wx3(isnan(b_Wy3)) = 0;
                    b_Wx4(isnan(b_Wy4)) = 0;
                    b_Wx5(isnan(b_Wy5)) = 0;
                    b_Wy(isnan(b_Wx)) = 0;
                B(n,[(n-(max_y*max_x)) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) = [b_Wz b_Wy b_P b_Ey b_Ex b_Ez];
                    
                elseif j==max_x
%                     B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) = zeros;
                    b_Wx1 = (-L_ek(i,j,k)*L_ek(a(1),b(1),c(1))*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(a(1),b(1),c(1))*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1))))/1;
                    if c(2) == 0
                        b_Wx2 = 0;
                    else
                        b_Wx2 = (-L_ek(i,j,k)*L_ek(a(2),b(2),c(2))*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(a(2),b(2),c(2))*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1))))/1;
                    end
                    if c(3) == 0 
                        b_Wx3 = 0;
                    else
                    b_Wx3 = (-L_ek(i,j,k)*L_ek(a(3),b(3),c(3))*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(a(3),b(3),c(3))*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1))))/1;
                    end
                    if c(4) == 0 
                        b_Wx4 = 0;
                    else
                    b_Wx4 = (-L_ek(i,j,k)*L_ek(a(4),b(4),c(4))*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(a(4),b(4),c(4))*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1))))/1;
                    end
                    if c(5) == 0 
                        b_Wx5 = 0;
                    else
                    b_Wx5 = (-L_ek(i,j,k)*L_ek(a(5),b(5),c(5))*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(a(5),b(5),c(5))*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1))))/1;
                    end
                    b_Ey = -wsize(1)*L_ek(i,j,k)*L_ek(i+1,j,k)*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(i+1,j,k)*0.5*(y(i+1)-y(i))+L_ek(i,j,k)*0.5*(y(i+2)-y(i+1)));
                    b_Wy = -wsize(1)*L_ek(i,j,k)*L_ek(i-1,j,k)*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(i,j,k)*0.5*(y(i)-y(i-1))+L_ek(i-1,j,k)*0.5*(y(i+1)-y(i)));
                    b_Ez = -wsize(1)*L_ek(i,j,k)*L_ek(i,j,k+1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(L_ek(i,j,k+1)*0.5*(z(k+1)-z(k))+L_ek(i,j,k)*0.5*(z(k+2)-z(k+1)));
                    b_Wz = -wsize(1)*L_ek(i,j,k)*L_ek(i,j,k-1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(L_ek(i,j,k)*0.5*(z(k)-z(k-1))+L_ek(i,j,k-1)*0.5*(z(k+1)-z(k)));
                    b_P  =  -wsize(1)*(b_Wx1  +b_Wx2 +b_Wx3 +b_Wx4 +b_Wx5 + b_Ey +b_Wy + b_Ez + b_Wz);
                    b_Ey(isnan(b_Ey)) = 0;
                    b_Ez(isnan(b_Ez)) = 0;
                    b_P(isnan(b_P)) = 0;
                    b_Wz(isnan(b_Wz)) = 0;
                    b_Wx1(isnan(b_Wy1)) = 0;
                    b_Wx2(isnan(b_Wy2)) = 0;
                    b_Wx3(isnan(b_Wy3)) = 0;
                    b_Wx4(isnan(b_Wy4)) = 0;
                    b_Wx5(isnan(b_Wy5)) = 0;
                    b_Wy(isnan(b_Wx)) = 0;
                B(n,[(n-(max_y*max_x)) nWx1 nWx2 nWx3 nWx4 nWx5 (n-1) n (n+1) (n+(max_y*max_x))]) = [b_Wz b_Wx1 b_Wx2 b_Wx3 b_Wx4 b_Wx5 b_Wy b_P b_Ey b_Ez];
                    
                elseif k==1
%                     B(n,[(n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) = zeros;
                    b_Ex = -wsize(1)*L_ek(i,j,k)*L_ek(i,j+1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(i,j+1,k)*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1)));
                    b_Wx1 = (-L_ek(i,j,k)*L_ek(a(1),b(1),c(1))*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(a(1),b(1),c(1))*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1))))/1;
                    if c(2) == 0
                        b_Wx2 = 0;
                    else
                        b_Wx2 = (-L_ek(i,j,k)*L_ek(a(2),b(2),c(2))*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(a(2),b(2),c(2))*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1))))/1;
                    end
                    if c(3) == 0 
                        b_Wx3 = 0;
                    else
                    b_Wx3 = (-L_ek(i,j,k)*L_ek(a(3),b(3),c(3))*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(a(3),b(3),c(3))*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1))))/1;
                    end
                    if c(4) == 0 
                        b_Wx4 = 0;
                    else
                    b_Wx4 = (-L_ek(i,j,k)*L_ek(a(4),b(4),c(4))*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(a(4),b(4),c(4))*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1))))/1;
                    end
                    if c(5) == 0 
                        b_Wx5 = 0;
                    else
                    b_Wx5 = (-L_ek(i,j,k)*L_ek(a(5),b(5),c(5))*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(a(5),b(5),c(5))*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1))))/1;
                    end
                    b_Ey = -wsize(1)*L_ek(i,j,k)*L_ek(i+1,j,k)*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(i+1,j,k)*0.5*(y(i+1)-y(i))+L_ek(i,j,k)*0.5*(y(i+2)-y(i+1)));
                    b_Wy = -wsize(1)*L_ek(i,j,k)*L_ek(i-1,j,k)*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(i,j,k)*0.5*(y(i)-y(i-1))+L_ek(i-1,j,k)*0.5*(y(i+1)-y(i)));
                    b_Ez = -wsize(1)*L_ek(i,j,k)*L_ek(i,j,k+1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(L_ek(i,j,k+1)*0.5*(z(k+1)-z(k))+L_ek(i,j,k)*0.5*(z(k+2)-z(k+1)));
                    b_P  =  -wsize(1)*(b_Ex +b_Wx1  +b_Wx2 +b_Wx3 +b_Wx4 +b_Wx5 + b_Ey +b_Wy + b_Ez);
                    b_Ex(isnan(b_Ex)) = 0;
                    b_Ey(isnan(b_Ey)) = 0;
                    b_Ez(isnan(b_Ez)) = 0;
                    b_P(isnan(b_P)) = 0;
                    b_Wx1(isnan(b_Wy1)) = 0;
                    b_Wx2(isnan(b_Wy2)) = 0;
                    b_Wx3(isnan(b_Wy3)) = 0;
                    b_Wx4(isnan(b_Wy4)) = 0;
                    b_Wx5(isnan(b_Wy5)) = 0;
                    b_Wy(isnan(b_Wx)) = 0;
                B(n,[nWx1 nWx2 nWx3 nWx4 nWx5 (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) = [b_Wx1 b_Wx2 b_Wx3 b_Wx4 b_Wx5 b_Wy b_P b_Ey b_Ex b_Ez];
                
                elseif k==max_z
%                     B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) = zeros;
                    b_Ex = -wsize(1)*L_ek(i,j,k)*L_ek(i,j+1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(i,j+1,k)*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1)));
                    b_Wx1 = (-L_ek(i,j,k)*L_ek(a(1),b(1),c(1))*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(a(1),b(1),c(1))*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1))))/1;
                    if c(2) == 0
                        b_Wx2 = 0;
                    else
                        b_Wx2 = (-L_ek(i,j,k)*L_ek(a(2),b(2),c(2))*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(a(2),b(2),c(2))*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1))))/1;
                    end
                    if c(3) == 0 
                        b_Wx3 = 0;
                    else
                    b_Wx3 = (-L_ek(i,j,k)*L_ek(a(3),b(3),c(3))*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(a(3),b(3),c(3))*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1))))/1;
                    end
                    if c(4) == 0 
                        b_Wx4 = 0;
                    else
                    b_Wx4 = (-L_ek(i,j,k)*L_ek(a(4),b(4),c(4))*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(a(4),b(4),c(4))*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1))))/1;
                    end
                    if c(5) == 0 
                        b_Wx5 = 0;
                    else
                    b_Wx5 = (-L_ek(i,j,k)*L_ek(a(5),b(5),c(5))*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(a(5),b(5),c(5))*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1))))/1;
                    end
                    b_Ey = -wsize(1)*L_ek(i,j,k)*L_ek(i+1,j,k)*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(i+1,j,k)*0.5*(y(i+1)-y(i))+L_ek(i,j,k)*0.5*(y(i+2)-y(i+1)));
                    b_Wy = -wsize(1)*L_ek(i,j,k)*L_ek(i-1,j,k)*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(i,j,k)*0.5*(y(i)-y(i-1))+L_ek(i-1,j,k)*0.5*(y(i+1)-y(i)));
                    b_Wz = -wsize(1)*L_ek(i,j,k)*L_ek(i,j,k-1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(L_ek(i,j,k)*0.5*(z(k)-z(k-1))+L_ek(i,j,k-1)*0.5*(z(k+1)-z(k)));
                    b_P  =  -wsize(1)*(b_Ex + b_Wx1  +b_Wx2 +b_Wx3 +b_Wx4 +b_Wx5 + b_Ey +b_Wy + b_Wz);
                    b_Ex(isnan(b_Ex)) = 0;
                    b_Ey(isnan(b_Ey)) = 0;
                    b_P(isnan(b_P)) = 0; 
                    b_Wz(isnan(b_Wz)) = 0;
                    b_Wx1(isnan(b_Wy1)) = 0;
                    b_Wx2(isnan(b_Wy2)) = 0;
                    b_Wx3(isnan(b_Wy3)) = 0;
                    b_Wx4(isnan(b_Wy4)) = 0;
                    b_Wx5(isnan(b_Wy5)) = 0;
                    b_Wy(isnan(b_Wx)) = 0;
                B(n,[(n-(max_y*max_x)) nWx1 nWx2 nWx3 nWx4 nWx5 (n-1) n (n+1) (n+max_y)]) = [b_Wz b_Wx1 b_Wx2 b_Wx3 b_Wx4 b_Wx5 b_Wy b_P b_Ey b_Ex];
                
               elseif (i>1 && j>1 && k>1) && (i<max_y && j<max_x && k<max_z)
%                    B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) = zeros;
                    b_Ex = -wsize(1)*L_ek(i,j,k)*L_ek(i,j+1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(i,j+1,k)*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1)));
                    b_Ey = -wsize(1)*L_ek(i,j,k)*L_ek(i+1,j,k)*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(i+1,j,k)*0.5*(y(i+1)-y(i))+L_ek(i,j,k)*0.5*(y(i+2)-y(i+1)));
                    
                    b_Wx1 = (-L_ek(i,j,k)*L_ek(a(1),b(1),c(1))*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(a(1),b(1),c(1))*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1))))/1;
                    if c(2) == 0
                        b_Wx2 = 0;
                    else
                        b_Wx2 = (-L_ek(i,j,k)*L_ek(a(2),b(2),c(2))*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(a(2),b(2),c(2))*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1))))/1;
                    end
                    if c(3) == 0 
                        b_Wx3 = 0;
                    else
                    b_Wx3 = (-L_ek(i,j,k)*L_ek(a(3),b(3),c(3))*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(a(3),b(3),c(3))*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1))))/1;
                    end
                    if c(4) == 0 
                        b_Wx4 = 0;
                    else
                    b_Wx4 = (-L_ek(i,j,k)*L_ek(a(4),b(4),c(4))*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(a(4),b(4),c(4))*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1))))/1;
                    end
                    if c(5) == 0 
                        b_Wx5 = 0;
                    else
                    b_Wx5 = (-L_ek(i,j,k)*L_ek(a(5),b(5),c(5))*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(a(5),b(5),c(5))*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1))))/1;
                    end
                    b_Wy = -wsize(1)*L_ek(i,j,k)*L_ek(i-1,j,k)*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(i,j,k)*0.5*(y(i)-y(i-1))+L_ek(i-1,j,k)*0.5*(y(i+1)-y(i)));
                    b_Ez = -wsize(1)*L_ek(i,j,k)*L_ek(i,j,k+1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(L_ek(i,j,k+1)*0.5*(z(k+1)-z(k))+L_ek(i,j,k)*0.5*(z(k+2)-z(k+1)));
                    b_Wz = -wsize(1)*L_ek(i,j,k)*L_ek(i,j,k-1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(L_ek(i,j,k)*0.5*(z(k)-z(k-1))+L_ek(i,j,k-1)*0.5*(z(k+1)-z(k)));
                    b_P  =  -wsize(1)*(b_Ex + b_Wy + b_Ey +b_Wx1  +b_Wx2 +b_Wx3 +b_Wx4 +b_Wx5 + b_Ez + b_Wz);
                    b_Ex(isnan(b_Ex)) = 0;
                    b_Ey(isnan(b_Ey)) = 0;
                    b_Ez(isnan(b_Ez)) = 0;
                    b_P(isnan(b_P)) = 0;
                    b_Wz(isnan(b_Wz)) = 0;
                    b_Wx1(isnan(b_Wy1)) = 0;
                    b_Wx2(isnan(b_Wy2)) = 0;
                    b_Wx3(isnan(b_Wy3)) = 0;
                    b_Wx4(isnan(b_Wy4)) = 0;
                    b_Wx5(isnan(b_Wy5)) = 0;
                    b_Wy(isnan(b_Wx)) = 0;
                    
                B(n,[(n-(max_y*max_x)) nWx1 nWx2 nWx3 nWx4 nWx5 (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) = [b_Wz b_Wx1 b_Wx2 b_Wx3 b_Wx4 b_Wx5  b_Wy b_P b_Ey b_Ex b_Ez];
                    
                end           
    
    
end
