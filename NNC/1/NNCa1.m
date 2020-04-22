[max_y, max_x, max_z]=size(sigma); 
w = 1;
len = length(NNC1)+1;
while w < len
    
%     NNCC1(n,:) = [ i j k ];
%     
%     i = NNCC1(n,1);
%     j = NNCC1(n,2);
%     k = NNCC1(n,3);
    
    check = find(NNC.NNC1 == NNC.NNC1(w));
       
    w = check(end);
    wsize = size(check);
    
    NNCO;
    
    i = NNCC1(w,1);
    j = NNCC1(w,2);
    k = NNCC1(w,3);  
    
    a = zeros(3,1);
    b = zeros(3,1);
    c = zeros(3,1);
    
        
    if wsize(1) > 1
        
        for wcheck = 1:wsize(1)
                    
            a(wcheck) = NNCC2(w,1);
            b(wcheck) = NNCC2(w,2);
            c(wcheck) = NNCC2(w,3); 
            w = w+1;            
        end
        
    else
            a(1) = NNCC2(w,1);
            b(1) = NNCC2(w,2);
            c(1) = NNCC2(w,3); 
            
            w = w+1;
    end
    
    n=(k-1)*(max_y*max_x)+(j-1)*max_y+i;
    nEy1 = (c(1)-1)*(max_y*max_x)+(b(1)-1)*max_y+a(1); 
    nEy2 = (c(2)-1)*(max_y*max_x)+(b(2)-1)*max_y+a(2);
    nEy2(nEy2<0)=n+1;
    nEy3 = (c(3)-1)*(max_y*max_x)+(b(3)-1)*max_y+a(3);
    nEy3(nEy3<0)=n+1;
                
                %----------------defines corners of resistance matrix ---------
                 
                
     if (j==1 && i==1) && k==1
%                  B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) =zeros;
                    b_Ex = -L_ek(i,j,k)*L_ek(i,j+1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(i,j+1,k)*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1)));
                    b_Ey1 = (-L_ek(i,j,k)*L_ek(a(1),b(1),c(1))*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(a(1),b(1),c(1))*0.5*(y(i+1)-y(i))+L_ek(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    if c(2) == 0
                        b_Ey2 = 0;
                    else
                        b_Ey2 = (-L_ek(i,j,k)*L_ek(a(2),b(2),c(2))*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(a(2),b(2),c(2))*0.5*(y(i+1)-y(i))+L_ek(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    if c(3) == 0 
                        b_Ey3 = 0;
                    else
                    b_Ey3 = (-L_ek(i,j,k)*L_ek(a(3),b(3),c(3))*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(a(3),b(3),c(3))*0.5*(y(i+1)-y(i))+L_ek(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    b_Ez = -L_ek(i,j,k)*L_ek(i,j,k+1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(L_ek(i,j,k+1)*0.5*(z(k+1)-z(k))+L_ek(i,j,k)*0.5*(z(k+2)-z(k+1)));
                    b_P  =  -(b_Ex + (b_Ey1+b_Ey2+b_Ey3) + b_Ez);
                    b_Ex(isnan(b_Ex)) = 0;
                    b_Ey1(isnan(b_Ey1)) = 0;
                    b_Ey2(isnan(b_Ey2)) = 0;
                    b_Ey3(isnan(b_Ey3)) = 0;
                    b_Ez(isnan(b_Ez)) = 0;
                    b_P(isnan(b_P)) = 0;
                B(n,[n (nEy1) (nEy2) (nEy3) (n+max_y) (n+(max_y*max_x))]) = [b_P b_Ey1 b_Ey2 b_Ey3 b_Ex b_Ez];
                    
                elseif (j==max_x && i==max_y) && k==max_z
%                   B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) =zeros;  
                    b_Wx = -L_ek(i,j,k)*L_ek(i,j-1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(i,j,k)*0.5*(x(j)-x(j-1))+L_ek(i,j-1,k)*0.5*(x(j+1)-x(j)));
                    b_Wy = -L_ek(i,j,k)*L_ek(i-1,j,k)*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(i,j,k)*0.5*(y(i)-y(i-1))+L_ek(i-1,j,k)*0.5*(y(i+1)-y(i)));
                    b_Wz = -L_ek(i,j,k)*L_ek(i,j,k-1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(L_ek(i,j,k)*0.5*(z(k)-z(k-1))+L_ek(i,j,k-1)*0.5*(z(k+1)-z(k)));
                    b_P  =  -(b_Wx + b_Wy + b_Wz);
                    b_P(isnan(b_P)) = 0;
                    b_Wz(isnan(b_Wz)) = 0;
                    b_Wy(isnan(b_Wy)) = 0;
                    b_Wx(isnan(b_Wx)) = 0;
                B(n,[n-(max_y*max_x) n-max_y n-1 n]) = [b_Wz b_Wx b_Wy b_P];     
                    
                elseif (j==1 && i==1) && k==max_z
%                   B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) =zeros; 
                    b_Ex = -L_ek(i,j,k)*L_ek(i,j+1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(i,j+1,k)*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1)));
                    b_Ey1 = (-L_ek(i,j,k)*L_ek(a(1),b(1),c(1))*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(a(1),b(1),c(1))*0.5*(y(i+1)-y(i))+L_ek(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    if c(2) == 0
                        b_Ey2 = 0;
                    else
                        b_Ey2 = (-L_ek(i,j,k)*L_ek(a(2),b(2),c(2))*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(a(2),b(2),c(2))*0.5*(y(i+1)-y(i))+L_ek(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    if c(3) == 0 
                        b_Ey3 = 0;
                    else
                    b_Ey3 = (-L_ek(i,j,k)*L_ek(a(3),b(3),c(3))*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(a(3),b(3),c(3))*0.5*(y(i+1)-y(i))+L_ek(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    b_Wz = -L_ek(i,j,k)*L_ek(i,j,k-1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(L_ek(i,j,k)*0.5*(z(k)-z(k-1))+L_ek(i,j,k-1)*0.5*(z(k+1)-z(k)));
                    b_P  =  -(b_Ex + b_Ey1 +b_Ey2 +b_Ey3 + b_Wz);
                    b_Ex(isnan(b_Ex)) = 0;
                    b_Ey1(isnan(b_Ey1)) = 0;
                    b_Ey2(isnan(b_Ey2)) = 0;
                    b_Ey3(isnan(b_Ey3)) = 0;
                    b_P(isnan(b_P)) = 0;
                    b_Wz(isnan(b_Wz)) = 0;
                B(n,[(n-(max_y*max_x)) n (nEy1) (nEy2) (nEy3) (n+max_y)]) = [b_Wz b_P b_Ey1 b_Ey2 b_Ey3 b_Ex]; 
                
                elseif (j==max_x && i==1) && k==max_z
%                    B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) =zeros;
                    b_Wx = -L_ek(i,j,k)*L_ek(i,j-1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(i,j,k)*0.5*(x(j)-x(j-1))+L_ek(i,j-1,k)*0.5*(x(j+1)-x(j)));
                    b_Ey1 = (-L_ek(i,j,k)*L_ek(a(1),b(1),c(1))*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(a(1),b(1),c(1))*0.5*(y(i+1)-y(i))+L_ek(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    if c(2) == 0
                        b_Ey2 = 0;
                    else
                        b_Ey2 = (-L_ek(i,j,k)*L_ek(a(2),b(2),c(2))*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(a(2),b(2),c(2))*0.5*(y(i+1)-y(i))+L_ek(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    if c(3) == 0 
                        b_Ey3 = 0;
                    else
                    b_Ey3 = (-L_ek(i,j,k)*L_ek(a(3),b(3),c(3))*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(a(3),b(3),c(3))*0.5*(y(i+1)-y(i))+L_ek(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    b_Wz = -L_ek(i,j,k)*L_ek(i,j,k-1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(L_ek(i,j,k)*0.5*(z(k)-z(k-1))+L_ek(i,j,k-1)*0.5*(z(k+1)-z(k)));
                    b_P  =  -(b_Wx + b_Ey1 +b_Ey2 +b_Ey3 + b_Wz);
                    b_Ey1(isnan(b_Ey1)) = 0;
                    b_Ey2(isnan(b_Ey2)) = 0;
                    b_Ey3(isnan(b_Ey3)) = 0;
                    b_P(isnan(b_P)) = 0;
                    b_Wz(isnan(b_Wz)) = 0;
                    b_Wx(isnan(b_Wx)) = 0;
                B(n,[(n-(max_y*max_x)) (n-max_y) n (nEy1) (nEy2) (nEy3)]) = [b_Wz b_Wx b_P b_Ey1 b_Ey2 b_Ey3];
                    
                elseif (j==max_x && i==1) && k==1
%                     B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) =zeros;
                    b_Wx = -L_ek(i,j,k)*L_ek(i,j-1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(i,j,k)*0.5*(x(j)-x(j-1))+L_ek(i,j-1,k)*0.5*(x(j+1)-x(j)));
                    b_Ey1 = (-L_ek(i,j,k)*L_ek(a(1),b(1),c(1))*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(a(1),b(1),c(1))*0.5*(y(i+1)-y(i))+L_ek(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    if c(2) == 0
                        b_Ey2 = 0;
                    else
                        b_Ey2 = (-L_ek(i,j,k)*L_ek(a(2),b(2),c(2))*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(a(2),b(2),c(2))*0.5*(y(i+1)-y(i))+L_ek(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    if c(3) == 0 
                        b_Ey3 = 0;
                    else
                    b_Ey3 = (-L_ek(i,j,k)*L_ek(a(3),b(3),c(3))*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(a(3),b(3),c(3))*0.5*(y(i+1)-y(i))+L_ek(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    b_Ez = -L_ek(i,j,k)*L_ek(i,j,k+1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(L_ek(i,j,k+1)*0.5*(z(k+1)-z(k))+L_ek(i,j,k)*0.5*(z(k+2)-z(k+1)));
                    b_P  =  -(b_Wx + b_Ey1 +b_Ey2 +b_Ey3 + b_Ez);
                    b_Ey1(isnan(b_Ey1)) = 0;
                    b_Ey2(isnan(b_Ey2)) = 0;
                    b_Ey3(isnan(b_Ey3)) = 0;
                    b_Ez(isnan(b_Ez)) = 0;
                    b_P(isnan(b_P)) = 0;
                    b_Wx(isnan(b_Wx)) = 0;
                B(n,[(n-max_y) n (n+max_y) (n+(max_y*max_x))]) = [b_Wx b_P b_Ey1 b_Ey2 b_Ey3 b_Ez];
                    
                elseif (j==max_x && i==max_y) && k==1
%                     B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) =zeros;
                    b_Wx = -L_ek(i,j,k)*L_ek(i,j-1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(i,j,k)*0.5*(x(j)-x(j-1))+L_ek(i,j-1,k)*0.5*(x(j+1)-x(j)));
                    b_Wy = -L_ek(i,j,k)*L_ek(i-1,j,k)*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(i,j,k)*0.5*(y(i)-y(i-1))+L_ek(i-1,j,k)*0.5*(y(i+1)-y(i)));
                    b_Ez = -L_ek(i,j,k)*L_ek(i,j,k+1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(L_ek(i,j,k+1)*0.5*(z(k+1)-z(k))+L_ek(i,j,k)*0.5*(z(k+2)-z(k+1)));
                    b_P  =  -(b_Wx + b_Wy + b_Ez);
                    b_Ez(isnan(b_Ez)) = 0;
                    b_P(isnan(b_P)) = 0;
                    b_Wy(isnan(b_Wy)) = 0;
                    b_Wx(isnan(b_Wx)) = 0;
                B(n,[(n-max_y) (n-1) n (n+(max_y*max_x))]) = [b_Wx b_Wy b_P b_Ez];
                    
                elseif (j==1 && i==max_y) && k==max_z
%                     B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) =zeros;
                    b_Ex = -L_ek(i,j,k)*L_ek(i,j+1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(i,j+1,k)*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1)));
                    b_Wy = -L_ek(i,j,k)*L_ek(i-1,j,k)*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(i,j,k)*0.5*(y(i)-y(i-1))+L_ek(i-1,j,k)*0.5*(y(i+1)-y(i)));
                    b_Wz = -L_ek(i,j,k)*L_ek(i,j,k-1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(L_ek(i,j,k)*0.5*(z(k)-z(k-1))+L_ek(i,j,k-1)*0.5*(z(k+1)-z(k)));
                    b_P  =  -(b_Ex + b_Wy + b_Wz);
                    b_Ex(isnan(b_Ex)) = 0;
                    b_P(isnan(b_P)) = 0;
                    b_Wz(isnan(b_Wz)) = 0;
                    b_Wy(isnan(b_Wy)) = 0;
                    B(n,[(n-(max_y*max_x)) (n-1) n (n+max_y)]) = [b_Wz b_Wy b_P b_Ex];
                    
                elseif (j==1 && i==max_y) && k==1
%                     B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) =zeros;
                    b_Ex = -L_ek(i,j,k)*L_ek(i,j+1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(i,j+1,k)*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1)));
                    b_Wy = -L_ek(i,j,k)*L_ek(i-1,j,k)*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(i,j,k)*0.5*(y(i)-y(i-1))+L_ek(i-1,j,k)*0.5*(y(i+1)-y(i)));
                    b_Ez = -L_ek(i,j,k)*L_ek(i,j,k+1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(L_ek(i,j,k+1)*0.5*(z(k+1)-z(k))+L_ek(i,j,k)*0.5*(z(k+2)-z(k+1)));
                    b_P  =  -(b_Ex + b_Wy + b_Ez);
                    b_Ex(isnan(b_Ex)) = 0;
                    b_Ez(isnan(b_Ez)) = 0;
                    b_P(isnan(b_P)) = 0;
                    b_Wy(isnan(b_Wy)) = 0;
                B(n,[(n-1) n (n+max_y) (n+(max_y*max_x))]) = [b_Wy b_P b_Ex b_Ez];
                    
                    %----------------defines edges of resistance matrix -----------
                elseif (i==1 && k==1)
%                     B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) =zeros;
                    b_Ex = -L_ek(i,j,k)*L_ek(i,j+1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(i,j+1,k)*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1)));
                    b_Wx = -L_ek(i,j,k)*L_ek(i,j-1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(i,j,k)*0.5*(x(j)-x(j-1))+L_ek(i,j-1,k)*0.5*(x(j+1)-x(j)));
                    b_Ey1 = (-L_ek(i,j,k)*L_ek(a(1),b(1),c(1))*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(a(1),b(1),c(1))*0.5*(y(i+1)-y(i))+L_ek(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    if c(2) == 0
                        b_Ey2 = 0;
                    else
                        b_Ey2 = (-L_ek(i,j,k)*L_ek(a(2),b(2),c(2))*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(a(2),b(2),c(2))*0.5*(y(i+1)-y(i))+L_ek(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    if c(3) == 0 
                        b_Ey3 = 0;
                    else
                    b_Ey3 = (-L_ek(i,j,k)*L_ek(a(3),b(3),c(3))*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(a(3),b(3),c(3))*0.5*(y(i+1)-y(i))+L_ek(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    b_Ez = -L_ek(i,j,k)*L_ek(i,j,k+1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(L_ek(i,j,k+1)*0.5*(z(k+1)-z(k))+L_ek(i,j,k)*0.5*(z(k+2)-z(k+1)));
                    b_P  =  -(b_Ex + b_Wx + b_Ey1 +b_Ey2 +b_Ey3 + b_Ez);
                    b_Ex(isnan(b_Ex)) = 0;
                    b_Ey1(isnan(b_Ey1)) = 0;
                    b_Ey2(isnan(b_Ey2)) = 0;
                    b_Ey3(isnan(b_Ey3)) = 0;
                    b_Ez(isnan(b_Ez)) = 0;
                    b_P(isnan(b_P)) = 0;
                B(n,[(n-max_y) n (nEy1) (nEy2) (nEy3) (n+max_y) (n+(max_y*max_x))]) = [b_Wx b_P b_Ey1 b_Ey2 b_Ey3 b_Ex b_Ez];
                
                elseif (i==1 && j==1)
%                     B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) =zeros;
                    b_Ex = -L_ek(i,j,k)*L_ek(i,j+1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(i,j+1,k)*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1)));
                    b_Ey1 = (-L_ek(i,j,k)*L_ek(a(1),b(1),c(1))*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(a(1),b(1),c(1))*0.5*(y(i+1)-y(i))+L_ek(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    if c(2) == 0
                        b_Ey2 = 0;
                    else
                        b_Ey2 = (-L_ek(i,j,k)*L_ek(a(2),b(2),c(2))*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(a(2),b(2),c(2))*0.5*(y(i+1)-y(i))+L_ek(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    if c(3) == 0 
                        b_Ey3 = 0;
                    else
                    b_Ey3 = (-L_ek(i,j,k)*L_ek(a(3),b(3),c(3))*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(a(3),b(3),c(3))*0.5*(y(i+1)-y(i))+L_ek(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    b_Ez = -L_ek(i,j,k)*L_ek(i,j,k+1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(L_ek(i,j,k+1)*0.5*(z(k+1)-z(k))+L_ek(i,j,k)*0.5*(z(k+2)-z(k+1)));
                    b_Wz = -L_ek(i,j,k)*L_ek(i,j,k-1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(L_ek(i,j,k)*0.5*(z(k)-z(k-1))+L_ek(i,j,k-1)*0.5*(z(k+1)-z(k)));
                    b_P  =  -(b_Ex + b_Ey1 +b_Ey2 +b_Ey3 + b_Ez + b_Wz);
                    b_Ex(isnan(b_Ex)) = 0;
                    b_Ey1(isnan(b_Ey1)) = 0;
                    b_Ey2(isnan(b_Ey2)) = 0;
                    b_Ey3(isnan(b_Ey3)) = 0;
                    b_Ez(isnan(b_Ez)) = 0;
                    b_P(isnan(b_P)) = 0;
                    b_Wz(isnan(b_Wz)) = 0;
                B(n,[(n-(max_y*max_x)) n (nEy1) (nEy2) (nEy3) (n+max_y) (n+(max_y*max_x))]) = [b_Wz b_P b_Ey1 b_Ey2 b_Ey3 b_Ex b_Ez];
                    
                elseif (i==1 && k==max_z)
%                     B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) =zeros;
                    b_Ex = -L_ek(i,j,k)*L_ek(i,j+1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(i,j+1,k)*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1)));
                    b_Wx = -L_ek(i,j,k)*L_ek(i,j-1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(i,j,k)*0.5*(x(j)-x(j-1))+L_ek(i,j-1,k)*0.5*(x(j+1)-x(j)));
                    b_Ey1 = (-L_ek(i,j,k)*L_ek(a(1),b(1),c(1))*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(a(1),b(1),c(1))*0.5*(y(i+1)-y(i))+L_ek(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    if c(2) == 0
                        b_Ey2 = 0;
                    else
                        b_Ey2 = (-L_ek(i,j,k)*L_ek(a(2),b(2),c(2))*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(a(2),b(2),c(2))*0.5*(y(i+1)-y(i))+L_ek(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    if c(3) == 0 
                        b_Ey3 = 0;
                    else
                    b_Ey3 = (-L_ek(i,j,k)*L_ek(a(3),b(3),c(3))*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(a(3),b(3),c(3))*0.5*(y(i+1)-y(i))+L_ek(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    b_Wz = -L_ek(i,j,k)*L_ek(i,j,k-1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(L_ek(i,j,k)*0.5*(z(k)-z(k-1))+L_ek(i,j,k-1)*0.5*(z(k+1)-z(k)));
                    b_P  =  -(b_Ex + b_Wx + b_Ey1 +b_Ey2 +b_Ey3 + b_Wz);
                    b_Ex(isnan(b_Ex)) = 0;
                    b_Ey1(isnan(b_Ey1)) = 0;
                    b_Ey2(isnan(b_Ey2)) = 0;
                    b_Ey3(isnan(b_Ey3)) = 0;
                    b_P(isnan(b_P)) = 0;
                    b_Wz(isnan(b_Wz)) = 0;
                    b_Wx(isnan(b_Wx)) = 0;
                B(n,[(n-(max_y*max_x)) (n-max_y) n (nEy1) (nEy2) (nEy3) (n+max_y)]) = [b_Wz b_Wx b_P b_Ey1 b_Ey2 b_Ey3 b_Ex];
                    
                elseif (i==1 && j==max_x)
%                     B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) =zeros;
                    b_Wx = -L_ek(i,j,k)*L_ek(i,j-1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(i,j,k)*0.5*(x(j)-x(j-1))+L_ek(i,j-1,k)*0.5*(x(j+1)-x(j)));
                    b_Ey1 = (-L_ek(i,j,k)*L_ek(a(1),b(1),c(1))*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(a(1),b(1),c(1))*0.5*(y(i+1)-y(i))+L_ek(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    if c(2) == 0
                        b_Ey2 = 0;
                    else
                        b_Ey2 = (-L_ek(i,j,k)*L_ek(a(2),b(2),c(2))*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(a(2),b(2),c(2))*0.5*(y(i+1)-y(i))+L_ek(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    if c(3) == 0 
                        b_Ey3 = 0;
                    else
                    b_Ey3 = (-L_ek(i,j,k)*L_ek(a(3),b(3),c(3))*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(a(3),b(3),c(3))*0.5*(y(i+1)-y(i))+L_ek(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    b_Ez = -L_ek(i,j,k)*L_ek(i,j,k+1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(L_ek(i,j,k+1)*0.5*(z(k+1)-z(k))+L_ek(i,j,k)*0.5*(z(k+2)-z(k+1)));
                    b_Wz = -L_ek(i,j,k)*L_ek(i,j,k-1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(L_ek(i,j,k)*0.5*(z(k)-z(k-1))+L_ek(i,j,k-1)*0.5*(z(k+1)-z(k)));
                    b_P  =  -(b_Wx + b_Ey1 +b_Ey2 +b_Ey3 + b_Ez + b_Wz);
                    b_Ey1(isnan(b_Ey1)) = 0;
                    b_Ey2(isnan(b_Ey2)) = 0;
                    b_Ey3(isnan(b_Ey3)) = 0;
                    b_Ez(isnan(b_Ez)) = 0;
                    b_P(isnan(b_P)) = 0;
                    b_Wz(isnan(b_Wz)) = 0;
                    b_Wx(isnan(b_Wx)) = 0;
                B(n,[(n-(max_y*max_x)) (n-max_y) n (nEy1) (nEy2) (nEy3) (n+(max_y*max_x))]) = [b_Wz b_Wx b_P b_Ey1 b_Ey2 b_Ey3 b_Ez];
                   
                elseif (i==max_y && k==1)
%                     B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) =zeros;
                    b_Ex = -L_ek(i,j,k)*L_ek(i,j+1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(i,j+1,k)*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1)));
                    b_Wx = -L_ek(i,j,k)*L_ek(i,j-1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(i,j,k)*0.5*(x(j)-x(j-1))+L_ek(i,j-1,k)*0.5*(x(j+1)-x(j)));
                    b_Wy = -L_ek(i,j,k)*L_ek(i-1,j,k)*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(i,j,k)*0.5*(y(i)-y(i-1))+L_ek(i-1,j,k)*0.5*(y(i+1)-y(i)));
                    b_Ez = -L_ek(i,j,k)*L_ek(i,j,k+1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(L_ek(i,j,k+1)*0.5*(z(k+1)-z(k))+L_ek(i,j,k)*0.5*(z(k+2)-z(k+1)));
                    b_P  =  -(b_Ex + b_Wx + b_Wy + b_Ez);
                    b_Ex(isnan(b_Ex)) = 0;
                    b_Ez(isnan(b_Ez)) = 0;
                    b_P(isnan(b_P)) = 0;
                    b_Wy(isnan(b_Wy)) = 0;
                    b_Wx(isnan(b_Wx)) = 0;
                B(n,[(n-max_y) (n-1) n (n+max_y) (n+(max_y*max_x))]) = [b_Wx b_Wy b_P b_Ex b_Ez];
                    
                elseif (i==max_y && j==max_x)
%                     B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) =zeros;
                    b_Wx = -L_ek(i,j,k)*L_ek(i,j-1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(i,j,k)*0.5*(x(j)-x(j-1))+L_ek(i,j-1,k)*0.5*(x(j+1)-x(j)));
                    b_Wy = -L_ek(i,j,k)*L_ek(i-1,j,k)*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(i,j,k)*0.5*(y(i)-y(i-1))+L_ek(i-1,j,k)*0.5*(y(i+1)-y(i)));
                    b_Ez = -L_ek(i,j,k)*L_ek(i,j,k+1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(L_ek(i,j,k+1)*0.5*(z(k+1)-z(k))+L_ek(i,j,k)*0.5*(z(k+2)-z(k+1)));
                    b_Wz = -L_ek(i,j,k)*L_ek(i,j,k-1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(L_ek(i,j,k)*0.5*(z(k)-z(k-1))+L_ek(i,j,k-1)*0.5*(z(k+1)-z(k)));
                    b_P  =  -(b_Wx + b_Wy + b_Ez + b_Wz);
                    b_Ez(isnan(b_Ez)) = 0;
                    b_P(isnan(b_P)) = 0;
                    b_Wz(isnan(b_Wz)) = 0;
                    b_Wy(isnan(b_Wy)) = 0;
                    b_Wx(isnan(b_Wx)) = 0;
                B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+(max_y*max_x))]) = [b_Wz b_Wx b_Wy b_P b_Ez];
                    
                elseif (i==max_y && j==1)
%                     B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) =zeros;
                    b_Ex = -L_ek(i,j,k)*L_ek(i,j+1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(i,j+1,k)*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1)));
                    b_Wy = -L_ek(i,j,k)*L_ek(i-1,j,k)*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(i,j,k)*0.5*(y(i)-y(i-1))+L_ek(i-1,j,k)*0.5*(y(i+1)-y(i)));
                    b_Ez = -L_ek(i,j,k)*L_ek(i,j,k+1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(L_ek(i,j,k+1)*0.5*(z(k+1)-z(k))+L_ek(i,j,k)*0.5*(z(k+2)-z(k+1)));
                    b_Wz = -L_ek(i,j,k)*L_ek(i,j,k-1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(L_ek(i,j,k)*0.5*(z(k)-z(k-1))+L_ek(i,j,k-1)*0.5*(z(k+1)-z(k)));
                    b_P  =  -(b_Ex + b_Wy + b_Ez + b_Wz);
                    b_Ex(isnan(b_Ex)) = 0;
                    b_Ez(isnan(b_Ez)) = 0;
                    b_P(isnan(b_P)) = 0;
                    b_Wz(isnan(b_Wz)) = 0;
                    b_Wy(isnan(b_Wy)) = 0;
                B(n,[(n-(max_y*max_x)) (n-1) n (n+max_y) (n+(max_y*max_x))]) = [b_Wz b_Wy b_P b_Ex b_Ez];
                    
                elseif (i==max_y && k==max_z)
%                     B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) =zeros;
                    b_Ex = -L_ek(i,j,k)*L_ek(i,j+1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(i,j+1,k)*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1)));
                    b_Wx = -L_ek(i,j,k)*L_ek(i,j-1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(i,j,k)*0.5*(x(j)-x(j-1))+L_ek(i,j-1,k)*0.5*(x(j+1)-x(j)));
                    b_Wy = -L_ek(i,j,k)*L_ek(i-1,j,k)*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(i,j,k)*0.5*(y(i)-y(i-1))+L_ek(i-1,j,k)*0.5*(y(i+1)-y(i)));
                    b_Wz = -L_ek(i,j,k)*L_ek(i,j,k-1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(L_ek(i,j,k)*0.5*(z(k)-z(k-1))+L_ek(i,j,k-1)*0.5*(z(k+1)-z(k)));
                    b_P  =  -(b_Ex + b_Wx + b_Wy + b_Wz);
                    b_Ex(isnan(b_Ex)) = 0;
                    b_P(isnan(b_P)) = 0;
                    b_Wz(isnan(b_Wz)) = 0;
                    b_Wy(isnan(b_Wy)) = 0;
                    b_Wx(isnan(b_Wx)) = 0;
                B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+max_y)]) = [b_Wz b_Wx b_Wy b_P b_Ex];
                    
                    
                elseif (j==max_x && k==1)
%                     B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) =zeros;
                    b_Wx = -L_ek(i,j,k)*L_ek(i,j-1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(i,j,k)*0.5*(x(j)-x(j-1))+L_ek(i,j-1,k)*0.5*(x(j+1)-x(j)));
                    b_Ey1 = (-L_ek(i,j,k)*L_ek(a(1),b(1),c(1))*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(a(1),b(1),c(1))*0.5*(y(i+1)-y(i))+L_ek(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    if c(2) == 0
                        b_Ey2 = 0;
                    else
                        b_Ey2 = (-L_ek(i,j,k)*L_ek(a(2),b(2),c(2))*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(a(2),b(2),c(2))*0.5*(y(i+1)-y(i))+L_ek(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    if c(3) == 0 
                        b_Ey3 = 0;
                    else
                    b_Ey3 = (-L_ek(i,j,k)*L_ek(a(3),b(3),c(3))*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(a(3),b(3),c(3))*0.5*(y(i+1)-y(i))+L_ek(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    b_Wy = -L_ek(i,j,k)*L_ek(i+1,j,k)*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(i,j,k)*0.5*(y(i)-y(i-1))+L_ek(i-1,j,k)*0.5*(y(i+1)-y(i)));
                    b_Ez = -L_ek(i,j,k)*L_ek(i,j,k+1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(L_ek(i,j,k+1)*0.5*(z(k+1)-z(k))+L_ek(i,j,k)*0.5*(z(k+2)-z(k+1)));
                    b_P  =  -(b_Wx + b_Ey1 +b_Ey2 +b_Ey3 + b_Wy + b_Ez);
                    b_Ey1(isnan(b_Ey1)) = 0;
                    b_Ey2(isnan(b_Ey2)) = 0;
                    b_Ey3(isnan(b_Ey3)) = 0;
                    b_Ez(isnan(b_Ez)) = 0;
                    b_P(isnan(b_P)) = 0;
                    b_Wy(isnan(b_Wy)) = 0;
                    b_Wx(isnan(b_Wx)) = 0;
                B(n,[(n-max_y) (n-1) n (nEy1) (nEy2) (nEy3) (n+(max_y*max_x))]) = [b_Wx b_Wy b_P b_Ey1 b_Ey2 b_Ey3 b_Ez];
                    
                elseif (j==max_x && k==max_z)
%                     B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) =zeros;
                    b_Wx = -L_ek(i,j,k)*L_ek(i,j-1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(i,j,k)*0.5*(x(j)-x(j-1))+L_ek(i,j-1,k)*0.5*(x(j+1)-x(j)));
                    b_Ey1 = (-L_ek(i,j,k)*L_ek(a(1),b(1),c(1))*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(a(1),b(1),c(1))*0.5*(y(i+1)-y(i))+L_ek(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    if c(2) == 0
                        b_Ey2 = 0;
                    else
                        b_Ey2 = (-L_ek(i,j,k)*L_ek(a(2),b(2),c(2))*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(a(2),b(2),c(2))*0.5*(y(i+1)-y(i))+L_ek(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    if c(3) == 0 
                        b_Ey3 = 0;
                    else
                    b_Ey3 = (-L_ek(i,j,k)*L_ek(a(3),b(3),c(3))*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(a(3),b(3),c(3))*0.5*(y(i+1)-y(i))+L_ek(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    b_Wy = -L_ek(i,j,k)*L_ek(i+1,j,k)*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(i,j,k)*0.5*(y(i)-y(i-1))+L_ek(i-1,j,k)*0.5*(y(i+1)-y(i)));
                    b_Wz = -L_ek(i,j,k)*L_ek(i,j,k-1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(L_ek(i,j,k)*0.5*(z(k)-z(k-1))+L_ek(i,j,k-1)*0.5*(z(k+1)-z(k)));
                    b_P  =  -(b_Wx + b_Ey1 +b_Ey2 +b_Ey3 + b_Wy + b_Wz);
                    b_Ey1(isnan(b_Ey1)) = 0;
                    b_Ey2(isnan(b_Ey2)) = 0;
                    b_Ey3(isnan(b_Ey3)) = 0;
                    b_P(isnan(b_P)) = 0;
                    b_Wz(isnan(b_Wz)) = 0;
                    b_Wy(isnan(b_Wy)) = 0;
                    b_Wx(isnan(b_Wx)) = 0;
                B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (nEy1) (nEy2) (nEy3)]) = [b_Wz b_Wx b_Wy b_P b_Ey1 b_Ey2 b_Ey3];
                    
                elseif (j==1 && k==max_z)
%                     B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) =zeros;
                    b_Ex = -L_ek(i,j,k)*L_ek(i,j+1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(i,j+1,k)*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1)));
                    b_Ey1 = (-L_ek(i,j,k)*L_ek(a(1),b(1),c(1))*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(a(1),b(1),c(1))*0.5*(y(i+1)-y(i))+L_ek(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    if c(2) == 0
                        b_Ey2 = 0;
                    else
                        b_Ey2 = (-L_ek(i,j,k)*L_ek(a(2),b(2),c(2))*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(a(2),b(2),c(2))*0.5*(y(i+1)-y(i))+L_ek(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    if c(3) == 0 
                        b_Ey3 = 0;
                    else
                    b_Ey3 = (-L_ek(i,j,k)*L_ek(a(3),b(3),c(3))*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(a(3),b(3),c(3))*0.5*(y(i+1)-y(i))+L_ek(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    b_Wy = -L_ek(i,j,k)*L_ek(i+1,j,k)*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(i,j,k)*0.5*(y(i)-y(i-1))+L_ek(i-1,j,k)*0.5*(y(i+1)-y(i)));
                    b_Wz = -L_ek(i,j,k)*L_ek(i,j,k-1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(L_ek(i,j,k)*0.5*(z(k)-z(k-1))+L_ek(i,j,k-1)*0.5*(z(k+1)-z(k)));
                    b_P  =  -(b_Ex + b_Ey1 +b_Ey2 +b_Ey3 + b_Wy + b_Wz);
                    b_Ex(isnan(b_Ex)) = 0;
                    b_Ey1(isnan(b_Ey1)) = 0;
                    b_Ey2(isnan(b_Ey2)) = 0;
                    b_Ey3(isnan(b_Ey3)) = 0;
                    b_P(isnan(b_P)) = 0;
                    b_Wz(isnan(b_Wz)) = 0;
                    b_Wy(isnan(b_Wy)) = 0;
                B(n,[(n-(max_y*max_x)) (n-1) n (nEy1) (nEy2) (nEy3) (n+max_y)]) = [b_Wz b_Wy b_P b_Ey1 b_Ey2 b_Ey3 b_Ex];
                    
                elseif (j==1 && k==1)
%                     B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) =zeros;
                    b_Ex = -L_ek(i,j,k)*L_ek(i,j+1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(i,j+1,k)*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1)));
                    b_Ey1 = (-L_ek(i,j,k)*L_ek(a(1),b(1),c(1))*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(a(1),b(1),c(1))*0.5*(y(i+1)-y(i))+L_ek(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    if c(2) == 0
                        b_Ey2 = 0;
                    else
                        b_Ey2 = (-L_ek(i,j,k)*L_ek(a(2),b(2),c(2))*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(a(2),b(2),c(2))*0.5*(y(i+1)-y(i))+L_ek(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    if c(3) == 0 
                        b_Ey3 = 0;
                    else
                    b_Ey3 = (-L_ek(i,j,k)*L_ek(a(3),b(3),c(3))*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(a(3),b(3),c(3))*0.5*(y(i+1)-y(i))+L_ek(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    b_Wy = -L_ek(i,j,k)*L_ek(i+1,j,k)*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(i,j,k)*0.5*(y(i)-y(i-1))+L_ek(i-1,j,k)*0.5*(y(i+1)-y(i)));
                    b_Ez = -L_ek(i,j,k)*L_ek(i,j,k+1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(L_ek(i,j,k+1)*0.5*(z(k+1)-z(k))+L_ek(i,j,k)*0.5*(z(k+2)-z(k+1)));
                    b_P  =  -(b_Ex + b_Ey1 +b_Ey2 +b_Ey3 + b_Wy + b_Ez);
                    b_Ex(isnan(b_Ex)) = 0;
                    b_Ey1(isnan(b_Ey1)) = 0;
                    b_Ey2(isnan(b_Ey2)) = 0;
                    b_Ey3(isnan(b_Ey3)) = 0;
                    b_Ez(isnan(b_Ez)) = 0;
                    b_P(isnan(b_P)) = 0;
                    b_Wy(isnan(b_Wy)) = 0;
                B(n,[(n-1) n (nEy1) (nEy2) (nEy3) (n+max_y) (n+(max_y*max_x))]) = [b_Wy b_P b_Ey1 b_Ey2 b_Ey3 b_Ex b_Ez];
                    
                    %----------------this part defines the surfaces
                    %defining the surface boundaries---------------
                elseif i==1
%                     B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) =zeros;
                    b_Ex = -L_ek(i,j,k)*L_ek(i,j+1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(i,j+1,k)*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1)));
                    b_Wx = -L_ek(i,j,k)*L_ek(i,j-1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(i,j,k)*0.5*(x(j)-x(j-1))+L_ek(i,j-1,k)*0.5*(x(j+1)-x(j)));
                    b_Ey1 = (-L_ek(i,j,k)*L_ek(a(1),b(1),c(1))*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(a(1),b(1),c(1))*0.5*(y(i+1)-y(i))+L_ek(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    if c(2) == 0
                        b_Ey2 = 0;
                    else
                        b_Ey2 = (-L_ek(i,j,k)*L_ek(a(2),b(2),c(2))*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(a(2),b(2),c(2))*0.5*(y(i+1)-y(i))+L_ek(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    if c(3) == 0 
                        b_Ey3 = 0;
                    else
                    b_Ey3 = (-L_ek(i,j,k)*L_ek(a(3),b(3),c(3))*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(a(3),b(3),c(3))*0.5*(y(i+1)-y(i))+L_ek(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    b_Ez = -L_ek(i,j,k)*L_ek(i,j,k+1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(L_ek(i,j,k+1)*0.5*(z(k+1)-z(k))+L_ek(i,j,k)*0.5*(z(k+2)-z(k+1)));
                    b_Wz = -L_ek(i,j,k)*L_ek(i,j,k-1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(L_ek(i,j,k)*0.5*(z(k)-z(k-1))+L_ek(i,j,k-1)*0.5*(z(k+1)-z(k)));
                    b_P  =  -(b_Ex + b_Wx + b_Ey1 +b_Ey2 +b_Ey3 + b_Ez + b_Wz);
                    b_Ex(isnan(b_Ex)) = 0;
                    b_Ey1(isnan(b_Ey1)) = 0;
                    b_Ey2(isnan(b_Ey2)) = 0;
                    b_Ey3(isnan(b_Ey3)) = 0;
                    b_Ez(isnan(b_Ez)) = 0;
                    b_P(isnan(b_P)) = 0;
                    b_Wz(isnan(b_Wz)) = 0;
                    b_Wx(isnan(b_Wx)) = 0;
                B(n,[(n-(max_y*max_x)) (n-max_y) n (nEy1) (nEy2) (nEy3) (n+max_y) (n+(max_y*max_x))]) = [b_Wz b_Wx b_P b_Ey1 b_Ey2 b_Ey3 b_Ex b_Ez];
                    
                elseif i==max_y
%                     B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) =zeros;
                    b_Ex = -L_ek(i,j,k)*L_ek(i,j+1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(i,j+1,k)*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1)));
                    b_Wx = -L_ek(i,j,k)*L_ek(i,j-1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(i,j,k)*0.5*(x(j)-x(j-1))+L_ek(i,j-1,k)*0.5*(x(j+1)-x(j)));
                    b_Wy = -L_ek(i,j,k)*L_ek(i-1,j,k)*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(i,j,k)*0.5*(y(i)-y(i-1))+L_ek(i-1,j,k)*0.5*(y(i+1)-y(i)));
                    b_Ez = -L_ek(i,j,k)*L_ek(i,j,k+1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(L_ek(i,j,k+1)*0.5*(z(k+1)-z(k))+L_ek(i,j,k)*0.5*(z(k+2)-z(k+1)));
                    b_Wz = -L_ek(i,j,k)*L_ek(i,j,k-1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(L_ek(i,j,k)*0.5*(z(k)-z(k-1))+L_ek(i,j,k-1)*0.5*(z(k+1)-z(k)));
                    b_P  =  -(b_Ex + b_Wx + b_Wy + b_Ez + b_Wz);
                    b_Ex(isnan(b_Ex)) = 0;
                    b_Ez(isnan(b_Ez)) = 0;
                    b_P(isnan(b_P)) = 0;
                    b_Wz(isnan(b_Wz)) = 0;
                    b_Wy(isnan(b_Wy)) = 0;
                    b_Wx(isnan(b_Wx)) = 0;
                B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+max_y) (n+(max_y*max_x))]) = [b_Wz b_Wx b_Wy b_P b_Ex b_Ez];
                    
                elseif j==1
%                     B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) =zeros;
                    b_Ex = -L_ek(i,j,k)*L_ek(i,j+1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(i,j+1,k)*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1)));
                    b_Ey1 = (-L_ek(i,j,k)*L_ek(a(1),b(1),c(1))*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(a(1),b(1),c(1))*0.5*(y(i+1)-y(i))+L_ek(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    if c(2) == 0
                        b_Ey2 = 0;
                    else
                        b_Ey2 = (-L_ek(i,j,k)*L_ek(a(2),b(2),c(2))*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(a(2),b(2),c(2))*0.5*(y(i+1)-y(i))+L_ek(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    if c(3) == 0 
                        b_Ey3 = 0;
                    else
                    b_Ey3 = (-L_ek(i,j,k)*L_ek(a(3),b(3),c(3))*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(a(3),b(3),c(3))*0.5*(y(i+1)-y(i))+L_ek(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    b_Wy = -L_ek(i,j,k)*L_ek(i-1,j,k)*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(i,j,k)*0.5*(y(i)-y(i-1))+L_ek(i-1,j,k)*0.5*(y(i+1)-y(i)));
                    b_Ez = -L_ek(i,j,k)*L_ek(i,j,k+1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(L_ek(i,j,k+1)*0.5*(z(k+1)-z(k))+L_ek(i,j,k)*0.5*(z(k+2)-z(k+1)));
                    b_Wz = -L_ek(i,j,k)*L_ek(i,j,k-1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(L_ek(i,j,k)*0.5*(z(k)-z(k-1))+L_ek(i,j,k-1)*0.5*(z(k+1)-z(k)));
                    b_P  =  -(b_Ex + b_Ey1 +b_Ey2 +b_Ey3 + b_Wy + b_Ez + b_Wz);
                    b_Ex(isnan(b_Ex)) = 0;
                    b_Ey1(isnan(b_Ey1)) = 0;
                    b_Ey2(isnan(b_Ey2)) = 0;
                    b_Ey3(isnan(b_Ey3)) = 0;
                    b_Ez(isnan(b_Ez)) = 0;
                    b_P(isnan(b_P)) = 0;
                    b_Wz(isnan(b_Wz)) = 0;
                    b_Wy(isnan(b_Wy)) = 0;
                B(n,[(n-(max_y*max_x)) (n-1) n (nEy1) (nEy2) (nEy3) (n+max_y) (n+(max_y*max_x))]) = [b_Wz b_Wy b_P b_Ey1 b_Ey2 b_Ey3 b_Ex b_Ez];
                    
                elseif j==max_x
%                     B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) =zeros;
                    b_Wx = -L_ek(i,j,k)*L_ek(i,j-1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(i,j,k)*0.5*(x(j)-x(j-1))+L_ek(i,j-1,k)*0.5*(x(j+1)-x(j)));
                    b_Ey1 = (-L_ek(i,j,k)*L_ek(a(1),b(1),c(1))*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(a(1),b(1),c(1))*0.5*(y(i+1)-y(i))+L_ek(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    if c(2) == 0
                        b_Ey2 = 0;
                    else
                        b_Ey2 = (-L_ek(i,j,k)*L_ek(a(2),b(2),c(2))*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(a(2),b(2),c(2))*0.5*(y(i+1)-y(i))+L_ek(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    if c(3) == 0 
                        b_Ey3 = 0;
                    else
                    b_Ey3 = (-L_ek(i,j,k)*L_ek(a(3),b(3),c(3))*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(a(3),b(3),c(3))*0.5*(y(i+1)-y(i))+L_ek(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    b_Wy = -L_ek(i,j,k)*L_ek(i-1,j,k)*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(i,j,k)*0.5*(y(i)-y(i-1))+L_ek(i-1,j,k)*0.5*(y(i+1)-y(i)));
                    b_Ez = -L_ek(i,j,k)*L_ek(i,j,k+1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(L_ek(i,j,k+1)*0.5*(z(k+1)-z(k))+L_ek(i,j,k)*0.5*(z(k+2)-z(k+1)));
                    b_Wz = -L_ek(i,j,k)*L_ek(i,j,k-1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(L_ek(i,j,k)*0.5*(z(k)-z(k-1))+L_ek(i,j,k-1)*0.5*(z(k+1)-z(k)));
                    b_P  =  -(b_Wx + b_Ey1 +b_Ey2 +b_Ey3 + b_Wy + b_Ez + b_Wz);
                    b_Ey1(isnan(b_Ey1)) = 0;
                    b_Ey2(isnan(b_Ey2)) = 0;
                    b_Ey3(isnan(b_Ey3)) = 0;
                    b_Ez(isnan(b_Ez)) = 0;
                    b_P(isnan(b_P)) = 0;
                    b_Wz(isnan(b_Wz)) = 0;
                    b_Wy(isnan(b_Wy)) = 0;
                    b_Wx(isnan(b_Wx)) = 0;
                B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (nEy1) (nEy2) (nEy3) (n+(max_y*max_x))]) = [b_Wz b_Wx b_Wy b_P b_Ey1 b_Ey2 b_Ey3 b_Ez];
                    
                elseif k==1
%                     B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) =zeros;
                    b_Ex = -L_ek(i,j,k)*L_ek(i,j+1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(i,j+1,k)*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1)));
                    b_Wx = -L_ek(i,j,k)*L_ek(i,j-1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(i,j,k)*0.5*(x(j)-x(j-1))+L_ek(i,j-1,k)*0.5*(x(j+1)-x(j)));
                    b_Ey1 = (-L_ek(i,j,k)*L_ek(a(1),b(1),c(1))*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(a(1),b(1),c(1))*0.5*(y(i+1)-y(i))+L_ek(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    if c(2) == 0
                        b_Ey2 = 0;
                    else
                        b_Ey2 = (-L_ek(i,j,k)*L_ek(a(2),b(2),c(2))*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(a(2),b(2),c(2))*0.5*(y(i+1)-y(i))+L_ek(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    if c(3) == 0 
                        b_Ey3 = 0;
                    else
                    b_Ey3 = (-L_ek(i,j,k)*L_ek(a(3),b(3),c(3))*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(a(3),b(3),c(3))*0.5*(y(i+1)-y(i))+L_ek(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    b_Wy = -L_ek(i,j,k)*L_ek(i-1,j,k)*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(i,j,k)*0.5*(y(i)-y(i-1))+L_ek(i-1,j,k)*0.5*(y(i+1)-y(i)));
                    b_Ez = -L_ek(i,j,k)*L_ek(i,j,k+1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(L_ek(i,j,k+1)*0.5*(z(k+1)-z(k))+L_ek(i,j,k)*0.5*(z(k+2)-z(k+1)));
                    b_P  =  -(b_Ex + b_Wx + b_Ey1 +b_Ey2 +b_Ey3 + b_Wy + b_Ez);
                    b_Ex(isnan(b_Ex)) = 0;
                    b_Ey1(isnan(b_Ey1)) = 0;
                    b_Ey2(isnan(b_Ey2)) = 0;
                    b_Ey3(isnan(b_Ey3)) = 0;
                    b_Ez(isnan(b_Ez)) = 0;
                    b_P(isnan(b_P)) = 0;
                    b_Wy(isnan(b_Wy)) = 0;
                    b_Wx(isnan(b_Wx)) = 0;
                B(n,[(n-max_y) (n-1) n (nEy1) (nEy2) (nEy3) (n+max_y) (n+(max_y*max_x))]) = [b_Wx b_Wy b_P b_Ey1 b_Ey2 b_Ey3 b_Ex b_Ez];
                
                elseif k==max_z
%                     B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) =zeros;
                    b_Ex = -L_ek(i,j,k)*L_ek(i,j+1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(i,j+1,k)*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1)));
                    b_Wx = -L_ek(i,j,k)*L_ek(i,j-1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(i,j,k)*0.5*(x(j)-x(j-1))+L_ek(i,j-1,k)*0.5*(x(j+1)-x(j)));
                    b_Ey1 = (-L_ek(i,j,k)*L_ek(a(1),b(1),c(1))*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(a(1),b(1),c(1))*0.5*(y(i+1)-y(i))+L_ek(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    if c(2) == 0
                        b_Ey2 = 0;
                    else
                        b_Ey2 = (-L_ek(i,j,k)*L_ek(a(2),b(2),c(2))*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(a(2),b(2),c(2))*0.5*(y(i+1)-y(i))+L_ek(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    if c(3) == 0 
                        b_Ey3 = 0;
                    else
                    b_Ey3 = (-L_ek(i,j,k)*L_ek(a(3),b(3),c(3))*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(a(3),b(3),c(3))*0.5*(y(i+1)-y(i))+L_ek(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    b_Wy = -L_ek(i,j,k)*L_ek(i-1,j,k)*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(i,j,k)*0.5*(y(i)-y(i-1))+L_ek(i-1,j,k)*0.5*(y(i+1)-y(i)));
                    b_Wz = -L_ek(i,j,k)*L_ek(i,j,k-1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(L_ek(i,j,k)*0.5*(z(k)-z(k-1))+L_ek(i,j,k-1)*0.5*(z(k+1)-z(k)));
                    b_P  =  -(b_Ex + b_Wx + b_Ey1 +b_Ey2 +b_Ey3 + b_Wy + b_Wz);
                    b_Ex(isnan(b_Ex)) = 0;
                    b_Ey1(isnan(b_Ey1)) = 0;
                    b_Ey2(isnan(b_Ey2)) = 0;
                    b_Ey3(isnan(b_Ey3)) = 0;
                    b_P(isnan(b_P)) = 0;
                    b_Wz(isnan(b_Wz)) = 0;
                    b_Wy(isnan(b_Wy)) = 0;
                    b_Wx(isnan(b_Wx)) = 0;
                B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (nEy1) (nEy2) (nEy3) (n+max_y)]) = [b_Wz b_Wx b_Wy b_P b_Ey1 b_Ey2 b_Ey3 b_Ex];
                
               elseif (i>1 && j>1 && k>1) && (i<max_y && j<max_x && k<max_z)
%                    B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) =zeros;
                    b_Ex = -wsize(1)*L_ek(i,j,k)*L_ek(i,j+1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(i,j+1,k)*0.5*(x(j+1)-x(j))+L_ek(i,j,k)*0.5*(x(j+2)-x(j+1)));
                    b_Wx = -wsize(1)*L_ek(i,j,k)*L_ek(i,j-1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(L_ek(i,j,k)*0.5*(x(j)-x(j-1))+L_ek(i,j-1,k)*0.5*(x(j+1)-x(j)));
                    b_Ey1 = (-L_ek(i,j,k)*L_ek(a(1),b(1),c(1))*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(a(1),b(1),c(1))*0.5*(y(i+1)-y(i))+L_ek(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    if c(2) == 0
                        b_Ey2 = 0;
                    else
                        b_Ey2 = (-L_ek(i,j,k)*L_ek(a(2),b(2),c(2))*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(a(2),b(2),c(2))*0.5*(y(i+1)-y(i))+L_ek(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    if c(3) == 0 
                        b_Ey3 = 0;
                    else
                    b_Ey3 = (-L_ek(i,j,k)*L_ek(a(3),b(3),c(3))*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(a(3),b(3),c(3))*0.5*(y(i+1)-y(i))+L_ek(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    b_Wy = -wsize(1)*L_ek(i,j,k)*L_ek(i-1,j,k)*(x(j+1)-x(j))*(z(k+1)-z(k))/(L_ek(i,j,k)*0.5*(y(i)-y(i-1))+L_ek(i-1,j,k)*0.5*(y(i+1)-y(i)));
                    b_Ez = -wsize(1)*L_ek(i,j,k)*L_ek(i,j,k+1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(L_ek(i,j,k+1)*0.5*(z(k+1)-z(k))+L_ek(i,j,k)*0.5*(z(k+2)-z(k+1)));
                    b_Wz = -wsize(1)*L_ek(i,j,k)*L_ek(i,j,k-1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(L_ek(i,j,k)*0.5*(z(k)-z(k-1))+L_ek(i,j,k-1)*0.5*(z(k+1)-z(k)));
                    b_P  =  -wsize(1)*(b_Ex + b_Wx + b_Ey1 +b_Ey2 +b_Ey3 + b_Wy + b_Ez + b_Wz);
                    b_Ex(isnan(b_Ex)) = 0;
                    b_Ey(isnan(b_Ey)) = 0;
                    b_Ez(isnan(b_Ez)) = 0;
                    b_P(isnan(b_P)) = 0;
                    b_Wz(isnan(b_Wz)) = 0;
                    b_Wy(isnan(b_Wy)) = 0;
                    b_Wx(isnan(b_Wx)) = 0;
                    
%                 B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) = [b_Wz b_Wx b_Wy b_P b_Ey b_Ex b_Ez];
                B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (nEy1) (nEy2) (nEy3) (n+max_y) (n+(max_y*max_x))]) = [b_Wz b_Wx b_Wy b_P b_Ey1 b_Ey2 b_Ey3 b_Ex b_Ez];
                    
                end           
    
    
end
