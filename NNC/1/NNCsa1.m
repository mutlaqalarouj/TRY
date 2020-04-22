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
                    a_Ex = -sigma(i,j,k)*sigma(i,j+1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(sigma(i,j+1,k)*0.5*(x(j+1)-x(j))+sigma(i,j,k)*0.5*(x(j+2)-x(j+1)));
                    a_Ey1 = (-sigma(i,j,k)*sigma(a(1),b(1),c(1))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(1),b(1),c(1))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    if c(2) == 0
                        a_Ey2 = 0;
                    else
                        a_Ey2 = (-sigma(i,j,k)*sigma(a(2),b(2),c(2))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(2),b(2),c(2))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    if c(3) == 0 
                        a_Ey3 = 0;
                    else
                    a_Ey3 = (-sigma(i,j,k)*sigma(a(3),b(3),c(3))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(3),b(3),c(3))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    a_Ez = -sigma(i,j,k)*sigma(i,j,k+1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(sigma(i,j,k+1)*0.5*(z(k+1)-z(k))+sigma(i,j,k)*0.5*(z(k+2)-z(k+1)));
                    a_P  =  -(a_Ex + (a_Ey1+a_Ey2+a_Ey3) + a_Ez);
                    a_Ex(isnan(a_Ex)) = 0;
                    a_Ey1(isnan(a_Ey1)) = 0;
                    a_Ey2(isnan(a_Ey2)) = 0;
                    a_Ey3(isnan(a_Ey3)) = 0;
                    a_Ez(isnan(a_Ez)) = 0;
                    a_P(isnan(a_P)) = 0;
                A(n,[n (nEy1) (nEy2) (nEy2) (n+max_y) (n+(max_y*max_x))]) = [a_P a_Ey1 a_Ey2 a_Ey3 a_Ex a_Ez];
                    
                elseif (j==max_x && i==max_y) && k==max_z
%                   B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) =zeros;  
                    a_Wx = -sigma(i,j,k)*sigma(i,j-1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(sigma(i,j,k)*0.5*(x(j)-x(j-1))+sigma(i,j-1,k)*0.5*(x(j+1)-x(j)));
                    a_Wy = -sigma(i,j,k)*sigma(i-1,j,k)*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(i,j,k)*0.5*(y(i)-y(i-1))+sigma(i-1,j,k)*0.5*(y(i+1)-y(i)));
                    a_Wz = -sigma(i,j,k)*sigma(i,j,k-1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(sigma(i,j,k)*0.5*(z(k)-z(k-1))+sigma(i,j,k-1)*0.5*(z(k+1)-z(k)));
                    a_P  =  -(a_Wx + a_Wy + a_Wz);
                    a_P(isnan(a_P)) = 0;
                    a_Wz(isnan(a_Wz)) = 0;
                    a_Wy(isnan(a_Wy)) = 0;
                    a_Wx(isnan(a_Wx)) = 0;
                A(n,[n-(max_y*max_x) n-max_y n-1 n]) = [a_Wz a_Wx a_Wy a_P];     
                    
                elseif (j==1 && i==1) && k==max_z
%                   B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) =zeros; 
                    a_Ex = -sigma(i,j,k)*sigma(i,j+1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(sigma(i,j+1,k)*0.5*(x(j+1)-x(j))+sigma(i,j,k)*0.5*(x(j+2)-x(j+1)));
                    a_Ey1 = (-sigma(i,j,k)*sigma(a(1),b(1),c(1))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(1),b(1),c(1))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    if c(2) == 0
                        a_Ey2 = 0;
                    else
                        a_Ey2 = (-sigma(i,j,k)*sigma(a(2),b(2),c(2))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(2),b(2),c(2))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    if c(3) == 0 
                        a_Ey3 = 0;
                    else
                    a_Ey3 = (-sigma(i,j,k)*sigma(a(3),b(3),c(3))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(3),b(3),c(3))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    a_Wz = -sigma(i,j,k)*sigma(i,j,k-1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(sigma(i,j,k)*0.5*(z(k)-z(k-1))+sigma(i,j,k-1)*0.5*(z(k+1)-z(k)));
                    a_P  =  -(a_Ex + a_Ey1 +a_Ey2 +a_Ey3 + a_Wz);
                    a_Ex(isnan(a_Ex)) = 0;
                    a_Ey1(isnan(a_Ey1)) = 0;
                    a_Ey2(isnan(a_Ey2)) = 0;
                    a_Ey3(isnan(a_Ey3)) = 0;
                    a_P(isnan(a_P)) = 0;
                    a_Wz(isnan(a_Wz)) = 0;
                A(n,[(n-(max_y*max_x)) n (nEy1) (nEy2) (nEy2) (n+max_y)]) = [a_Wz a_P a_Ey1 a_Ey2 a_Ey3 a_Ex]; 
                
                elseif (j==max_x && i==1) && k==max_z
%                    B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) =zeros;
                    a_Wx = -sigma(i,j,k)*sigma(i,j-1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(sigma(i,j,k)*0.5*(x(j)-x(j-1))+sigma(i,j-1,k)*0.5*(x(j+1)-x(j)));
                    a_Ey1 = (-sigma(i,j,k)*sigma(a(1),b(1),c(1))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(1),b(1),c(1))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    if c(2) == 0
                        a_Ey2 = 0;
                    else
                        a_Ey2 = (-sigma(i,j,k)*sigma(a(2),b(2),c(2))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(2),b(2),c(2))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    if c(3) == 0 
                        a_Ey3 = 0;
                    else
                    a_Ey3 = (-sigma(i,j,k)*sigma(a(3),b(3),c(3))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(3),b(3),c(3))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    a_Wz = -sigma(i,j,k)*sigma(i,j,k-1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(sigma(i,j,k)*0.5*(z(k)-z(k-1))+sigma(i,j,k-1)*0.5*(z(k+1)-z(k)));
                    a_P  =  -(a_Wx + a_Ey1 +a_Ey2 +a_Ey3 + a_Wz);
                    a_Ey1(isnan(a_Ey1)) = 0;
                    a_Ey2(isnan(a_Ey2)) = 0;
                    a_Ey3(isnan(a_Ey3)) = 0;
                    a_P(isnan(a_P)) = 0;
                    a_Wz(isnan(a_Wz)) = 0;
                    a_Wx(isnan(a_Wx)) = 0;
                A(n,[(n-(max_y*max_x)) (n-max_y) n (nEy1) (nEy2) (nEy2)]) = [a_Wz a_Wx a_P a_Ey1 a_Ey2 a_Ey3];
                    
                elseif (j==max_x && i==1) && k==1
%                     B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) =zeros;
                    a_Wx = -sigma(i,j,k)*sigma(i,j-1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(sigma(i,j,k)*0.5*(x(j)-x(j-1))+sigma(i,j-1,k)*0.5*(x(j+1)-x(j)));
                    a_Ey1 = (-sigma(i,j,k)*sigma(a(1),b(1),c(1))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(1),b(1),c(1))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    if c(2) == 0
                        a_Ey2 = 0;
                    else
                        a_Ey2 = (-sigma(i,j,k)*sigma(a(2),b(2),c(2))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(2),b(2),c(2))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    if c(3) == 0 
                        a_Ey3 = 0;
                    else
                    a_Ey3 = (-sigma(i,j,k)*sigma(a(3),b(3),c(3))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(3),b(3),c(3))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    a_Ez = -sigma(i,j,k)*sigma(i,j,k+1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(sigma(i,j,k+1)*0.5*(z(k+1)-z(k))+sigma(i,j,k)*0.5*(z(k+2)-z(k+1)));
                    a_P  =  -(a_Wx + a_Ey1 +a_Ey2 +a_Ey3 + a_Ez);
                    a_Ey1(isnan(a_Ey1)) = 0;
                    a_Ey2(isnan(a_Ey2)) = 0;
                    a_Ey3(isnan(a_Ey3)) = 0;
                    a_Ez(isnan(a_Ez)) = 0;
                    a_P(isnan(a_P)) = 0;
                    a_Wx(isnan(a_Wx)) = 0;
                A(n,[(n-max_y) n (n+max_y) (n+(max_y*max_x))]) = [a_Wx a_P a_Ey1 a_Ey2 a_Ey3 a_Ez];
                    
                elseif (j==max_x && i==max_y) && k==1
%                     B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) =zeros;
                    a_Wx = -sigma(i,j,k)*sigma(i,j-1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(sigma(i,j,k)*0.5*(x(j)-x(j-1))+sigma(i,j-1,k)*0.5*(x(j+1)-x(j)));
                    a_Wy = -sigma(i,j,k)*sigma(i-1,j,k)*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(i,j,k)*0.5*(y(i)-y(i-1))+sigma(i-1,j,k)*0.5*(y(i+1)-y(i)));
                    a_Ez = -sigma(i,j,k)*sigma(i,j,k+1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(sigma(i,j,k+1)*0.5*(z(k+1)-z(k))+sigma(i,j,k)*0.5*(z(k+2)-z(k+1)));
                    a_P  =  -(a_Wx + a_Wy + a_Ez);
                    a_Ez(isnan(a_Ez)) = 0;
                    a_P(isnan(a_P)) = 0;
                    a_Wy(isnan(a_Wy)) = 0;
                    a_Wx(isnan(a_Wx)) = 0;
                A(n,[(n-max_y) (n-1) n (n+(max_y*max_x))]) = [a_Wx a_Wy a_P a_Ez];
                    
                elseif (j==1 && i==max_y) && k==max_z
%                     B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) =zeros;
                    a_Ex = -sigma(i,j,k)*sigma(i,j+1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(sigma(i,j+1,k)*0.5*(x(j+1)-x(j))+sigma(i,j,k)*0.5*(x(j+2)-x(j+1)));
                    a_Wy = -sigma(i,j,k)*sigma(i-1,j,k)*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(i,j,k)*0.5*(y(i)-y(i-1))+sigma(i-1,j,k)*0.5*(y(i+1)-y(i)));
                    a_Wz = -sigma(i,j,k)*sigma(i,j,k-1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(sigma(i,j,k)*0.5*(z(k)-z(k-1))+sigma(i,j,k-1)*0.5*(z(k+1)-z(k)));
                    a_P  =  -(a_Ex + a_Wy + a_Wz);
                    a_Ex(isnan(a_Ex)) = 0;
                    a_P(isnan(a_P)) = 0;
                    a_Wz(isnan(a_Wz)) = 0;
                    a_Wy(isnan(a_Wy)) = 0;
                    A(n,[(n-(max_y*max_x)) (n-1) n (n+max_y)]) = [a_Wz a_Wy a_P a_Ex];
                    
                elseif (j==1 && i==max_y) && k==1
%                     B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) =zeros;
                    a_Ex = -sigma(i,j,k)*sigma(i,j+1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(sigma(i,j+1,k)*0.5*(x(j+1)-x(j))+sigma(i,j,k)*0.5*(x(j+2)-x(j+1)));
                    a_Wy = -sigma(i,j,k)*sigma(i-1,j,k)*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(i,j,k)*0.5*(y(i)-y(i-1))+sigma(i-1,j,k)*0.5*(y(i+1)-y(i)));
                    a_Ez = -sigma(i,j,k)*sigma(i,j,k+1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(sigma(i,j,k+1)*0.5*(z(k+1)-z(k))+sigma(i,j,k)*0.5*(z(k+2)-z(k+1)));
                    a_P  =  -(a_Ex + a_Wy + a_Ez);
                    a_Ex(isnan(a_Ex)) = 0;
                    a_Ez(isnan(a_Ez)) = 0;
                    a_P(isnan(a_P)) = 0;
                    a_Wy(isnan(a_Wy)) = 0;
                A(n,[(n-1) n (n+max_y) (n+(max_y*max_x))]) = [a_Wy a_P a_Ex a_Ez];
                    
                    %----------------defines edges of resistance matrix -----------
                elseif (i==1 && k==1)
%                     B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) =zeros;
                    a_Ex = -sigma(i,j,k)*sigma(i,j+1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(sigma(i,j+1,k)*0.5*(x(j+1)-x(j))+sigma(i,j,k)*0.5*(x(j+2)-x(j+1)));
                    a_Wx = -sigma(i,j,k)*sigma(i,j-1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(sigma(i,j,k)*0.5*(x(j)-x(j-1))+sigma(i,j-1,k)*0.5*(x(j+1)-x(j)));
                    a_Ey1 = (-sigma(i,j,k)*sigma(a(1),b(1),c(1))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(1),b(1),c(1))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    if c(2) == 0
                        a_Ey2 = 0;
                    else
                        a_Ey2 = (-sigma(i,j,k)*sigma(a(2),b(2),c(2))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(2),b(2),c(2))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    if c(3) == 0 
                        a_Ey3 = 0;
                    else
                    a_Ey3 = (-sigma(i,j,k)*sigma(a(3),b(3),c(3))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(3),b(3),c(3))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    a_Ez = -sigma(i,j,k)*sigma(i,j,k+1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(sigma(i,j,k+1)*0.5*(z(k+1)-z(k))+sigma(i,j,k)*0.5*(z(k+2)-z(k+1)));
                    a_P  =  -(a_Ex + a_Wx + a_Ey1 +a_Ey2 +a_Ey3 + a_Ez);
                    a_Ex(isnan(a_Ex)) = 0;
                    a_Ey1(isnan(a_Ey1)) = 0;
                    a_Ey2(isnan(a_Ey2)) = 0;
                    a_Ey3(isnan(a_Ey3)) = 0;
                    a_Ez(isnan(a_Ez)) = 0;
                    a_P(isnan(a_P)) = 0;
                A(n,[(n-max_y) n (nEy1) (nEy2) (nEy2) (n+max_y) (n+(max_y*max_x))]) = [a_Wx a_P a_Ey1 a_Ey2 a_Ey3 a_Ex a_Ez];
                
                elseif (i==1 && j==1)
%                     B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) =zeros;
                    a_Ex = -sigma(i,j,k)*sigma(i,j+1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(sigma(i,j+1,k)*0.5*(x(j+1)-x(j))+sigma(i,j,k)*0.5*(x(j+2)-x(j+1)));
                    a_Ey1 = (-sigma(i,j,k)*sigma(a(1),b(1),c(1))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(1),b(1),c(1))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    if c(2) == 0
                        a_Ey2 = 0;
                    else
                        a_Ey2 = (-sigma(i,j,k)*sigma(a(2),b(2),c(2))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(2),b(2),c(2))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    if c(3) == 0 
                        a_Ey3 = 0;
                    else
                    a_Ey3 = (-sigma(i,j,k)*sigma(a(3),b(3),c(3))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(3),b(3),c(3))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    a_Ez = -sigma(i,j,k)*sigma(i,j,k+1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(sigma(i,j,k+1)*0.5*(z(k+1)-z(k))+sigma(i,j,k)*0.5*(z(k+2)-z(k+1)));
                    a_Wz = -sigma(i,j,k)*sigma(i,j,k-1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(sigma(i,j,k)*0.5*(z(k)-z(k-1))+sigma(i,j,k-1)*0.5*(z(k+1)-z(k)));
                    a_P  =  -(a_Ex + a_Ey1 +a_Ey2 +a_Ey3 + a_Ez + a_Wz);
                    a_Ex(isnan(a_Ex)) = 0;
                    a_Ey1(isnan(a_Ey1)) = 0;
                    a_Ey2(isnan(a_Ey2)) = 0;
                    a_Ey3(isnan(a_Ey3)) = 0;
                    a_Ez(isnan(a_Ez)) = 0;
                    a_P(isnan(a_P)) = 0;
                    a_Wz(isnan(a_Wz)) = 0;
                A(n,[(n-(max_y*max_x)) n (nEy1) (nEy2) (nEy2) (n+max_y) (n+(max_y*max_x))]) = [a_Wz a_P a_Ey1 a_Ey2 a_Ey3 a_Ex a_Ez];
                    
                elseif (i==1 && k==max_z)
%                     B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) =zeros;
                    a_Ex = -sigma(i,j,k)*sigma(i,j+1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(sigma(i,j+1,k)*0.5*(x(j+1)-x(j))+sigma(i,j,k)*0.5*(x(j+2)-x(j+1)));
                    a_Wx = -sigma(i,j,k)*sigma(i,j-1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(sigma(i,j,k)*0.5*(x(j)-x(j-1))+sigma(i,j-1,k)*0.5*(x(j+1)-x(j)));
                    a_Ey1 = (-sigma(i,j,k)*sigma(a(1),b(1),c(1))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(1),b(1),c(1))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    if c(2) == 0
                        a_Ey2 = 0;
                    else
                        a_Ey2 = (-sigma(i,j,k)*sigma(a(2),b(2),c(2))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(2),b(2),c(2))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    if c(3) == 0 
                        a_Ey3 = 0;
                    else
                    a_Ey3 = (-sigma(i,j,k)*sigma(a(3),b(3),c(3))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(3),b(3),c(3))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    a_Wz = -sigma(i,j,k)*sigma(i,j,k-1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(sigma(i,j,k)*0.5*(z(k)-z(k-1))+sigma(i,j,k-1)*0.5*(z(k+1)-z(k)));
                    a_P  =  -(a_Ex + a_Wx + a_Ey1 +a_Ey2 +a_Ey3 + a_Wz);
                    a_Ex(isnan(a_Ex)) = 0;
                    a_Ey1(isnan(a_Ey1)) = 0;
                    a_Ey2(isnan(a_Ey2)) = 0;
                    a_Ey3(isnan(a_Ey3)) = 0;
                    a_P(isnan(a_P)) = 0;
                    a_Wz(isnan(a_Wz)) = 0;
                    a_Wx(isnan(a_Wx)) = 0;
                A(n,[(n-(max_y*max_x)) (n-max_y) n (nEy1) (nEy2) (nEy2) (n+max_y)]) = [a_Wz a_Wx a_P a_Ey1 a_Ey2 a_Ey3 a_Ex];
                    
                elseif (i==1 && j==max_x)
%                     B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) =zeros;
                    a_Wx = -sigma(i,j,k)*sigma(i,j-1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(sigma(i,j,k)*0.5*(x(j)-x(j-1))+sigma(i,j-1,k)*0.5*(x(j+1)-x(j)));
                    a_Ey1 = (-sigma(i,j,k)*sigma(a(1),b(1),c(1))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(1),b(1),c(1))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    if c(2) == 0
                        a_Ey2 = 0;
                    else
                        a_Ey2 = (-sigma(i,j,k)*sigma(a(2),b(2),c(2))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(2),b(2),c(2))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    if c(3) == 0 
                        a_Ey3 = 0;
                    else
                    a_Ey3 = (-sigma(i,j,k)*sigma(a(3),b(3),c(3))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(3),b(3),c(3))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    a_Ez = -sigma(i,j,k)*sigma(i,j,k+1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(sigma(i,j,k+1)*0.5*(z(k+1)-z(k))+sigma(i,j,k)*0.5*(z(k+2)-z(k+1)));
                    a_Wz = -sigma(i,j,k)*sigma(i,j,k-1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(sigma(i,j,k)*0.5*(z(k)-z(k-1))+sigma(i,j,k-1)*0.5*(z(k+1)-z(k)));
                    a_P  =  -(a_Wx + a_Ey1 +a_Ey2 +a_Ey3 + a_Ez + a_Wz);
                    a_Ey1(isnan(a_Ey1)) = 0;
                    a_Ey2(isnan(a_Ey2)) = 0;
                    a_Ey3(isnan(a_Ey3)) = 0;
                    a_Ez(isnan(a_Ez)) = 0;
                    a_P(isnan(a_P)) = 0;
                    a_Wz(isnan(a_Wz)) = 0;
                    a_Wx(isnan(a_Wx)) = 0;
                A(n,[(n-(max_y*max_x)) (n-max_y) n (nEy1) (nEy2) (nEy2) (n+(max_y*max_x))]) = [a_Wz a_Wx a_P a_Ey1 a_Ey2 a_Ey3 a_Ez];
                   
                elseif (i==max_y && k==1)
%                     B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) =zeros;
                    a_Ex = -sigma(i,j,k)*sigma(i,j+1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(sigma(i,j+1,k)*0.5*(x(j+1)-x(j))+sigma(i,j,k)*0.5*(x(j+2)-x(j+1)));
                    a_Wx = -sigma(i,j,k)*sigma(i,j-1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(sigma(i,j,k)*0.5*(x(j)-x(j-1))+sigma(i,j-1,k)*0.5*(x(j+1)-x(j)));
                    a_Wy = -sigma(i,j,k)*sigma(i-1,j,k)*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(i,j,k)*0.5*(y(i)-y(i-1))+sigma(i-1,j,k)*0.5*(y(i+1)-y(i)));
                    a_Ez = -sigma(i,j,k)*sigma(i,j,k+1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(sigma(i,j,k+1)*0.5*(z(k+1)-z(k))+sigma(i,j,k)*0.5*(z(k+2)-z(k+1)));
                    a_P  =  -(a_Ex + a_Wx + a_Wy + a_Ez);
                    a_Ex(isnan(a_Ex)) = 0;
                    a_Ez(isnan(a_Ez)) = 0;
                    a_P(isnan(a_P)) = 0;
                    a_Wy(isnan(a_Wy)) = 0;
                    a_Wx(isnan(a_Wx)) = 0;
                A(n,[(n-max_y) (n-1) n (n+max_y) (n+(max_y*max_x))]) = [a_Wx a_Wy a_P a_Ex a_Ez];
                    
                elseif (i==max_y && j==max_x)
%                     B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) =zeros;
                    a_Wx = -sigma(i,j,k)*sigma(i,j-1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(sigma(i,j,k)*0.5*(x(j)-x(j-1))+sigma(i,j-1,k)*0.5*(x(j+1)-x(j)));
                    a_Wy = -sigma(i,j,k)*sigma(i-1,j,k)*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(i,j,k)*0.5*(y(i)-y(i-1))+sigma(i-1,j,k)*0.5*(y(i+1)-y(i)));
                    a_Ez = -sigma(i,j,k)*sigma(i,j,k+1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(sigma(i,j,k+1)*0.5*(z(k+1)-z(k))+sigma(i,j,k)*0.5*(z(k+2)-z(k+1)));
                    a_Wz = -sigma(i,j,k)*sigma(i,j,k-1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(sigma(i,j,k)*0.5*(z(k)-z(k-1))+sigma(i,j,k-1)*0.5*(z(k+1)-z(k)));
                    a_P  =  -(a_Wx + a_Wy + a_Ez + a_Wz);
                    a_Ez(isnan(a_Ez)) = 0;
                    a_P(isnan(a_P)) = 0;
                    a_Wz(isnan(a_Wz)) = 0;
                    a_Wy(isnan(a_Wy)) = 0;
                    a_Wx(isnan(a_Wx)) = 0;
                A(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+(max_y*max_x))]) = [a_Wz a_Wx a_Wy a_P a_Ez];
                    
                elseif (i==max_y && j==1)
%                     B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) =zeros;
                    a_Ex = -sigma(i,j,k)*sigma(i,j+1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(sigma(i,j+1,k)*0.5*(x(j+1)-x(j))+sigma(i,j,k)*0.5*(x(j+2)-x(j+1)));
                    a_Wy = -sigma(i,j,k)*sigma(i-1,j,k)*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(i,j,k)*0.5*(y(i)-y(i-1))+sigma(i-1,j,k)*0.5*(y(i+1)-y(i)));
                    a_Ez = -sigma(i,j,k)*sigma(i,j,k+1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(sigma(i,j,k+1)*0.5*(z(k+1)-z(k))+sigma(i,j,k)*0.5*(z(k+2)-z(k+1)));
                    a_Wz = -sigma(i,j,k)*sigma(i,j,k-1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(sigma(i,j,k)*0.5*(z(k)-z(k-1))+sigma(i,j,k-1)*0.5*(z(k+1)-z(k)));
                    a_P  =  -(a_Ex + a_Wy + a_Ez + a_Wz);
                    a_Ex(isnan(a_Ex)) = 0;
                    a_Ez(isnan(a_Ez)) = 0;
                    a_P(isnan(a_P)) = 0;
                    a_Wz(isnan(a_Wz)) = 0;
                    a_Wy(isnan(a_Wy)) = 0;
                A(n,[(n-(max_y*max_x)) (n-1) n (n+max_y) (n+(max_y*max_x))]) = [a_Wz a_Wy a_P a_Ex a_Ez];
                    
                elseif (i==max_y && k==max_z)
%                     B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) =zeros;
                    a_Ex = -sigma(i,j,k)*sigma(i,j+1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(sigma(i,j+1,k)*0.5*(x(j+1)-x(j))+sigma(i,j,k)*0.5*(x(j+2)-x(j+1)));
                    a_Wx = -sigma(i,j,k)*sigma(i,j-1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(sigma(i,j,k)*0.5*(x(j)-x(j-1))+sigma(i,j-1,k)*0.5*(x(j+1)-x(j)));
                    a_Wy = -sigma(i,j,k)*sigma(i-1,j,k)*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(i,j,k)*0.5*(y(i)-y(i-1))+sigma(i-1,j,k)*0.5*(y(i+1)-y(i)));
                    a_Wz = -sigma(i,j,k)*sigma(i,j,k-1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(sigma(i,j,k)*0.5*(z(k)-z(k-1))+sigma(i,j,k-1)*0.5*(z(k+1)-z(k)));
                    a_P  =  -(a_Ex + a_Wx + a_Wy + a_Wz);
                    a_Ex(isnan(a_Ex)) = 0;
                    a_P(isnan(a_P)) = 0;
                    a_Wz(isnan(a_Wz)) = 0;
                    a_Wy(isnan(a_Wy)) = 0;
                    a_Wx(isnan(a_Wx)) = 0;
                A(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+max_y)]) = [a_Wz a_Wx a_Wy a_P a_Ex];
                    
                    
                elseif (j==max_x && k==1)
%                     B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) =zeros;
                    a_Wx = -sigma(i,j,k)*sigma(i,j-1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(sigma(i,j,k)*0.5*(x(j)-x(j-1))+sigma(i,j-1,k)*0.5*(x(j+1)-x(j)));
                    a_Ey1 = (-sigma(i,j,k)*sigma(a(1),b(1),c(1))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(1),b(1),c(1))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    if c(2) == 0
                        a_Ey2 = 0;
                    else
                        a_Ey2 = (-sigma(i,j,k)*sigma(a(2),b(2),c(2))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(2),b(2),c(2))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    if c(3) == 0 
                        a_Ey3 = 0;
                    else
                    a_Ey3 = (-sigma(i,j,k)*sigma(a(3),b(3),c(3))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(3),b(3),c(3))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    a_Wy = -sigma(i,j,k)*sigma(i+1,j,k)*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(i,j,k)*0.5*(y(i)-y(i-1))+sigma(i-1,j,k)*0.5*(y(i+1)-y(i)));
                    a_Ez = -sigma(i,j,k)*sigma(i,j,k+1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(sigma(i,j,k+1)*0.5*(z(k+1)-z(k))+sigma(i,j,k)*0.5*(z(k+2)-z(k+1)));
                    a_P  =  -(a_Wx + a_Ey1 +a_Ey2 +a_Ey3 + a_Wy + a_Ez);
                    a_Ey1(isnan(a_Ey1)) = 0;
                    a_Ey2(isnan(a_Ey2)) = 0;
                    a_Ey3(isnan(a_Ey3)) = 0;
                    a_Ez(isnan(a_Ez)) = 0;
                    a_P(isnan(a_P)) = 0;
                    a_Wy(isnan(a_Wy)) = 0;
                    a_Wx(isnan(a_Wx)) = 0;
                A(n,[(n-max_y) (n-1) n (nEy1) (nEy2) (nEy2) (n+(max_y*max_x))]) = [a_Wx a_Wy a_P a_Ey1 a_Ey2 a_Ey3 a_Ez];
                    
                elseif (j==max_x && k==max_z)
%                     B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) =zeros;
                    a_Wx = -sigma(i,j,k)*sigma(i,j-1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(sigma(i,j,k)*0.5*(x(j)-x(j-1))+sigma(i,j-1,k)*0.5*(x(j+1)-x(j)));
                    a_Ey1 = (-sigma(i,j,k)*sigma(a(1),b(1),c(1))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(1),b(1),c(1))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    if c(2) == 0
                        a_Ey2 = 0;
                    else
                        a_Ey2 = (-sigma(i,j,k)*sigma(a(2),b(2),c(2))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(2),b(2),c(2))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    if c(3) == 0 
                        a_Ey3 = 0;
                    else
                    a_Ey3 = (-sigma(i,j,k)*sigma(a(3),b(3),c(3))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(3),b(3),c(3))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    a_Wy = -sigma(i,j,k)*sigma(i+1,j,k)*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(i,j,k)*0.5*(y(i)-y(i-1))+sigma(i-1,j,k)*0.5*(y(i+1)-y(i)));
                    a_Wz = -sigma(i,j,k)*sigma(i,j,k-1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(sigma(i,j,k)*0.5*(z(k)-z(k-1))+sigma(i,j,k-1)*0.5*(z(k+1)-z(k)));
                    a_P  =  -(a_Wx + a_Ey1 +a_Ey2 +a_Ey3 + a_Wy + a_Wz);
                    a_Ey1(isnan(a_Ey1)) = 0;
                    a_Ey2(isnan(a_Ey2)) = 0;
                    a_Ey3(isnan(a_Ey3)) = 0;
                    a_P(isnan(a_P)) = 0;
                    a_Wz(isnan(a_Wz)) = 0;
                    a_Wy(isnan(a_Wy)) = 0;
                    a_Wx(isnan(a_Wx)) = 0;
                A(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (nEy1) (nEy2) (nEy2)]) = [a_Wz a_Wx a_Wy a_P a_Ey1 a_Ey2 a_Ey3];
                    
                elseif (j==1 && k==max_z)
%                     B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) =zeros;
                    a_Ex = -sigma(i,j,k)*sigma(i,j+1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(sigma(i,j+1,k)*0.5*(x(j+1)-x(j))+sigma(i,j,k)*0.5*(x(j+2)-x(j+1)));
                    a_Ey1 = (-sigma(i,j,k)*sigma(a(1),b(1),c(1))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(1),b(1),c(1))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    if c(2) == 0
                        a_Ey2 = 0;
                    else
                        a_Ey2 = (-sigma(i,j,k)*sigma(a(2),b(2),c(2))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(2),b(2),c(2))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    if c(3) == 0 
                        a_Ey3 = 0;
                    else
                    a_Ey3 = (-sigma(i,j,k)*sigma(a(3),b(3),c(3))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(3),b(3),c(3))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    a_Wy = -sigma(i,j,k)*sigma(i+1,j,k)*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(i,j,k)*0.5*(y(i)-y(i-1))+sigma(i-1,j,k)*0.5*(y(i+1)-y(i)));
                    a_Wz = -sigma(i,j,k)*sigma(i,j,k-1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(sigma(i,j,k)*0.5*(z(k)-z(k-1))+sigma(i,j,k-1)*0.5*(z(k+1)-z(k)));
                    a_P  =  -(a_Ex + a_Ey1 +a_Ey2 +a_Ey3 + a_Wy + a_Wz);
                    a_Ex(isnan(a_Ex)) = 0;
                    a_Ey1(isnan(a_Ey1)) = 0;
                    a_Ey2(isnan(a_Ey2)) = 0;
                    a_Ey3(isnan(a_Ey3)) = 0;
                    a_P(isnan(a_P)) = 0;
                    a_Wz(isnan(a_Wz)) = 0;
                    a_Wy(isnan(a_Wy)) = 0;
                A(n,[(n-(max_y*max_x)) (n-1) n (nEy1) (nEy2) (nEy2) (n+max_y)]) = [a_Wz a_Wy a_P a_Ey1 a_Ey2 a_Ey3 a_Ex];
                    
                elseif (j==1 && k==1)
%                     B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) =zeros;
                    a_Ex = -sigma(i,j,k)*sigma(i,j+1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(sigma(i,j+1,k)*0.5*(x(j+1)-x(j))+sigma(i,j,k)*0.5*(x(j+2)-x(j+1)));
                    a_Ey1 = (-sigma(i,j,k)*sigma(a(1),b(1),c(1))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(1),b(1),c(1))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    if c(2) == 0
                        a_Ey2 = 0;
                    else
                        a_Ey2 = (-sigma(i,j,k)*sigma(a(2),b(2),c(2))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(2),b(2),c(2))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    if c(3) == 0 
                        a_Ey3 = 0;
                    else
                    a_Ey3 = (-sigma(i,j,k)*sigma(a(3),b(3),c(3))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(3),b(3),c(3))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    a_Wy = -sigma(i,j,k)*sigma(i+1,j,k)*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(i,j,k)*0.5*(y(i)-y(i-1))+sigma(i-1,j,k)*0.5*(y(i+1)-y(i)));
                    a_Ez = -sigma(i,j,k)*sigma(i,j,k+1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(sigma(i,j,k+1)*0.5*(z(k+1)-z(k))+sigma(i,j,k)*0.5*(z(k+2)-z(k+1)));
                    a_P  =  -(a_Ex + a_Ey1 +a_Ey2 +a_Ey3 + a_Wy + a_Ez);
                    a_Ex(isnan(a_Ex)) = 0;
                    a_Ey1(isnan(a_Ey1)) = 0;
                    a_Ey2(isnan(a_Ey2)) = 0;
                    a_Ey3(isnan(a_Ey3)) = 0;
                    a_Ez(isnan(a_Ez)) = 0;
                    a_P(isnan(a_P)) = 0;
                    a_Wy(isnan(a_Wy)) = 0;
                A(n,[(n-1) n (nEy1) (nEy2) (nEy2) (n+max_y) (n+(max_y*max_x))]) = [a_Wy a_P a_Ey1 a_Ey2 a_Ey3 a_Ex a_Ez];
                    
                    %----------------this part defines the surfaces
                    %defining the surface boundaries---------------
                elseif i==1
%                     B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) =zeros;
                    a_Ex = -sigma(i,j,k)*sigma(i,j+1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(sigma(i,j+1,k)*0.5*(x(j+1)-x(j))+sigma(i,j,k)*0.5*(x(j+2)-x(j+1)));
                    a_Wx = -sigma(i,j,k)*sigma(i,j-1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(sigma(i,j,k)*0.5*(x(j)-x(j-1))+sigma(i,j-1,k)*0.5*(x(j+1)-x(j)));
                    a_Ey1 = (-sigma(i,j,k)*sigma(a(1),b(1),c(1))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(1),b(1),c(1))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    if c(2) == 0
                        a_Ey2 = 0;
                    else
                        a_Ey2 = (-sigma(i,j,k)*sigma(a(2),b(2),c(2))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(2),b(2),c(2))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    if c(3) == 0 
                        a_Ey3 = 0;
                    else
                    a_Ey3 = (-sigma(i,j,k)*sigma(a(3),b(3),c(3))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(3),b(3),c(3))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    a_Ez = -sigma(i,j,k)*sigma(i,j,k+1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(sigma(i,j,k+1)*0.5*(z(k+1)-z(k))+sigma(i,j,k)*0.5*(z(k+2)-z(k+1)));
                    a_Wz = -sigma(i,j,k)*sigma(i,j,k-1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(sigma(i,j,k)*0.5*(z(k)-z(k-1))+sigma(i,j,k-1)*0.5*(z(k+1)-z(k)));
                    a_P  =  -(a_Ex + a_Wx + a_Ey1 +a_Ey2 +a_Ey3 + a_Ez + a_Wz);
                    a_Ex(isnan(a_Ex)) = 0;
                    a_Ey1(isnan(a_Ey1)) = 0;
                    a_Ey2(isnan(a_Ey2)) = 0;
                    a_Ey3(isnan(a_Ey3)) = 0;
                    a_Ez(isnan(a_Ez)) = 0;
                    a_P(isnan(a_P)) = 0;
                    a_Wz(isnan(a_Wz)) = 0;
                    a_Wx(isnan(a_Wx)) = 0;
                A(n,[(n-(max_y*max_x)) (n-max_y) n (nEy1) (nEy2) (nEy2) (n+max_y) (n+(max_y*max_x))]) = [a_Wz a_Wx a_P a_Ey1 a_Ey2 a_Ey3 a_Ex a_Ez];
                    
                elseif i==max_y
%                     B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) =zeros;
                    a_Ex = -sigma(i,j,k)*sigma(i,j+1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(sigma(i,j+1,k)*0.5*(x(j+1)-x(j))+sigma(i,j,k)*0.5*(x(j+2)-x(j+1)));
                    a_Wx = -sigma(i,j,k)*sigma(i,j-1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(sigma(i,j,k)*0.5*(x(j)-x(j-1))+sigma(i,j-1,k)*0.5*(x(j+1)-x(j)));
                    a_Wy = -sigma(i,j,k)*sigma(i-1,j,k)*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(i,j,k)*0.5*(y(i)-y(i-1))+sigma(i-1,j,k)*0.5*(y(i+1)-y(i)));
                    a_Ez = -sigma(i,j,k)*sigma(i,j,k+1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(sigma(i,j,k+1)*0.5*(z(k+1)-z(k))+sigma(i,j,k)*0.5*(z(k+2)-z(k+1)));
                    a_Wz = -sigma(i,j,k)*sigma(i,j,k-1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(sigma(i,j,k)*0.5*(z(k)-z(k-1))+sigma(i,j,k-1)*0.5*(z(k+1)-z(k)));
                    a_P  =  -(a_Ex + a_Wx + a_Wy + a_Ez + a_Wz);
                    a_Ex(isnan(a_Ex)) = 0;
                    a_Ez(isnan(a_Ez)) = 0;
                    a_P(isnan(a_P)) = 0;
                    a_Wz(isnan(a_Wz)) = 0;
                    a_Wy(isnan(a_Wy)) = 0;
                    a_Wx(isnan(a_Wx)) = 0;
                A(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+max_y) (n+(max_y*max_x))]) = [a_Wz a_Wx a_Wy a_P a_Ex a_Ez];
                    
                elseif j==1
%                     B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) =zeros;
                    a_Ex = -sigma(i,j,k)*sigma(i,j+1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(sigma(i,j+1,k)*0.5*(x(j+1)-x(j))+sigma(i,j,k)*0.5*(x(j+2)-x(j+1)));
                    a_Ey1 = (-sigma(i,j,k)*sigma(a(1),b(1),c(1))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(1),b(1),c(1))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    if c(2) == 0
                        a_Ey2 = 0;
                    else
                        a_Ey2 = (-sigma(i,j,k)*sigma(a(2),b(2),c(2))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(2),b(2),c(2))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    if c(3) == 0 
                        a_Ey3 = 0;
                    else
                    a_Ey3 = (-sigma(i,j,k)*sigma(a(3),b(3),c(3))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(3),b(3),c(3))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    a_Wy = -sigma(i,j,k)*sigma(i-1,j,k)*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(i,j,k)*0.5*(y(i)-y(i-1))+sigma(i-1,j,k)*0.5*(y(i+1)-y(i)));
                    a_Ez = -sigma(i,j,k)*sigma(i,j,k+1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(sigma(i,j,k+1)*0.5*(z(k+1)-z(k))+sigma(i,j,k)*0.5*(z(k+2)-z(k+1)));
                    a_Wz = -sigma(i,j,k)*sigma(i,j,k-1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(sigma(i,j,k)*0.5*(z(k)-z(k-1))+sigma(i,j,k-1)*0.5*(z(k+1)-z(k)));
                    a_P  =  -(a_Ex + a_Ey1 +a_Ey2 +a_Ey3 + a_Wy + a_Ez + a_Wz);
                    a_Ex(isnan(a_Ex)) = 0;
                    a_Ey1(isnan(a_Ey1)) = 0;
                    a_Ey2(isnan(a_Ey2)) = 0;
                    a_Ey3(isnan(a_Ey3)) = 0;
                    a_Ez(isnan(a_Ez)) = 0;
                    a_P(isnan(a_P)) = 0;
                    a_Wz(isnan(a_Wz)) = 0;
                    a_Wy(isnan(a_Wy)) = 0;
                A(n,[(n-(max_y*max_x)) (n-1) n (nEy1) (nEy2) (nEy2) (n+max_y) (n+(max_y*max_x))]) = [a_Wz a_Wy a_P a_Ey1 a_Ey2 a_Ey3 a_Ex a_Ez];
                    
                elseif j==max_x
%                     B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) =zeros;
                    a_Wx = -sigma(i,j,k)*sigma(i,j-1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(sigma(i,j,k)*0.5*(x(j)-x(j-1))+sigma(i,j-1,k)*0.5*(x(j+1)-x(j)));
                    a_Ey1 = (-sigma(i,j,k)*sigma(a(1),b(1),c(1))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(1),b(1),c(1))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    if c(2) == 0
                        a_Ey2 = 0;
                    else
                        a_Ey2 = (-sigma(i,j,k)*sigma(a(2),b(2),c(2))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(2),b(2),c(2))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    if c(3) == 0 
                        a_Ey3 = 0;
                    else
                    a_Ey3 = (-sigma(i,j,k)*sigma(a(3),b(3),c(3))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(3),b(3),c(3))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    a_Wy = -sigma(i,j,k)*sigma(i-1,j,k)*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(i,j,k)*0.5*(y(i)-y(i-1))+sigma(i-1,j,k)*0.5*(y(i+1)-y(i)));
                    a_Ez = -sigma(i,j,k)*sigma(i,j,k+1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(sigma(i,j,k+1)*0.5*(z(k+1)-z(k))+sigma(i,j,k)*0.5*(z(k+2)-z(k+1)));
                    a_Wz = -sigma(i,j,k)*sigma(i,j,k-1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(sigma(i,j,k)*0.5*(z(k)-z(k-1))+sigma(i,j,k-1)*0.5*(z(k+1)-z(k)));
                    a_P  =  -(a_Wx + a_Ey1 +a_Ey2 +a_Ey3 + a_Wy + a_Ez + a_Wz);
                    a_Ey1(isnan(a_Ey1)) = 0;
                    a_Ey2(isnan(a_Ey2)) = 0;
                    a_Ey3(isnan(a_Ey3)) = 0;
                    a_Ez(isnan(a_Ez)) = 0;
                    a_P(isnan(a_P)) = 0;
                    a_Wz(isnan(a_Wz)) = 0;
                    a_Wy(isnan(a_Wy)) = 0;
                    a_Wx(isnan(a_Wx)) = 0;
                A(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (nEy1) (nEy2) (nEy2) (n+(max_y*max_x))]) = [a_Wz a_Wx a_Wy a_P a_Ey1 a_Ey2 a_Ey3 a_Ez];
                    
                elseif k==1
%                     B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) =zeros;
                    a_Ex = -sigma(i,j,k)*sigma(i,j+1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(sigma(i,j+1,k)*0.5*(x(j+1)-x(j))+sigma(i,j,k)*0.5*(x(j+2)-x(j+1)));
                    a_Wx = -sigma(i,j,k)*sigma(i,j-1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(sigma(i,j,k)*0.5*(x(j)-x(j-1))+sigma(i,j-1,k)*0.5*(x(j+1)-x(j)));
                    a_Ey1 = (-sigma(i,j,k)*sigma(a(1),b(1),c(1))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(1),b(1),c(1))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    if c(2) == 0
                        a_Ey2 = 0;
                    else
                        a_Ey2 = (-sigma(i,j,k)*sigma(a(2),b(2),c(2))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(2),b(2),c(2))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    if c(3) == 0 
                        a_Ey3 = 0;
                    else
                    a_Ey3 = (-sigma(i,j,k)*sigma(a(3),b(3),c(3))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(3),b(3),c(3))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    a_Wy = -sigma(i,j,k)*sigma(i-1,j,k)*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(i,j,k)*0.5*(y(i)-y(i-1))+sigma(i-1,j,k)*0.5*(y(i+1)-y(i)));
                    a_Ez = -sigma(i,j,k)*sigma(i,j,k+1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(sigma(i,j,k+1)*0.5*(z(k+1)-z(k))+sigma(i,j,k)*0.5*(z(k+2)-z(k+1)));
                    a_P  =  -(a_Ex + a_Wx + a_Ey1 +a_Ey2 +a_Ey3 + a_Wy + a_Ez);
                    a_Ex(isnan(a_Ex)) = 0;
                    a_Ey1(isnan(a_Ey1)) = 0;
                    a_Ey2(isnan(a_Ey2)) = 0;
                    a_Ey3(isnan(a_Ey3)) = 0;
                    a_Ez(isnan(a_Ez)) = 0;
                    a_P(isnan(a_P)) = 0;
                    a_Wy(isnan(a_Wy)) = 0;
                    a_Wx(isnan(a_Wx)) = 0;
                A(n,[(n-max_y) (n-1) n (nEy1) (nEy2) (nEy2) (n+max_y) (n+(max_y*max_x))]) = [a_Wx a_Wy a_P a_Ey1 a_Ey2 a_Ey3 a_Ex a_Ez];
                
                elseif k==max_z
%                     B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) =zeros;
                    a_Ex = -sigma(i,j,k)*sigma(i,j+1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(sigma(i,j+1,k)*0.5*(x(j+1)-x(j))+sigma(i,j,k)*0.5*(x(j+2)-x(j+1)));
                    a_Wx = -sigma(i,j,k)*sigma(i,j-1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(sigma(i,j,k)*0.5*(x(j)-x(j-1))+sigma(i,j-1,k)*0.5*(x(j+1)-x(j)));
                    a_Ey1 = (-sigma(i,j,k)*sigma(a(1),b(1),c(1))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(1),b(1),c(1))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    if c(2) == 0
                        a_Ey2 = 0;
                    else
                        a_Ey2 = (-sigma(i,j,k)*sigma(a(2),b(2),c(2))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(2),b(2),c(2))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    if c(3) == 0 
                        a_Ey3 = 0;
                    else
                    a_Ey3 = (-sigma(i,j,k)*sigma(a(3),b(3),c(3))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(3),b(3),c(3))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    a_Wy = -sigma(i,j,k)*sigma(i-1,j,k)*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(i,j,k)*0.5*(y(i)-y(i-1))+sigma(i-1,j,k)*0.5*(y(i+1)-y(i)));
                    a_Wz = -sigma(i,j,k)*sigma(i,j,k-1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(sigma(i,j,k)*0.5*(z(k)-z(k-1))+sigma(i,j,k-1)*0.5*(z(k+1)-z(k)));
                    a_P  =  -(a_Ex + a_Wx + a_Ey1 +a_Ey2 +a_Ey3 + a_Wy + a_Wz);
                    a_Ex(isnan(a_Ex)) = 0;
                    a_Ey1(isnan(a_Ey1)) = 0;
                    a_Ey2(isnan(a_Ey2)) = 0;
                    a_Ey3(isnan(a_Ey3)) = 0;
                    a_P(isnan(a_P)) = 0;
                    a_Wz(isnan(a_Wz)) = 0;
                    a_Wy(isnan(a_Wy)) = 0;
                    a_Wx(isnan(a_Wx)) = 0;
                A(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (nEy1) (nEy2) (nEy2) (n+max_y)]) = [a_Wz a_Wx a_Wy a_P a_Ey1 a_Ey2 a_Ey3 a_Ex];
                
%                 elseif (i>1 && j>1 && k>1) && (i<max_y && j<max_x && k<max_z)
% %                    B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) =zeros;
%                     a_Ex = -sigma(i,j,k)*sigma(i,j+1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(sigma(i,j+1,k)*0.5*(x(j+1)-x(j))+sigma(i,j,k)*0.5*(x(j+2)-x(j+1)));
%                     a_Wx = -sigma(i,j,k)*sigma(i,j-1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(sigma(i,j,k)*0.5*(x(j)-x(j-1))+sigma(i,j-1,k)*0.5*(x(j+1)-x(j)));
%                     a_Ey = -sigma(i,j,k)*sigma(a(1),b(1),c(1))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(1),b(1),c(1))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1)));
%                     a_Wy = -sigma(i,j,k)*sigma(i-1,j,k)*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(i,j,k)*0.5*(y(i)-y(i-1))+sigma(i-1,j,k)*0.5*(y(i+1)-y(i)));
%                     a_Ez = -sigma(i,j,k)*sigma(i,j,k+1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(sigma(i,j,k+1)*0.5*(z(k+1)-z(k))+sigma(i,j,k)*0.5*(z(k+2)-z(k+1)));
%                     a_Wz = -sigma(i,j,k)*sigma(i,j,k-1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(sigma(i,j,k)*0.5*(z(k)-z(k-1))+sigma(i,j,k-1)*0.5*(z(k+1)-z(k)));
%                     a_P  =  -(a_Ex + a_Wx + a_Ey + a_Wy + a_Ez + a_Wz);
%                     a_Ex(isnan(a_Ex)) = 0;
%                     a_Ey(isnan(a_Ey)) = 0;
%                     a_Ez(isnan(a_Ez)) = 0;
%                     a_P(isnan(a_P)) = 0;
%                     a_Wz(isnan(a_Wz)) = 0;
%                     a_Wy(isnan(a_Wy)) = 0;
%                     a_Wx(isnan(a_Wx)) = 0;
%                     
% %                 B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) = [a_Wz a_Wx a_Wy a_P a_Ey a_Ex a_Ez];
%                 A(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (nEy1) (n+max_y) (n+(max_y*max_x))]) = [a_Wz a_Wx a_Wy a_P a_Ey a_Ex a_Ez];
%                 
               elseif (i>1 && j>1 && k>1) && (i<max_y && j<max_x && k<max_z)
%                    B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) =zeros;
                    a_Ex = -wsize(1)*sigma(i,j,k)*sigma(i,j+1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(sigma(i,j+1,k)*0.5*(x(j+1)-x(j))+sigma(i,j,k)*0.5*(x(j+2)-x(j+1)));
                    a_Wx = -wsize(1)*sigma(i,j,k)*sigma(i,j-1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(sigma(i,j,k)*0.5*(x(j)-x(j-1))+sigma(i,j-1,k)*0.5*(x(j+1)-x(j)));
                    a_Ey1 = (-sigma(i,j,k)*sigma(a(1),b(1),c(1))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(1),b(1),c(1))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    if c(2) == 0
                        a_Ey2 = 0;
                    else
                        a_Ey2 = (-sigma(i,j,k)*sigma(a(2),b(2),c(2))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(2),b(2),c(2))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    if c(3) == 0 
                        a_Ey3 = 0;
                    else
                    a_Ey3 = (-sigma(i,j,k)*sigma(a(3),b(3),c(3))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(3),b(3),c(3))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    a_Wy = -wsize(1)*sigma(i,j,k)*sigma(i-1,j,k)*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(i,j,k)*0.5*(y(i)-y(i-1))+sigma(i-1,j,k)*0.5*(y(i+1)-y(i)));
                    a_Ez = -wsize(1)*sigma(i,j,k)*sigma(i,j,k+1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(sigma(i,j,k+1)*0.5*(z(k+1)-z(k))+sigma(i,j,k)*0.5*(z(k+2)-z(k+1)));
                    a_Wz = -wsize(1)*sigma(i,j,k)*sigma(i,j,k-1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(sigma(i,j,k)*0.5*(z(k)-z(k-1))+sigma(i,j,k-1)*0.5*(z(k+1)-z(k)));
                    a_P  =  -wsize(1)*(a_Ex + a_Wx + a_Ey1 +a_Ey2 +a_Ey3 + a_Wy + a_Ez + a_Wz);
                    a_Ex(isnan(a_Ex)) = 0;
                    a_Ey(isnan(a_Ey)) = 0;
                    a_Ez(isnan(a_Ez)) = 0;
                    a_P(isnan(a_P)) = 0;
                    a_Wz(isnan(a_Wz)) = 0;
                    a_Wy(isnan(a_Wy)) = 0;
                    a_Wx(isnan(a_Wx)) = 0;
                    
%                 B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) = [a_Wz a_Wx a_Wy a_P a_Ey a_Ex a_Ez];
                A(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (nEy1) (nEy2) (nEy2) (n+max_y) (n+(max_y*max_x))]) = [a_Wz a_Wx a_Wy a_P a_Ey1 a_Ey2 a_Ey3 a_Ex a_Ez];
                    
                end           
    
    
end
