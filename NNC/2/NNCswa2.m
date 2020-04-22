[max_y, max_x, max_z]=size(sigma); 
% ne = max_y*max_x*max_z;
% B = sparse(ne,ne);
w = 1;
len = length(NNCI1a)+1;


while w < len
    
    
%     check = find(NNC.NNC2 == NNC.NNC2(w));
    check = find(NNCI2a == NNCI2a(w));
       
    w = check(end);
    wsize = size(check);
    
    NNCOV;
    
    i = NNCI2a_dex(w,1);
    j = NNCI2a_dex(w,2);
    k = NNCI2a_dex(w,3);  
    
    a = zeros(5,1);
    b = zeros(5,1);
    c = zeros(5,1);
    
        
    if wsize(1) > 1
        
        for wcheck = 1:wsize(1)
                    
            a(wcheck) = NNCI1a_dex(w,1);
            b(wcheck) = NNCI1a_dex(w,2);
            c(wcheck) = NNCI1a_dex(w,3); 
            w = w+1;            
        end
        
    else
            a(1) = NNCI1a_dex(w,1);
            b(1) = NNCI1a_dex(w,2);
            c(1) = NNCI1a_dex(w,3); 
            
            w = w+1;
    end
    
    n=(k-1)*(max_y*max_x)+(j-1)*max_y+i;
    nWy1 = (c(1)-1)*(max_y*max_x)+(b(1)-1)*max_y+a(1); 
    nWy2 = (c(2)-1)*(max_y*max_x)+(b(2)-1)*max_y+a(2);
    nWy2(nWy2<0)=n-1;
    nWy3 = (c(3)-1)*(max_y*max_x)+(b(3)-1)*max_y+a(3);
    nWy3(nWy3<0)=n-1;
    nWy4 = (c(4)-1)*(max_y*max_x)+(b(4)-1)*max_y+a(4);
    nWy4(nWy4<0)=n-1;
    nWy5 = (c(5)-1)*(max_y*max_x)+(b(5)-1)*max_y+a(5);
    nWy5(nWy5<0)=n-1;
                
                %----------------defines corners of resistance matrix ---------
                 
                
     if (j==1 && i==1) && k==1
%                     B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) = zeros;
                    a_Ex = -sigma(i,j,k)*sigma(i,j+1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(sigma(i,j+1,k)*0.5*(x(j+1)-x(j))+sigma(i,j,k)*0.5*(x(j+2)-x(j+1)));
                    a_Ey = -sigma(i,j,k)*sigma(i+1,j,k)*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(i+1,j,k)*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1)));
                    a_Ez = -sigma(i,j,k)*sigma(i,j,k+1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(sigma(i,j,k+1)*0.5*(z(k+1)-z(k))+sigma(i,j,k)*0.5*(z(k+2)-z(k+1)));
                    a_P  =  -(a_Ex + (a_Ey+a_Ey2+a_Ey3) + a_Ez);
                    a_Ex(isnan(a_Ex)) = 0;
                    a_Ey(isnan(a_Ey)) = 0;
                    a_Ez(isnan(a_Ez)) = 0;
                    a_P(isnan(a_P)) = 0;
                A(n,[n (n+1) (n+max_y) (n+(max_y*max_x))]) = [a_P a_Ey a_Ex a_Ez];
                    
                elseif (j==max_x && i==max_y) && k==max_z
%                   B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) = zeros;
                    a_Wx = -wsize(1)*sigma(i,j,k)*sigma(i,j-1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(sigma(i,j,k)*0.5*(x(j)-x(j-1))+sigma(i,j-1,k)*0.5*(x(j+1)-x(j)));
                    a_Wy1 = (-sigma(i,j,k)*sigma(a(1),b(1),c(1))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(1),b(1),c(1))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    if c(2) == 0
                        a_Wy2 = 0;
                    else
                        a_Wy2 = (-sigma(i,j,k)*sigma(a(2),b(2),c(2))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(2),b(2),c(2))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    if c(3) == 0 
                        a_Wy3 = 0;
                    else
                    a_Wy3 = (-sigma(i,j,k)*sigma(a(3),b(3),c(3))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(3),b(3),c(3))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    if c(4) == 0 
                        a_Wy4 = 0;
                    else
                    a_Wy4 = (-sigma(i,j,k)*sigma(a(4),b(4),c(4))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(4),b(4),c(4))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    if c(5) == 0 
                        a_Wy5 = 0;
                    else
                    a_Wy5 = (-sigma(i,j,k)*sigma(a(5),b(5),c(5))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(5),b(5),c(5))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    a_Wz = -wsize(1)*sigma(i,j,k)*sigma(i,j,k-1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(sigma(i,j,k)*0.5*(z(k)-z(k-1))+sigma(i,j,k-1)*0.5*(z(k+1)-z(k)));
                    a_P  =  -wsize(1)*(a_Wx +a_Wy1  +a_Wy2 +a_Wy3 +a_Wy4 +a_Wy5 + a_Wz);
                    a_P(isnan(a_P)) = 0;
                    a_Wz(isnan(a_Wz)) = 0;
                    a_Wy1(isnan(a_Wy1)) = 0;
                    a_Wy2(isnan(a_Wy2)) = 0;
                    a_Wy3(isnan(a_Wy3)) = 0;
                    a_Wy4(isnan(a_Wy4)) = 0;
                    a_Wy5(isnan(a_Wy5)) = 0;
                    a_Wx(isnan(a_Wx)) = 0;
                A(n,[n-(max_y*max_x) n-max_y nWy1  nWy2 nWy3 nWy4 nWy5 n]) = [a_Wz a_Wx a_Wy1  a_Wy2 a_Wy3 a_Wy4 a_Wy5 a_P];     
                    
                elseif (j==1 && i==1) && k==max_z
%                   B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) = zeros; 
                    a_Ex = -sigma(i,j,k)*sigma(i,j+1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(sigma(i,j+1,k)*0.5*(x(j+1)-x(j))+sigma(i,j,k)*0.5*(x(j+2)-x(j+1)));
                    a_Ey = -sigma(i,j,k)*sigma(i+1,j,k)*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(i+1,j,k)*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1)));
                    a_Wz = -sigma(i,j,k)*sigma(i,j,k-1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(sigma(i,j,k)*0.5*(z(k)-z(k-1))+sigma(i,j,k-1)*0.5*(z(k+1)-z(k)));
                    a_P  =  -(a_Ex + a_Ey + a_Wz);
                    a_Ex(isnan(a_Ex)) = 0;
                    a_Ey(isnan(a_Ey)) = 0;
                    a_P(isnan(a_P)) = 0;
                    a_Wz(isnan(a_Wz)) = 0;
                A(n,[(n-(max_y*max_x)) n (n+1) (n+max_y)]) = [a_Wz a_P a_Ey a_Ex]; 
                
                elseif (j==max_x && i==1) && k==max_z
%                     B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) = zeros;
                    a_Wx = -sigma(i,j,k)*sigma(i,j-1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(sigma(i,j,k)*0.5*(x(j)-x(j-1))+sigma(i,j-1,k)*0.5*(x(j+1)-x(j)));
                    a_Ey = -sigma(i,j,k)*sigma(i+1,j,k)*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(i+1,j,k)*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1)));
                    a_Wz = -sigma(i,j,k)*sigma(i,j,k-1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(sigma(i,j,k)*0.5*(z(k)-z(k-1))+sigma(i,j,k-1)*0.5*(z(k+1)-z(k)));
                    a_P  =  -(a_Wx + a_Ey + a_Wz);
                    a_Ey(isnan(a_Ey)) = 0;
                    a_P(isnan(a_P)) = 0;
                    a_Wz(isnan(a_Wz)) = 0;
                    a_Wx(isnan(a_Wx)) = 0;
                A(n,[(n-(max_y*max_x)) (n-max_y) n (n+1)]) = [a_Wz a_Wx a_P a_Ey];
                    
                elseif (j==max_x && i==1) && k==1
%                     B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) = zeros;
                    a_Wx = -sigma(i,j,k)*sigma(i,j-1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(sigma(i,j,k)*0.5*(x(j)-x(j-1))+sigma(i,j-1,k)*0.5*(x(j+1)-x(j)));
                    a_Ey = -sigma(i,j,k)*sigma(i+1,j,k)*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(i+1,j,k)*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1)));
                    a_Ez = -sigma(i,j,k)*sigma(i,j,k+1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(sigma(i,j,k+1)*0.5*(z(k+1)-z(k))+sigma(i,j,k)*0.5*(z(k+2)-z(k+1)));
                    a_P  =  -(a_Wx + a_Ey + a_Ez);
                    a_Ey(isnan(a_Ey)) = 0;
                    a_Ez(isnan(a_Ez)) = 0;
                    a_P(isnan(a_P)) = 0;
                    a_Wx(isnan(a_Wx)) = 0;
                A(n,[(n-max_y) n (n+max_y) (n+(max_y*max_x))]) = [a_Wx a_P a_Ey a_Ez];
                    
                elseif (j==max_x && i==max_y) && k==1
%                     B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) = zeros;
                    a_Wx = -wsize(1)*sigma(i,j,k)*sigma(i,j-1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(sigma(i,j,k)*0.5*(x(j)-x(j-1))+sigma(i,j-1,k)*0.5*(x(j+1)-x(j)));
                    a_Wy1 = (-sigma(i,j,k)*sigma(a(1),b(1),c(1))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(1),b(1),c(1))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    if c(2) == 0
                        a_Wy2 = 0;
                    else
                        a_Wy2 = (-sigma(i,j,k)*sigma(a(2),b(2),c(2))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(2),b(2),c(2))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    if c(3) == 0 
                        a_Wy3 = 0;
                    else
                    a_Wy3 = (-sigma(i,j,k)*sigma(a(3),b(3),c(3))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(3),b(3),c(3))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    if c(4) == 0 
                        a_Wy4 = 0;
                    else
                    a_Wy4 = (-sigma(i,j,k)*sigma(a(4),b(4),c(4))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(4),b(4),c(4))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    if c(5) == 0 
                        a_Wy5 = 0;
                    else
                    a_Wy5 = (-sigma(i,j,k)*sigma(a(5),b(5),c(5))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(5),b(5),c(5))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    a_Ez = -wsize(1)*sigma(i,j,k)*sigma(i,j,k+1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(sigma(i,j,k+1)*0.5*(z(k+1)-z(k))+sigma(i,j,k)*0.5*(z(k+2)-z(k+1)));
                    a_P  =  -wsize(1)*(a_Wx +a_Wy1  +a_Wy2 +a_Wy3 +a_Wy4 +a_Wy5 + a_Ez);
                    a_Ez(isnan(a_Ez)) = 0;
                    a_P(isnan(a_P)) = 0;
                    a_Wy1(isnan(a_Wy1)) = 0;
                    a_Wy2(isnan(a_Wy2)) = 0;
                    a_Wy3(isnan(a_Wy3)) = 0;
                    a_Wy4(isnan(a_Wy4)) = 0;
                    a_Wy5(isnan(a_Wy5)) = 0;
                    a_Wx(isnan(a_Wx)) = 0;
                A(n,[(n-max_y) nWy1  nWy2 nWy3 nWy4 nWy5 n (n+(max_y*max_x))]) = [a_Wx a_Wy1  a_Wy2 a_Wy3 a_Wy4 a_Wy5 a_P a_Ez];
                    
                elseif (j==1 && i==max_y) && k==max_z
%                     B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) = zeros;
                    a_Ex = -wsize(1)*sigma(i,j,k)*sigma(i,j+1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(sigma(i,j+1,k)*0.5*(x(j+1)-x(j))+sigma(i,j,k)*0.5*(x(j+2)-x(j+1)));
                    a_Wy1 = (-sigma(i,j,k)*sigma(a(1),b(1),c(1))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(1),b(1),c(1))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    if c(2) == 0
                        a_Wy2 = 0;
                    else
                        a_Wy2 = (-sigma(i,j,k)*sigma(a(2),b(2),c(2))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(2),b(2),c(2))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    if c(3) == 0 
                        a_Wy3 = 0;
                    else
                    a_Wy3 = (-sigma(i,j,k)*sigma(a(3),b(3),c(3))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(3),b(3),c(3))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    if c(4) == 0 
                        a_Wy4 = 0;
                    else
                    a_Wy4 = (-sigma(i,j,k)*sigma(a(4),b(4),c(4))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(4),b(4),c(4))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    if c(5) == 0 
                        a_Wy5 = 0;
                    else
                    a_Wy5 = (-sigma(i,j,k)*sigma(a(5),b(5),c(5))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(5),b(5),c(5))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    a_Wz = -wsize(1)*sigma(i,j,k)*sigma(i,j,k-1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(sigma(i,j,k)*0.5*(z(k)-z(k-1))+sigma(i,j,k-1)*0.5*(z(k+1)-z(k)));
                    a_P  =  -wsize(1)*(a_Ex +a_Wy1  +a_Wy2 +a_Wy3 +a_Wy4 +a_Wy5 + a_Wz);
                    a_Ex(isnan(a_Ex)) = 0;
                    a_P(isnan(a_P)) = 0;
                    a_Wz(isnan(a_Wz)) = 0;
                    a_Wy1(isnan(a_Wy1)) = 0;
                    a_Wy2(isnan(a_Wy2)) = 0;
                    a_Wy3(isnan(a_Wy3)) = 0;
                    a_Wy4(isnan(a_Wy4)) = 0;
                    a_Wy5(isnan(a_Wy5)) = 0;
                    A(n,[(n-(max_y*max_x)) nWy1  nWy2 nWy3 nWy4 nWy5 n (n+max_y)]) = [a_Wz a_Wy1  a_Wy2 a_Wy3 a_Wy4 a_Wy5 a_P a_Ex];
                    
                elseif (j==1 && i==max_y) && k==1
%                     B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) = zeros;
                    a_Ex = -wsize(1)*sigma(i,j,k)*sigma(i,j+1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(sigma(i,j+1,k)*0.5*(x(j+1)-x(j))+sigma(i,j,k)*0.5*(x(j+2)-x(j+1)));
                    a_Wy1 = (-sigma(i,j,k)*sigma(a(1),b(1),c(1))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(1),b(1),c(1))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    if c(2) == 0
                        a_Wy2 = 0;
                    else
                        a_Wy2 = (-sigma(i,j,k)*sigma(a(2),b(2),c(2))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(2),b(2),c(2))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    if c(3) == 0 
                        a_Wy3 = 0;
                    else
                    a_Wy3 = (-sigma(i,j,k)*sigma(a(3),b(3),c(3))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(3),b(3),c(3))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    if c(4) == 0 
                        a_Wy4 = 0;
                    else
                    a_Wy4 = (-sigma(i,j,k)*sigma(a(4),b(4),c(4))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(4),b(4),c(4))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    if c(5) == 0 
                        a_Wy5 = 0;
                    else
                    a_Wy5 = (-sigma(i,j,k)*sigma(a(5),b(5),c(5))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(5),b(5),c(5))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    a_Ez = -wsize(1)*sigma(i,j,k)*sigma(i,j,k+1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(sigma(i,j,k+1)*0.5*(z(k+1)-z(k))+sigma(i,j,k)*0.5*(z(k+2)-z(k+1)));
                    a_P  =  -wsize(1)*(a_Ex +a_Wy1  +a_Wy2 +a_Wy3 +a_Wy4 +a_Wy5 + a_Ez);
                    a_Ex(isnan(a_Ex)) = 0;
                    a_Ez(isnan(a_Ez)) = 0;
                    a_P(isnan(a_P)) = 0;
                    a_Wy1(isnan(a_Wy1)) = 0;
                    a_Wy2(isnan(a_Wy2)) = 0;
                    a_Wy3(isnan(a_Wy3)) = 0;                    
                    a_Wy4(isnan(a_Wy4)) = 0;
                    a_Wy5(isnan(a_Wy5)) = 0;
                A(n,[nWy1  nWy2 nWy3 nWy4 nWy5 n (n+max_y) (n+(max_y*max_x))]) = [a_Wy1  a_Wy2 a_Wy3 a_Wy4 a_Wy5 a_P a_Ex a_Ez];
                    
                    %----------------defines edges of resistance matrix -----------
                elseif (i==1 && k==1)
%                     B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) = zeros;
                    a_Ex = -sigma(i,j,k)*sigma(i,j+1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(sigma(i,j+1,k)*0.5*(x(j+1)-x(j))+sigma(i,j,k)*0.5*(x(j+2)-x(j+1)));
                    a_Wx = -sigma(i,j,k)*sigma(i,j-1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(sigma(i,j,k)*0.5*(x(j)-x(j-1))+sigma(i,j-1,k)*0.5*(x(j+1)-x(j)));
                    a_Ey = -sigma(i,j,k)*sigma(i+1,j,k)*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(i+1,j,k)*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1)));
                    a_Ez = -sigma(i,j,k)*sigma(i,j,k+1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(sigma(i,j,k+1)*0.5*(z(k+1)-z(k))+sigma(i,j,k)*0.5*(z(k+2)-z(k+1)));
                    a_P  =  -(a_Ex + a_Wx + a_Ey + a_Ez);
                    a_Ex(isnan(a_Ex)) = 0;
                    a_Ey(isnan(a_Ey)) = 0;
                    a_Ez(isnan(a_Ez)) = 0;
                    a_P(isnan(a_P)) = 0;
                A(n,[(n-max_y) n (n+1) (n+max_y) (n+(max_y*max_x))]) = [a_Wx a_P a_Ey a_Ex a_Ez];
                
                elseif (i==1 && j==1)
%                     B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) = zeros;
                    a_Ex = -sigma(i,j,k)*sigma(i,j+1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(sigma(i,j+1,k)*0.5*(x(j+1)-x(j))+sigma(i,j,k)*0.5*(x(j+2)-x(j+1)));
                    a_Ey = -sigma(i,j,k)*sigma(i+1,j,k)*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(i+1,j,k)*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1)));
                    a_Ez = -sigma(i,j,k)*sigma(i,j,k+1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(sigma(i,j,k+1)*0.5*(z(k+1)-z(k))+sigma(i,j,k)*0.5*(z(k+2)-z(k+1)));
                    a_Wz = -sigma(i,j,k)*sigma(i,j,k-1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(sigma(i,j,k)*0.5*(z(k)-z(k-1))+sigma(i,j,k-1)*0.5*(z(k+1)-z(k)));
                    a_P  =  -(a_Ex + a_Ey + a_Ez + a_Wz);
                    a_Ex(isnan(a_Ex)) = 0;
                    a_Ey(isnan(a_Ey)) = 0;
                    a_Ez(isnan(a_Ez)) = 0;
                    a_P(isnan(a_P)) = 0;
                    a_Wz(isnan(a_Wz)) = 0;
                A(n,[(n-(max_y*max_x)) n (n+1) (n+max_y) (n+(max_y*max_x))]) = [a_Wz a_P a_Ey a_Ex a_Ez];
                    
                elseif (i==1 && k==max_z)
%                     B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) = zeros;
                    a_Ex = -sigma(i,j,k)*sigma(i,j+1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(sigma(i,j+1,k)*0.5*(x(j+1)-x(j))+sigma(i,j,k)*0.5*(x(j+2)-x(j+1)));
                    a_Wx = -sigma(i,j,k)*sigma(i,j-1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(sigma(i,j,k)*0.5*(x(j)-x(j-1))+sigma(i,j-1,k)*0.5*(x(j+1)-x(j)));
                    a_Ey = -sigma(i,j,k)*sigma(i+1,j,k)*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(i+1,j,k)*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1)));
                    a_Wz = -sigma(i,j,k)*sigma(i,j,k-1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(sigma(i,j,k)*0.5*(z(k)-z(k-1))+sigma(i,j,k-1)*0.5*(z(k+1)-z(k)));
                    a_P  =  -(a_Ex + a_Wx + a_Ey + a_Wz);
                    a_Ex(isnan(a_Ex)) = 0;
                    a_Ey(isnan(a_Ey)) = 0;
                    a_P(isnan(a_P)) = 0;
                    a_Wz(isnan(a_Wz)) = 0;
                    a_Wx(isnan(a_Wx)) = 0;
                A(n,[(n-(max_y*max_x)) (n-max_y) n (n+1) (n+max_y)]) = [a_Wz a_Wx a_P a_Ey a_Ex];
                    
                elseif (i==1 && j==max_x)
%                     B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) = zeros;
                    a_Wx = -sigma(i,j,k)*sigma(i,j-1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(sigma(i,j,k)*0.5*(x(j)-x(j-1))+sigma(i,j-1,k)*0.5*(x(j+1)-x(j)));
                    a_Ey = -sigma(i,j,k)*sigma(i+1,j,k)*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(i+1,j,k)*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1)));
                    a_Ez = -sigma(i,j,k)*sigma(i,j,k+1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(sigma(i,j,k+1)*0.5*(z(k+1)-z(k))+sigma(i,j,k)*0.5*(z(k+2)-z(k+1)));
                    a_Wz = -sigma(i,j,k)*sigma(i,j,k-1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(sigma(i,j,k)*0.5*(z(k)-z(k-1))+sigma(i,j,k-1)*0.5*(z(k+1)-z(k)));
                    a_P  =  -(a_Wx + a_Ey + a_Ez + a_Wz);
                    a_Ey(isnan(a_Ey)) = 0;
                    a_Ez(isnan(a_Ez)) = 0;
                    a_P(isnan(a_P)) = 0;
                    a_Wz(isnan(a_Wz)) = 0;
                    a_Wx(isnan(a_Wx)) = 0;
                A(n,[(n-(max_y*max_x)) (n-max_y) n (n+1) (n+(max_y*max_x))]) = [a_Wz a_Wx a_P a_Ey a_Ez];
                   
                elseif (i==max_y && k==1)
%                     B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) = zeros;
                    a_Ex = -wsize(1)*sigma(i,j,k)*sigma(i,j+1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(sigma(i,j+1,k)*0.5*(x(j+1)-x(j))+sigma(i,j,k)*0.5*(x(j+2)-x(j+1)));
                    a_Wx = -wsize(1)*sigma(i,j,k)*sigma(i,j-1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(sigma(i,j,k)*0.5*(x(j)-x(j-1))+sigma(i,j-1,k)*0.5*(x(j+1)-x(j)));
                    a_Wy1 = (-sigma(i,j,k)*sigma(a(1),b(1),c(1))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(1),b(1),c(1))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    if c(2) == 0
                        a_Wy2 = 0;
                    else
                        a_Wy2 = (-sigma(i,j,k)*sigma(a(2),b(2),c(2))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(2),b(2),c(2))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    if c(3) == 0 
                        a_Wy3 = 0;
                    else
                    a_Wy3 = (-sigma(i,j,k)*sigma(a(3),b(3),c(3))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(3),b(3),c(3))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    if c(4) == 0 
                        a_Wy4 = 0;
                    else
                    a_Wy4 = (-sigma(i,j,k)*sigma(a(4),b(4),c(4))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(4),b(4),c(4))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    if c(5) == 0 
                        a_Wy5 = 0;
                    else
                    a_Wy5 = (-sigma(i,j,k)*sigma(a(5),b(5),c(5))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(5),b(5),c(5))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    a_Ez = -wsize(1)*sigma(i,j,k)*sigma(i,j,k+1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(sigma(i,j,k+1)*0.5*(z(k+1)-z(k))+sigma(i,j,k)*0.5*(z(k+2)-z(k+1)));
                    a_P  =  -wsize(1)*(a_Ex + a_Wx +a_Wy1  +a_Wy2 +a_Wy3 +a_Wy4 +a_Wy5 + a_Ez);
                    a_Ex(isnan(a_Ex)) = 0;
                    a_Ez(isnan(a_Ez)) = 0;
                    a_P(isnan(a_P)) = 0;
                    a_Wy1(isnan(a_Wy1)) = 0;
                    a_Wy2(isnan(a_Wy2)) = 0;
                    a_Wy3(isnan(a_Wy3)) = 0;
                    a_Wy4(isnan(a_Wy4)) = 0;
                    a_Wy5(isnan(a_Wy5)) = 0;
                    a_Wx(isnan(a_Wx)) = 0;
                A(n,[(n-max_y) nWy1  nWy2 nWy3 nWy4 nWy5 n (n+max_y) (n+(max_y*max_x))]) = [a_Wx a_Wy1  a_Wy2 a_Wy3 a_Wy4 a_Wy5 a_P a_Ex a_Ez];
                    
                elseif (i==max_y && j==max_x)
%                     B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) = zeros;
                    a_Wx = -wsize(1)*sigma(i,j,k)*sigma(i,j-1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(sigma(i,j,k)*0.5*(x(j)-x(j-1))+sigma(i,j-1,k)*0.5*(x(j+1)-x(j)));
                    a_Wy1 = (-sigma(i,j,k)*sigma(a(1),b(1),c(1))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(1),b(1),c(1))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    if c(2) == 0
                        a_Wy2 = 0;
                    else
                        a_Wy2 = (-sigma(i,j,k)*sigma(a(2),b(2),c(2))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(2),b(2),c(2))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    if c(3) == 0 
                        a_Wy3 = 0;
                    else
                    a_Wy3 = (-sigma(i,j,k)*sigma(a(3),b(3),c(3))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(3),b(3),c(3))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    if c(4) == 0 
                        a_Wy4 = 0;
                    else
                    a_Wy4 = (-sigma(i,j,k)*sigma(a(4),b(4),c(4))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(4),b(4),c(4))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    if c(5) == 0 
                        a_Wy5 = 0;
                    else
                    a_Wy5 = (-sigma(i,j,k)*sigma(a(5),b(5),c(5))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(5),b(5),c(5))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    a_Ez = -wsize(1)*sigma(i,j,k)*sigma(i,j,k+1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(sigma(i,j,k+1)*0.5*(z(k+1)-z(k))+sigma(i,j,k)*0.5*(z(k+2)-z(k+1)));
                    a_Wz = -wsize(1)*sigma(i,j,k)*sigma(i,j,k-1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(sigma(i,j,k)*0.5*(z(k)-z(k-1))+sigma(i,j,k-1)*0.5*(z(k+1)-z(k)));
                    a_P  =  -wsize(1)*(a_Wx +a_Wy1  +a_Wy2 +a_Wy3 +a_Wy4 +a_Wy5 + a_Ez + a_Wz);
                    a_Ez(isnan(a_Ez)) = 0;
                    a_P(isnan(a_P)) = 0;
                    a_Wz(isnan(a_Wz)) = 0;
                    a_Wy1(isnan(a_Wy1)) = 0;
                    a_Wy2(isnan(a_Wy2)) = 0;
                    a_Wy3(isnan(a_Wy3)) = 0;
                    a_Wy4(isnan(a_Wy4)) = 0;
                    a_Wy5(isnan(a_Wy5)) = 0;
                    a_Wx(isnan(a_Wx)) = 0;
                A(n,[(n-(max_y*max_x)) (n-max_y) nWy1  nWy2 nWy3 nWy4 nWy5 n (n+(max_y*max_x))]) = [a_Wz a_Wx a_Wy1  a_Wy2 a_Wy3 a_Wy4 a_Wy5 a_P a_Ez];
                    
                elseif (i==max_y && j==1)
%                     B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) = zeros;
                    a_Ex = -wsize(1)*sigma(i,j,k)*sigma(i,j+1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(sigma(i,j+1,k)*0.5*(x(j+1)-x(j))+sigma(i,j,k)*0.5*(x(j+2)-x(j+1)));
                    a_Wy1 = (-sigma(i,j,k)*sigma(a(1),b(1),c(1))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(1),b(1),c(1))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    if c(2) == 0
                        a_Wy2 = 0;
                    else
                        a_Wy2 = (-sigma(i,j,k)*sigma(a(2),b(2),c(2))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(2),b(2),c(2))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    if c(3) == 0 
                        a_Wy3 = 0;
                    else
                    a_Wy3 = (-sigma(i,j,k)*sigma(a(3),b(3),c(3))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(3),b(3),c(3))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    if c(4) == 0 
                        a_Wy4 = 0;
                    else
                    a_Wy4 = (-sigma(i,j,k)*sigma(a(4),b(4),c(4))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(4),b(4),c(4))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    if c(5) == 0 
                        a_Wy5 = 0;
                    else
                    a_Wy5 = (-sigma(i,j,k)*sigma(a(5),b(5),c(5))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(5),b(5),c(5))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    a_Ez = -wsize(1)*sigma(i,j,k)*sigma(i,j,k+1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(sigma(i,j,k+1)*0.5*(z(k+1)-z(k))+sigma(i,j,k)*0.5*(z(k+2)-z(k+1)));
                    a_Wz = -wsize(1)*sigma(i,j,k)*sigma(i,j,k-1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(sigma(i,j,k)*0.5*(z(k)-z(k-1))+sigma(i,j,k-1)*0.5*(z(k+1)-z(k)));
                    a_P  =  -wsize(1)*(a_Ex +a_Wy1  +a_Wy2 +a_Wy3 +a_Wy4 +a_Wy5 + a_Ez + a_Wz);
                    a_Ex(isnan(a_Ex)) = 0;
                    a_Ez(isnan(a_Ez)) = 0;
                    a_P(isnan(a_P)) = 0;
                    a_Wz(isnan(a_Wz)) = 0;
                    a_Wy1(isnan(a_Wy1)) = 0;
                    a_Wy2(isnan(a_Wy2)) = 0;
                    a_Wy3(isnan(a_Wy3)) = 0;
                    a_Wy4(isnan(a_Wy4)) = 0;
                    a_Wy5(isnan(a_Wy5)) = 0;
                A(n,[(n-(max_y*max_x)) nWy1  nWy2 nWy3 nWy4 nWy5 n (n+max_y) (n+(max_y*max_x))]) = [a_Wz a_Wy1  a_Wy2 a_Wy3 a_Wy4 a_Wy5 a_P a_Ex a_Ez];
                    
                elseif (i==max_y && k==max_z)
%                     B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) = zeros;
                    a_Ex = -wsize(1)*sigma(i,j,k)*sigma(i,j+1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(sigma(i,j+1,k)*0.5*(x(j+1)-x(j))+sigma(i,j,k)*0.5*(x(j+2)-x(j+1)));
                    a_Wx = -wsize(1)*sigma(i,j,k)*sigma(i,j-1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(sigma(i,j,k)*0.5*(x(j)-x(j-1))+sigma(i,j-1,k)*0.5*(x(j+1)-x(j)));
                    a_Wy1 = (-sigma(i,j,k)*sigma(a(1),b(1),c(1))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(1),b(1),c(1))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    if c(2) == 0
                        a_Wy2 = 0;
                    else
                        a_Wy2 = (-sigma(i,j,k)*sigma(a(2),b(2),c(2))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(2),b(2),c(2))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    if c(3) == 0 
                        a_Wy3 = 0;
                    else
                    a_Wy3 = (-sigma(i,j,k)*sigma(a(3),b(3),c(3))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(3),b(3),c(3))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    if c(4) == 0 
                        a_Wy4 = 0;
                    else
                    a_Wy4 = (-sigma(i,j,k)*sigma(a(4),b(4),c(4))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(4),b(4),c(4))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    if c(5) == 0 
                        a_Wy5 = 0;
                    else
                    a_Wy5 = (-sigma(i,j,k)*sigma(a(5),b(5),c(5))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(5),b(5),c(5))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    a_Wz = -wsize(1)*sigma(i,j,k)*sigma(i,j,k-1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(sigma(i,j,k)*0.5*(z(k)-z(k-1))+sigma(i,j,k-1)*0.5*(z(k+1)-z(k)));
                    a_P  =  -wsize(1)*(a_Ex + a_Wx +a_Wy1  +a_Wy2 +a_Wy3 +a_Wy4 +a_Wy5 + a_Wz);
                    a_Ex(isnan(a_Ex)) = 0;
                    a_P(isnan(a_P)) = 0;
                    a_Wz(isnan(a_Wz)) = 0;
                    a_Wy1(isnan(a_Wy1)) = 0;
                    a_Wy2(isnan(a_Wy2)) = 0;
                    a_Wy3(isnan(a_Wy3)) = 0;
                    a_Wy4(isnan(a_Wy4)) = 0;
                    a_Wy5(isnan(a_Wy5)) = 0;
                    a_Wx(isnan(a_Wx)) = 0;
                A(n,[(n-(max_y*max_x)) (n-max_y) nWy1  nWy2 nWy3 nWy4 nWy5 n (n+max_y)]) = [a_Wz a_Wx a_Wy1  a_Wy2 a_Wy3 a_Wy4 a_Wy5 a_P a_Ex];
                    
                    
                elseif (j==max_x && k==1)
%                     B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) = zeros;
                    a_Wx = -wsize(1)*sigma(i,j,k)*sigma(i,j-1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(sigma(i,j,k)*0.5*(x(j)-x(j-1))+sigma(i,j-1,k)*0.5*(x(j+1)-x(j)));
                    a_Ey = -wsize(1)*sigma(i,j,k)*sigma(i+1,j,k)*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(i+1,j,k)*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1)));
                    a_Wy1 = (-sigma(i,j,k)*sigma(a(1),b(1),c(1))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(1),b(1),c(1))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    if c(2) == 0
                        a_Wy2 = 0;
                    else
                        a_Wy2 = (-sigma(i,j,k)*sigma(a(2),b(2),c(2))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(2),b(2),c(2))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    if c(3) == 0 
                        a_Wy3 = 0;
                    else
                    a_Wy3 = (-sigma(i,j,k)*sigma(a(3),b(3),c(3))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(3),b(3),c(3))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    if c(4) == 0 
                        a_Wy4 = 0;
                    else
                    a_Wy4 = (-sigma(i,j,k)*sigma(a(4),b(4),c(4))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(4),b(4),c(4))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    if c(5) == 0 
                        a_Wy5 = 0;
                    else
                    a_Wy5 = (-sigma(i,j,k)*sigma(a(5),b(5),c(5))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(5),b(5),c(5))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    a_Ez = -wsize(1)*sigma(i,j,k)*sigma(i,j,k+1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(sigma(i,j,k+1)*0.5*(z(k+1)-z(k))+sigma(i,j,k)*0.5*(z(k+2)-z(k+1)));
                    a_P  =  -wsize(1)*(a_Wx + a_Ey +a_Wy1  +a_Wy2 +a_Wy3 +a_Wy4 +a_Wy5 + a_Ez);
                    a_Ey(isnan(a_Ey)) = 0;
                    a_Ez(isnan(a_Ez)) = 0;
                    a_P(isnan(a_P)) = 0;
                    a_Wy1(isnan(a_Wy1)) = 0;
                    a_Wy2(isnan(a_Wy2)) = 0;
                    a_Wy3(isnan(a_Wy3)) = 0;
                    a_Wy4(isnan(a_Wy4)) = 0;
                    a_Wy5(isnan(a_Wy5)) = 0;
                    a_Wx(isnan(a_Wx)) = 0;
                A(n,[(n-max_y) nWy1  nWy2 nWy3 nWy4 nWy5 n (n+1) (n+(max_y*max_x))]) = [a_Wx a_Wy1  a_Wy2 a_Wy3 a_Wy4 a_Wy5 a_P a_Ey a_Ez];
                    
                elseif (j==max_x && k==max_z)
%                     B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) = zeros;
                    a_Wx = -wsize(1)*sigma(i,j,k)*sigma(i,j-1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(sigma(i,j,k)*0.5*(x(j)-x(j-1))+sigma(i,j-1,k)*0.5*(x(j+1)-x(j)));
                    a_Ey = -wsize(1)*sigma(i,j,k)*sigma(i+1,j,k)*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(i+1,j,k)*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1)));
                    a_Wy1 = (-sigma(i,j,k)*sigma(a(1),b(1),c(1))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(1),b(1),c(1))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    if c(2) == 0
                        a_Wy2 = 0;
                    else
                        a_Wy2 = (-sigma(i,j,k)*sigma(a(2),b(2),c(2))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(2),b(2),c(2))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    if c(3) == 0 
                        a_Wy3 = 0;
                    else
                    a_Wy3 = (-sigma(i,j,k)*sigma(a(3),b(3),c(3))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(3),b(3),c(3))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    if c(4) == 0 
                        a_Wy4 = 0;
                    else
                    a_Wy4 = (-sigma(i,j,k)*sigma(a(4),b(4),c(4))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(4),b(4),c(4))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    if c(5) == 0 
                        a_Wy5 = 0;
                    else
                    a_Wy5 = (-sigma(i,j,k)*sigma(a(5),b(5),c(5))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(5),b(5),c(5))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    a_Wz = -wsize(1)*sigma(i,j,k)*sigma(i,j,k-1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(sigma(i,j,k)*0.5*(z(k)-z(k-1))+sigma(i,j,k-1)*0.5*(z(k+1)-z(k)));
                    a_P  =  -wsize(1)*(a_Wx + a_Ey +a_Wy1  +a_Wy2 +a_Wy3 +a_Wy4 +a_Wy5 + a_Wz);
                    a_Ey(isnan(a_Ey)) = 0;
                    a_P(isnan(a_P)) = 0;
                    a_Wz(isnan(a_Wz)) = 0;
                    a_Wy1(isnan(a_Wy1)) = 0;
                    a_Wy2(isnan(a_Wy2)) = 0;
                    a_Wy3(isnan(a_Wy3)) = 0;
                    a_Wy4(isnan(a_Wy4)) = 0;
                    a_Wy5(isnan(a_Wy5)) = 0;
                    a_Wx(isnan(a_Wx)) = 0;
                A(n,[(n-(max_y*max_x)) (n-max_y) nWy1  nWy2 nWy3 nWy4 nWy5 n (n+1)]) = [a_Wz a_Wx a_Wy1  a_Wy2 a_Wy3 a_Wy4 a_Wy5 a_P a_Ey];
                    
                elseif (j==1 && k==max_z)
%                     B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) = zeros;
                    a_Ex = -wsize(1)*sigma(i,j,k)*sigma(i,j+1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(sigma(i,j+1,k)*0.5*(x(j+1)-x(j))+sigma(i,j,k)*0.5*(x(j+2)-x(j+1)));
                    a_Ey = -wsize(1)*sigma(i,j,k)*sigma(i+1,j,k)*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(i+1,j,k)*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1)));
                    a_Wy1 = (-sigma(i,j,k)*sigma(a(1),b(1),c(1))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(1),b(1),c(1))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    if c(2) == 0
                        a_Wy2 = 0;
                    else
                        a_Wy2 = (-sigma(i,j,k)*sigma(a(2),b(2),c(2))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(2),b(2),c(2))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    if c(3) == 0 
                        a_Wy3 = 0;
                    else
                    a_Wy3 = (-sigma(i,j,k)*sigma(a(3),b(3),c(3))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(3),b(3),c(3))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    if c(4) == 0 
                        a_Wy4 = 0;
                    else
                    a_Wy4 = (-sigma(i,j,k)*sigma(a(4),b(4),c(4))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(4),b(4),c(4))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    if c(5) == 0 
                        a_Wy5 = 0;
                    else
                    a_Wy5 = (-sigma(i,j,k)*sigma(a(5),b(5),c(5))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(5),b(5),c(5))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    a_Wz = -wsize(1)*sigma(i,j,k)*sigma(i,j,k-1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(sigma(i,j,k)*0.5*(z(k)-z(k-1))+sigma(i,j,k-1)*0.5*(z(k+1)-z(k)));
                    a_P  =  -wsize(1)*(a_Ex + a_Ey +a_Wy1  +a_Wy2 +a_Wy3 +a_Wy4 +a_Wy5 + a_Wz);
                    a_Ex(isnan(a_Ex)) = 0;
                    a_Ey(isnan(a_Ey)) = 0;
                    a_P(isnan(a_P)) = 0;
                    a_Wz(isnan(a_Wz)) = 0;
                    a_Wy1(isnan(a_Wy1)) = 0;
                    a_Wy2(isnan(a_Wy2)) = 0;
                    a_Wy3(isnan(a_Wy3)) = 0;
                    a_Wy4(isnan(a_Wy4)) = 0;
                    a_Wy5(isnan(a_Wy5)) = 0;
                A(n,[(n-(max_y*max_x)) nWy1  nWy2 nWy3 nWy4 nWy5 n (n+1) (n+max_y)]) = [a_Wz a_Wy1  a_Wy2 a_Wy3 a_Wy4 a_Wy5 a_P a_Ey a_Ex];
                    
                elseif (j==1 && k==1)
%                     B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) = zeros;
                    a_Ex = -wsize(1)*sigma(i,j,k)*sigma(i,j+1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(sigma(i,j+1,k)*0.5*(x(j+1)-x(j))+sigma(i,j,k)*0.5*(x(j+2)-x(j+1)));
                    a_Ey = -wsize(1)*sigma(i,j,k)*sigma(i+1,j,k)*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(i+1,j,k)*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1)));
                    a_Wy1 = (-sigma(i,j,k)*sigma(a(1),b(1),c(1))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(1),b(1),c(1))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    if c(2) == 0
                        a_Wy2 = 0;
                    else
                        a_Wy2 = (-sigma(i,j,k)*sigma(a(2),b(2),c(2))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(2),b(2),c(2))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    if c(3) == 0 
                        a_Wy3 = 0;
                    else
                    a_Wy3 = (-sigma(i,j,k)*sigma(a(3),b(3),c(3))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(3),b(3),c(3))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    if c(4) == 0 
                        a_Wy4 = 0;
                    else
                    a_Wy4 = (-sigma(i,j,k)*sigma(a(4),b(4),c(4))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(4),b(4),c(4))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    if c(5) == 0 
                        a_Wy5 = 0;
                    else
                    a_Wy5 = (-sigma(i,j,k)*sigma(a(5),b(5),c(5))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(5),b(5),c(5))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    a_Ez = -wsize(1)*sigma(i,j,k)*sigma(i,j,k+1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(sigma(i,j,k+1)*0.5*(z(k+1)-z(k))+sigma(i,j,k)*0.5*(z(k+2)-z(k+1)));
                    a_P  =  -wsize(1)*(a_Ex + a_Ey +a_Wy1  +a_Wy2 +a_Wy3 +a_Wy4 +a_Wy5 + a_Ez);
                    a_Ex(isnan(a_Ex)) = 0;
                    a_Ey(isnan(a_Ey)) = 0;
                    a_Ez(isnan(a_Ez)) = 0;
                    a_P(isnan(a_P)) = 0;
                    a_Wy1(isnan(a_Wy1)) = 0;
                    a_Wy2(isnan(a_Wy2)) = 0;
                    a_Wy3(isnan(a_Wy3)) = 0;
                    a_Wy4(isnan(a_Ey4)) = 0;
                    a_Wy5(isnan(a_Ey5)) = 0;
                A(n,[nWy1  nWy2 nWy3 nWy4 nWy5 n (n+1) (n+max_y) (n+(max_y*max_x))]) = [a_Wy1  a_Wy2 a_Wy3 a_Wy4 a_Wy5 a_P a_Ey a_Ex a_Ez];
                    
                    %----------------this part defines the surfaces
                    %defining the surface boundaries---------------
                elseif i==1
%                     B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) = zeros;
                    a_Ex = -sigma(i,j,k)*sigma(i,j+1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(sigma(i,j+1,k)*0.5*(x(j+1)-x(j))+sigma(i,j,k)*0.5*(x(j+2)-x(j+1)));
                    a_Wx = -sigma(i,j,k)*sigma(i,j-1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(sigma(i,j,k)*0.5*(x(j)-x(j-1))+sigma(i,j-1,k)*0.5*(x(j+1)-x(j)));
                    a_Ey = -sigma(i,j,k)*sigma(i+1,j,k)*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(i+1,j,k)*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1)));
                    a_Ez = -sigma(i,j,k)*sigma(i,j,k+1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(sigma(i,j,k+1)*0.5*(z(k+1)-z(k))+sigma(i,j,k)*0.5*(z(k+2)-z(k+1)));
                    a_Wz = -sigma(i,j,k)*sigma(i,j,k-1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(sigma(i,j,k)*0.5*(z(k)-z(k-1))+sigma(i,j,k-1)*0.5*(z(k+1)-z(k)));
                    a_P  =  -(a_Ex + a_Wx + a_Ey + a_Ez + a_Wz);
                    a_Ex(isnan(a_Ex)) = 0;
                    a_Ey(isnan(a_Ey)) = 0;
                    a_Ez(isnan(a_Ez)) = 0;
                    a_P(isnan(a_P)) = 0;
                    a_Wz(isnan(a_Wz)) = 0;
                    a_Wx(isnan(a_Wx)) = 0;
                A(n,[(n-(max_y*max_x)) (n-max_y) n (n+1) (n+max_y) (n+(max_y*max_x))]) = [a_Wz a_Wx a_P a_Ey a_Ex a_Ez];
                    
                elseif i==max_y
%                     B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) = zeros;
                    a_Ex = -wsize(1)*sigma(i,j,k)*sigma(i,j+1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(sigma(i,j+1,k)*0.5*(x(j+1)-x(j))+sigma(i,j,k)*0.5*(x(j+2)-x(j+1)));
                    a_Wx = -wsize(1)*sigma(i,j,k)*sigma(i,j-1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(sigma(i,j,k)*0.5*(x(j)-x(j-1))+sigma(i,j-1,k)*0.5*(x(j+1)-x(j)));
                    a_Wy1 = (-sigma(i,j,k)*sigma(a(1),b(1),c(1))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(1),b(1),c(1))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    if c(2) == 0
                        a_Wy2 = 0;
                    else
                        a_Wy2 = (-sigma(i,j,k)*sigma(a(2),b(2),c(2))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(2),b(2),c(2))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    if c(3) == 0 
                        a_Wy3 = 0;
                    else
                    a_Wy3 = (-sigma(i,j,k)*sigma(a(3),b(3),c(3))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(3),b(3),c(3))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    if c(4) == 0 
                        a_Wy4 = 0;
                    else
                    a_Wy4 = (-sigma(i,j,k)*sigma(a(4),b(4),c(4))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(4),b(4),c(4))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    if c(5) == 0 
                        a_Wy5 = 0;
                    else
                    a_Wy5 = (-sigma(i,j,k)*sigma(a(5),b(5),c(5))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(5),b(5),c(5))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    a_Ez = -wsize(1)*sigma(i,j,k)*sigma(i,j,k+1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(sigma(i,j,k+1)*0.5*(z(k+1)-z(k))+sigma(i,j,k)*0.5*(z(k+2)-z(k+1)));
                    a_Wz = -wsize(1)*sigma(i,j,k)*sigma(i,j,k-1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(sigma(i,j,k)*0.5*(z(k)-z(k-1))+sigma(i,j,k-1)*0.5*(z(k+1)-z(k)));
                    a_P  =  -wsize(1)*(a_Ex + a_Wx +a_Wy1  +a_Wy2 +a_Wy3 +a_Wy4 +a_Wy5 + a_Ez + a_Wz);
                    a_Ex(isnan(a_Ex)) = 0;
                    a_Ez(isnan(a_Ez)) = 0;
                    a_P(isnan(a_P)) = 0;
                    a_Wz(isnan(a_Wz)) = 0;
                    a_Wy1(isnan(a_Wy1)) = 0;
                    a_Wy2(isnan(a_Wy2)) = 0;
                    a_Wy3(isnan(a_Wy3)) = 0;
                    a_Wy4(isnan(a_Wy4)) = 0;
                    a_Wy5(isnan(a_Wy5)) = 0;
                    a_Wx(isnan(a_Wx)) = 0;
                A(n,[(n-(max_y*max_x)) (n-max_y) nWy1  nWy2 nWy3 nWy4 nWy5 n (n+max_y) (n+(max_y*max_x))]) = [a_Wz a_Wx a_Wy1  a_Wy2 a_Wy3 a_Wy4 a_Wy5 a_P a_Ex a_Ez];
                    
                elseif j==1
%                     B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) = zeros;
                    a_Ex = -wsize(1)*sigma(i,j,k)*sigma(i,j+1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(sigma(i,j+1,k)*0.5*(x(j+1)-x(j))+sigma(i,j,k)*0.5*(x(j+2)-x(j+1)));
                    a_Ey = -wsize(1)*sigma(i,j,k)*sigma(i+1,j,k)*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(i+1,j,k)*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1)));
                    a_Wy1 = (-sigma(i,j,k)*sigma(a(1),b(1),c(1))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(1),b(1),c(1))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    if c(2) == 0
                        a_Wy2 = 0;
                    else
                        a_Wy2 = (-sigma(i,j,k)*sigma(a(2),b(2),c(2))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(2),b(2),c(2))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    if c(3) == 0 
                        a_Wy3 = 0;
                    else
                    a_Wy3 = (-sigma(i,j,k)*sigma(a(3),b(3),c(3))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(3),b(3),c(3))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    if c(4) == 0 
                        a_Wy4 = 0;
                    else
                    a_Wy4 = (-sigma(i,j,k)*sigma(a(4),b(4),c(4))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(4),b(4),c(4))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    if c(5) == 0 
                        a_Wy5 = 0;
                    else
                    a_Wy5 = (-sigma(i,j,k)*sigma(a(5),b(5),c(5))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(5),b(5),c(5))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    a_Ez = -wsize(1)*sigma(i,j,k)*sigma(i,j,k+1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(sigma(i,j,k+1)*0.5*(z(k+1)-z(k))+sigma(i,j,k)*0.5*(z(k+2)-z(k+1)));
                    a_Wz = -wsize(1)*sigma(i,j,k)*sigma(i,j,k-1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(sigma(i,j,k)*0.5*(z(k)-z(k-1))+sigma(i,j,k-1)*0.5*(z(k+1)-z(k)));
                    a_P  =  -wsize(1)*(a_Ex + a_Ey +a_Wy1  +a_Wy2 +a_Wy3 +a_Wy4 +a_Wy5 + a_Ez + a_Wz);
                    a_Ex(isnan(a_Ex)) = 0;
                    a_Ey(isnan(a_Ey)) = 0;
                    a_Ez(isnan(a_Ez)) = 0;
                    a_P(isnan(a_P)) = 0;
                    a_Wz(isnan(a_Wz)) = 0;
                    a_Wy1(isnan(a_Wy1)) = 0;
                    a_Wy2(isnan(a_Wy2)) = 0;
                    a_Wy3(isnan(a_Wy3)) = 0;
                    a_Wy4(isnan(a_Wy4)) = 0;
                    a_Wy5(isnan(a_Wy5)) = 0;
                A(n,[(n-(max_y*max_x)) nWy1  nWy2 nWy3 nWy4 nWy5 n (n+1) (n+max_y) (n+(max_y*max_x))]) = [a_Wz a_Wy1  a_Wy2 a_Wy3 a_Wy4 a_Wy5 a_P a_Ey a_Ex a_Ez];
                    
                elseif j==max_x
%                     B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) = zeros;
                    a_Wx = -wsize(1)*sigma(i,j,k)*sigma(i,j-1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(sigma(i,j,k)*0.5*(x(j)-x(j-1))+sigma(i,j-1,k)*0.5*(x(j+1)-x(j)));
                    a_Ey = -wsize(1)*sigma(i,j,k)*sigma(i+1,j,k)*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(i+1,j,k)*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1)));
                    a_Wy1 = (-sigma(i,j,k)*sigma(a(1),b(1),c(1))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(1),b(1),c(1))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    if c(2) == 0
                        a_Wy2 = 0;
                    else
                        a_Wy2 = (-sigma(i,j,k)*sigma(a(2),b(2),c(2))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(2),b(2),c(2))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    if c(3) == 0 
                        a_Wy3 = 0;
                    else
                    a_Wy3 = (-sigma(i,j,k)*sigma(a(3),b(3),c(3))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(3),b(3),c(3))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    if c(4) == 0 
                        a_Wy4 = 0;
                    else
                    a_Wy4 = (-sigma(i,j,k)*sigma(a(4),b(4),c(4))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(4),b(4),c(4))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    if c(5) == 0 
                        a_Wy5 = 0;
                    else
                    a_Wy5 = (-sigma(i,j,k)*sigma(a(5),b(5),c(5))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(5),b(5),c(5))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    a_Ez = -wsize(1)*sigma(i,j,k)*sigma(i,j,k+1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(sigma(i,j,k+1)*0.5*(z(k+1)-z(k))+sigma(i,j,k)*0.5*(z(k+2)-z(k+1)));
                    a_Wz = -wsize(1)*sigma(i,j,k)*sigma(i,j,k-1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(sigma(i,j,k)*0.5*(z(k)-z(k-1))+sigma(i,j,k-1)*0.5*(z(k+1)-z(k)));
                    a_P  =  -wsize(1)*(a_Wx + a_Ey +a_Wy1  +a_Wy2 +a_Wy3 +a_Wy4 +a_Wy5 + a_Ez + a_Wz);
                    a_Ey(isnan(a_Ey)) = 0;
                    a_Ez(isnan(a_Ez)) = 0;
                    a_P(isnan(a_P)) = 0;
                    a_Wz(isnan(a_Wz)) = 0;
                    a_Wy1(isnan(a_Wy1)) = 0;
                    a_Wy2(isnan(a_Wy2)) = 0;
                    a_Wy3(isnan(a_Wy3)) = 0;
                    a_Wy4(isnan(a_Wy4)) = 0;
                    a_Wy5(isnan(a_Wy5)) = 0;
                    a_Wx(isnan(a_Wx)) = 0;
                A(n,[(n-(max_y*max_x)) (n-max_y) nWy1  nWy2 nWy3 nWy4 nWy5 n (n+1) (n+(max_y*max_x))]) = [a_Wz a_Wx a_Wy1  a_Wy2 a_Wy3 a_Wy4 a_Wy5 a_P a_Ey a_Ez];
                    
                elseif k==1
%                     B(n,[(n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) = zeros;
                    a_Ex = -wsize(1)*sigma(i,j,k)*sigma(i,j+1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(sigma(i,j+1,k)*0.5*(x(j+1)-x(j))+sigma(i,j,k)*0.5*(x(j+2)-x(j+1)));
                    a_Wx = -wsize(1)*sigma(i,j,k)*sigma(i,j-1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(sigma(i,j,k)*0.5*(x(j)-x(j-1))+sigma(i,j-1,k)*0.5*(x(j+1)-x(j)));
                    a_Ey = -wsize(1)*sigma(i,j,k)*sigma(i+1,j,k)*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(i+1,j,k)*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1)));
                    a_Wy1 = (-sigma(i,j,k)*sigma(a(1),b(1),c(1))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(1),b(1),c(1))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    if c(2) == 0
                        a_Wy2 = 0;
                    else
                        a_Wy2 = (-sigma(i,j,k)*sigma(a(2),b(2),c(2))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(2),b(2),c(2))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    if c(3) == 0 
                        a_Wy3 = 0;
                    else
                    a_Wy3 = (-sigma(i,j,k)*sigma(a(3),b(3),c(3))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(3),b(3),c(3))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    if c(4) == 0 
                        a_Wy4 = 0;
                    else
                    a_Wy4 = (-sigma(i,j,k)*sigma(a(4),b(4),c(4))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(4),b(4),c(4))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    if c(5) == 0 
                        a_Wy5 = 0;
                    else
                    a_Wy5 = (-sigma(i,j,k)*sigma(a(5),b(5),c(5))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(5),b(5),c(5))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    a_Ez = -wsize(1)*sigma(i,j,k)*sigma(i,j,k+1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(sigma(i,j,k+1)*0.5*(z(k+1)-z(k))+sigma(i,j,k)*0.5*(z(k+2)-z(k+1)));
                    a_P  =  -wsize(1)*(a_Ex + a_Wx + a_Ey +a_Wy1  +a_Wy2 +a_Wy3 +a_Wy4 +a_Wy5 + a_Ez);
                    a_Ex(isnan(a_Ex)) = 0;
                    a_Ey(isnan(a_Ey)) = 0;
                    a_Ez(isnan(a_Ez)) = 0;
                    a_P(isnan(a_P)) = 0;
                    a_Wy1(isnan(a_Wy1)) = 0;
                    a_Wy2(isnan(a_Wy2)) = 0;
                    a_Wy3(isnan(a_Wy3)) = 0;
                    a_Wy4(isnan(a_Wy4)) = 0;
                    a_Wy5(isnan(a_Wy5)) = 0;
                    a_Wx(isnan(a_Wx)) = 0;
                A(n,[(n-max_y) nWy1  nWy2 nWy3 nWy4 nWy5 n (n+1) (n+max_y) (n+(max_y*max_x))]) = [a_Wx a_Wy1  a_Wy2 a_Wy3 a_Wy4 a_Wy5 a_P a_Ey a_Ex a_Ez];
                
                elseif k==max_z
%                     B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) = zeros;
                    a_Ex = -wsize(1)*sigma(i,j,k)*sigma(i,j+1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(sigma(i,j+1,k)*0.5*(x(j+1)-x(j))+sigma(i,j,k)*0.5*(x(j+2)-x(j+1)));
                    a_Wx = -wsize(1)*sigma(i,j,k)*sigma(i,j-1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(sigma(i,j,k)*0.5*(x(j)-x(j-1))+sigma(i,j-1,k)*0.5*(x(j+1)-x(j)));
                    a_Ey = -wsize(1)*sigma(i,j,k)*sigma(i+1,j,k)*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(i+1,j,k)*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1)));
                    a_Wy1 = (-sigma(i,j,k)*sigma(a(1),b(1),c(1))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(1),b(1),c(1))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    if c(2) == 0
                        a_Wy2 = 0;
                    else
                        a_Wy2 = (-sigma(i,j,k)*sigma(a(2),b(2),c(2))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(2),b(2),c(2))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    if c(3) == 0 
                        a_Wy3 = 0;
                    else
                    a_Wy3 = (-sigma(i,j,k)*sigma(a(3),b(3),c(3))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(3),b(3),c(3))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    if c(4) == 0 
                        a_Wy4 = 0;
                    else
                    a_Wy4 = (-sigma(i,j,k)*sigma(a(4),b(4),c(4))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(4),b(4),c(4))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    if c(5) == 0 
                        a_Wy5 = 0;
                    else
                    a_Wy5 = (-sigma(i,j,k)*sigma(a(5),b(5),c(5))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(5),b(5),c(5))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    a_Wz = -wsize(1)*sigma(i,j,k)*sigma(i,j,k-1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(sigma(i,j,k)*0.5*(z(k)-z(k-1))+sigma(i,j,k-1)*0.5*(z(k+1)-z(k)));
                    a_P  =  -wsize(1)*(a_Ex + a_Wx + a_Ey +a_Wy1  +a_Wy2 +a_Wy3 +a_Wy4 +a_Wy5 + a_Wz);
                    a_Ex(isnan(a_Ex)) = 0;
                    a_Ey(isnan(a_Ey)) = 0;
                    a_P(isnan(a_P)) = 0;
                    a_Wz(isnan(a_Wz)) = 0;
                    a_Wy1(isnan(a_Wy1)) = 0;
                    a_Wy2(isnan(a_Wy2)) = 0;
                    a_Wy3(isnan(a_Wy3)) = 0;
                    a_Wy4(isnan(a_Wy4)) = 0;
                    a_Wy5(isnan(a_Wy5)) = 0;
                    a_Wx(isnan(a_Wx)) = 0;
                A(n,[(n-(max_y*max_x)) (n-max_y) nWy1  nWy2 nWy3 nWy4 nWy5 n (n+1) (n+max_y)]) = [a_Wz a_Wx a_Wy1  a_Wy2 a_Wy3 a_Wy4 a_Wy5 a_P a_Ey a_Ex];
                
               elseif (i>1 && j>1 && k>1) && (i<max_y && j<max_x && k<max_z)
%                    B(n,[(n-(max_y*max_x)) (n-max_y) (n-1) n (n+1) (n+max_y) (n+(max_y*max_x))]) = zeros;
                    a_Ex = -wsize(1)*sigma(i,j,k)*sigma(i,j+1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(sigma(i,j+1,k)*0.5*(x(j+1)-x(j))+sigma(i,j,k)*0.5*(x(j+2)-x(j+1)));
                    a_Wx = -wsize(1)*sigma(i,j,k)*sigma(i,j-1,k)*(y(i+1)-y(i))*(z(k+1)-z(k))/(sigma(i,j,k)*0.5*(x(j)-x(j-1))+sigma(i,j-1,k)*0.5*(x(j+1)-x(j)));
                    a_Ey = -wsize(1)*sigma(i,j,k)*sigma(i+1,j,k)*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(i+1,j,k)*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1)));
                    a_Wy1 = (-sigma(i,j,k)*sigma(a(1),b(1),c(1))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(1),b(1),c(1))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    if c(2) == 0
                        a_Wy2 = 0;
                    else
                        a_Wy2 = (-sigma(i,j,k)*sigma(a(2),b(2),c(2))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(2),b(2),c(2))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    if c(3) == 0 
                        a_Wy3 = 0;
                    else
                    a_Wy3 = (-sigma(i,j,k)*sigma(a(3),b(3),c(3))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(3),b(3),c(3))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    if c(4) == 0 
                        a_Wy4 = 0;
                    else
                    a_Wy4 = (-sigma(i,j,k)*sigma(a(4),b(4),c(4))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(4),b(4),c(4))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    if c(5) == 0 
                        a_Wy5 = 0;
                    else
                    a_Wy5 = (-sigma(i,j,k)*sigma(a(5),b(5),c(5))*(x(j+1)-x(j))*(z(k+1)-z(k))/(sigma(a(5),b(5),c(5))*0.5*(y(i+1)-y(i))+sigma(i,j,k)*0.5*(y(i+2)-y(i+1))))/1;
                    end
                    a_Ez = -wsize(1)*sigma(i,j,k)*sigma(i,j,k+1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(sigma(i,j,k+1)*0.5*(z(k+1)-z(k))+sigma(i,j,k)*0.5*(z(k+2)-z(k+1)));
                    a_Wz = -wsize(1)*sigma(i,j,k)*sigma(i,j,k-1)*(x(j+1)-x(j))*(y(i+1)-y(i))/(sigma(i,j,k)*0.5*(z(k)-z(k-1))+sigma(i,j,k-1)*0.5*(z(k+1)-z(k)));
                    a_P  =  -wsize(1)*(a_Ex + a_Wx + a_Ey +a_Wy1  +a_Wy2 +a_Wy3 +a_Wy4 +a_Wy5 + a_Ez + a_Wz);
                    a_Ex(isnan(a_Ex)) = 0;
                    a_Ey(isnan(a_Ey)) = 0;
                    a_Ez(isnan(a_Ez)) = 0;
                    a_P(isnan(a_P)) = 0;
                    a_Wz(isnan(a_Wz)) = 0;
                    a_Wy1(isnan(a_Wy1)) = 0;
                    a_Wy2(isnan(a_Wy2)) = 0;
                    a_Wy3(isnan(a_Wy3)) = 0;
                    a_Wy4(isnan(a_Wy4)) = 0;
                    a_Wy5(isnan(a_Wy5)) = 0;
                    a_Wx(isnan(a_Wx)) = 0;
                    
                A(n,[(n-(max_y*max_x)) (n-max_y) nWy1  nWy2 nWy3 nWy4 nWy5 n (n+1) (n+max_y) (n+(max_y*max_x))]) = [a_Wz a_Wx a_Wy1  a_Wy2 a_Wy3 a_Wy4 a_Wy5 a_P a_Ey a_Ex a_Ez];
                    
                end           
    
    
end
