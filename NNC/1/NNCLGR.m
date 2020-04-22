[max_ya, max_xa, max_za]=size(Pn);
max_y = 48;
max_x = 139;
max_z = 9;

diffy = max_ya - max_y;
diffx = max_xa - max_x;
diffz = max_za - max_z - 10;

FL = max_y * max_x;

for n = 1:length(NNCG)
      
    if NNCG(n) > FL
    
        k = 1 + floor(NNCG(n)/FL)+ diffz;
    
        i = 1 + floor(rem(NNCG(n),(FL*(k-1-diffz)))/max_x)+diffy;
        if (rem(NNCG(n),(FL*(k-1)))/max_x) - floor(rem(NNCG(n),(FL*(k-1)))/max_x) == 0
            i = i - 1 + diffy;
        end
            
        j = rem(NNCG(n),max_y*max_x)-(max_x*(i-1-diffy));
        if j == 0
            j = max_x;
        end
    
    else
    
        k = 1 + diffz;
   
        i = 1 + floor(NNCG(n)/max_x) + diffy;
        if (NNCG(n)/max_x) - floor((NNCG(n)/max_x)) == 0
            i = i - 1 + diffy;
        end
   
        j = NNCG(n) - (max_x*(i-1-diffy));   
        if j == 0
            j = max_x;
        end
    
    end
    
    NNCCG(n,:) = [ i j k ];

end

for n = 1:length(HOSTNUM)
      
    if HOSTNUM(n) > FL
    
        k = 1 + floor(HOSTNUM(n)/FL)+ diffz;
    
        i = 1 + floor(rem(HOSTNUM(n),(FL*(k-1-diffz)))/max_x)+diffy;
        if (rem(HOSTNUM(n),(FL*(k-1)))/max_x) - floor(rem(HOSTNUM(n),(FL*(k-1)))/max_x) == 0
            i = i - 1 + diffy;
        end
            
        j = rem(HOSTNUM(n),max_y*max_x)-(max_x*(i-1-diffy));
        if j == 0
            j = max_x;
        end
    
    else
    
        k = 1 + diffz;
   
        i = 1 + floor(HOSTNUM(n)/max_x) + diffy;
        if (HOSTNUM(n)/max_x) - floor((HOSTNUM(n)/max_x)) == 0
            i = i - 1 + diffy;
        end
   
        j = HOSTNUM(n) - (max_x*(i-1-diffy));   
        if j == 0
            j = max_x;
        end
    
    end
    
    HOSTNUMi(n,:) = [ i j k ];

end


[max_y, max_x, max_z]=size(sigma);

FL = max_y * max_x;

for n = 1:length(NNCL)
      
    if NNCL(n) > FL
    
        k = 1 + floor(NNCL(n)/FL);
        if (rem(NNCL(n),(FL))) == 0
            k = k - 1;
        end
    
        i = 1 + floor((NNCL(n)-(FL*(k-1)))/max_x);
        if rem(NNCL(n),max_x) == 0;
            i = i - 1;
        end
    
        j = floor(NNCL(n)-(FL*(k-1)))-(max_x*(i-1));
        if j == 0
            j = max_x;
        end
    
    else
    
        k = 1;
   
        i = 1 + floor(NNCL(n)/max_x);
        if (rem(NNCL(n),max_x)) == 0
            i = i - 1;
        end
           
        j = NNCL(n) - (max_x*(i-1)); 
        if j == 0
            j = max_x;
        end
    
    end
    
    NNCCL(n,:) = [ i j k ];

end