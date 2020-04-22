[max_ya, max_xa, max_za]=size(sigma);
max_y = DY;
max_x = DX;
max_z = DZ;

diffy = max_ya - max_y;
diffx = max_xa - max_x;
diffz = max_za - max_z - 10;

FL = max_y * max_x;
NNC1 = NNC.NNC1;
NNC2 = NNC.NNC2;

for n = 1:length(NNC1)
      
    if NNC1(n) > FL
    
        k = 1 + floor(NNC1(n)/FL)+ diffz;
    
        i = 1 + floor(rem(NNC1(n),(FL*(k-1-diffz)))/max_x)+diffy;
        if (rem(NNC1(n),(FL*(k-1)))/max_x) - floor(rem(NNC1(n),(FL*(k-1)))/max_x) == 0
            i = i - 1 + diffy;
        end
            
        j = rem(NNC1(n),max_y*max_x)-(max_x*(i-1-diffy));
        if j == 0
            j = max_x;
        end
    
    else
    
        k = 1 + diffz;
   
        i = 1 + floor(NNC1(n)/max_x) + diffy;
        if (NNC1(n)/max_x) - floor((NNC1(n)/max_x)) == 0
            i = i - 1 + diffy;
        end
   
        j = NNC1(n) - (max_x*(i-1-diffy));   
        if j == 0
            j = max_x;
        end
    
    end
    
    NNCC1(n,:) = [ i j k ];

end

for n = 1:length(NNC2)
      
    if NNC2(n) > FL
    
        k = 1 + floor(NNC2(n)/FL)+ diffz;
    
        i = 1 + floor(rem(NNC2(n),(FL*(k-1-diffz)))/max_x)+diffy;
        if (rem(NNC2(n),(FL*(k-1)))/max_x) - floor(rem(NNC2(n),(FL*(k-1)))/max_x) == 0
            i = i - 1 + diffy;
        end
            
        j = rem(NNC2(n),max_y*max_x)-(max_x*(i-1-diffy));
        if j == 0
            j = max_x;
        end
    
    else
    
        k = 1 + diffz;
   
        i = 1 + floor(NNC2(n)/max_x) + diffy;
        if (NNC2(n)/max_x) - floor((NNC2(n)/max_x)) == 0
            i = i - 1 + diffy;
        end
   
        j = NNC2(n) - (max_x*(i-1-diffy));   
        if j == 0
            j = max_x;
        end
    
    end
    
    NNCC2(n,:) = [ i j k ];

end