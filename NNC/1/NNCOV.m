[max_y, max_x, max_z]=size(sigma);

FL = max_y * max_x;


for n = 1:length(NNCI1a)
      
    if NNCI1a(n) > FL
    
        k = 1 + floor(NNCI1a(n)/FL);
    
        i = 1 + floor(rem(NNCI1a(n),(FL*(k-1)))/max_x);
        if (rem(NNCI1a(n),(FL*(k-1)))/max_x) - floor(rem(NNCI1a(n),(FL*(k-1)))/max_x) == 0
            i = i - 1;
        end
            
        j = rem(NNCI1a(n),max_y*max_x)-(max_x*(i-1));
        if j == 0
            j = max_x;
        end
    
    else
    
        k = 1;
   
        i = 1 + floor(NNCI1a(n)/max_x);
        if (NNCI1a(n)/max_x) - floor((NNCI1a(n)/max_x)) == 0
            i = i - 1;
        end
   
        j = NNCI1a(n) - (max_x*(i-1));   
        if j == 0
            j = max_x;
        end
    
    end
    
    NNCI1a_dex(n,:) = [ i j k ];

end

for n = 1:length(NNCI2a)
      
    if NNCI2a(n) > FL
    
        k = 1 + floor(NNCI2a(n)/FL);
    
        i = 1 + floor(rem(NNCI2a(n),(FL*(k-1)))/max_x);
        if (rem(NNCI2a(n),(FL*(k-1)))/max_x) - floor(rem(NNCI2a(n),(FL*(k-1)))/max_x) == 0
            i = i - 1;
        end
    
        j = rem(NNCI2a(n),max_y*max_x)-(max_x*(i-1));
        if j == 0
            j = max_x;
        end
    
    else
    
        k = 1;
   
        i = 1 + floor(NNCI2a(n)/max_x);
        if (NNCI2a(n)/max_x) - floor((NNCI2a(n)/max_x)) == 0
            i = i - 1;
        end
   
        j = NNCI2a(n) - (max_x*(i-1)); 
        if j == 0
            j = max_x;
        end
    
    end
    
    NNCI2a_dex(n,:) = [ i j k ];

end

for n = 1:length(NNCI1b)
      
    if NNCI1b(n) > FL
    
        k = 1 + floor(NNCI1b(n)/FL);
    
        i = 1 + floor(rem(NNCI1b(n),(FL*(k-1)))/max_x);
        if (rem(NNCI1b(n),(FL*(k-1)))/max_x) - floor(rem(NNCI1b(n),(FL*(k-1)))/max_x) == 0
            i = i - 1;
        end
            
        j = rem(NNCI1b(n),max_y*max_x)-(max_x*(i-1));
        if j == 0
            j = max_x;
        end
    
    else
    
        k = 1;
   
        i = 1 + floor(NNCI1b(n)/max_x);
        if (NNCI1b(n)/max_x) - floor((NNCI1b(n)/max_x)) == 0
            i = i - 1;
        end
   
        j = NNCI1b(n) - (max_x*(i-1));   
        if j == 0
            j = max_x;
        end
    
    end
    
    NNCI1b_dex(n,:) = [ i j k ];

end

for n = 1:length(NNCI2b)
      
    if NNCI2b(n) > FL
    
        k = 1 + floor(NNCI2b(n)/FL);
    
        i = 1 + floor(rem(NNCI2b(n),(FL*(k-1)))/max_x);
        if (rem(NNCI2b(n),(FL*(k-1)))/max_x) - floor(rem(NNCI2b(n),(FL*(k-1)))/max_x) == 0
            i = i - 1;
        end
    
        j = rem(NNCI2b(n),max_y*max_x)-(max_x*(i-1));
        if j == 0
            j = max_x;
        end
    
    else
    
        k = 1;
   
        i = 1 + floor(NNCI2b(n)/max_x);
        if (NNCI2b(n)/max_x) - floor((NNCI2b(n)/max_x)) == 0
            i = i - 1;
        end
   
        j = NNCI2b(n) - (max_x*(i-1)); 
        if j == 0
            j = max_x;
        end
    
    end
    
    NNCI2b_dex(n,:) = [ i j k ];

end