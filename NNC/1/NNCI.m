max_z = 28;
max_x = 139;
max_y = 53;
wa=1;
wb=1;


for k=1:max_z
        for j=1:max_x-1
            for i=1:max_y-1
                n=(k-1)*(max_y*max_x)+(i-1)*max_x+j;
                
                if Z(i*2,j*2,k) ~= Z(i*2+1,j*2,k)
                    v=0;
                    for ww=1:max_z
                        
                        if Z(i*2,j*2,k) < Z(i*2+1,j*2,ww) && Z(i*2+1,j*2,ww) < Z(i*2,j*2,k+1)
                            
                            NNCI1a(wa) = n;
                            NNCI2a(wa) = (ww-1)*(max_y*max_x)+(i)*max_x+j;
                            if k > 1 && v == 0
                                wa = wa+1;
                                NNCI1a(wa) = n;
                                NNCI2a(wa) = (ww-1)*(max_y*max_x)+(i)*max_x+j;
                            end
                            wa = wa+1;
                            v = 1;
                            
                        end
                        
                    end
                end
                    
                if Z(i*2,j*2,k) ~= Z(i*2,j*2+1,k)
                    v=0;
                    
                    for ww=1:max_z
                        if Z(i*2,j*2,k) < Z(i*2,j*2+1,ww) && Z(i*2,j*2+1,ww) < Z(i*2,j*2,k+1)
                            
                            NNCI1b(wb) = n;
                            NNCI2b(wb) = ((ww-1)*(max_y*max_x)+(i-1)*max_x+j)+1;
                            if k > 1 && v == 0
                                wb = wb+1;
                                NNCI1b(wb) = n;
                                NNCI2b(wb) = ((ww-1)*(max_y*max_x)+(i)*max_x+j)+1;
                            end
                            wb = wb+1;
                            v = 1;
                        end
                        
                    end 
                    
                end
            end
        end
end


