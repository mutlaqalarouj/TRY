function [SP,P] = SPbrugge2(DX,DY,DZ,n)
%'TP9D09F3FD_MODEL_93.EGRID'
%'TP9D09F3FD_MODEL_93.UNRST'
% 
% function SP = SPeagebound(DX,DY,DZ,i)
% 
% EGRID ='TP9D09F3FD_MODEL_93.EGRID';
% F=readFfile(EGRID);
% 
% if i <= 9
%     FFILE =['TP9D09F3FD_MODEL_93.F000', int2str(i)];
%     G = readFfile(FFILE);
%     
% elseif (9<i) && (i<=99)
%     FFILE =['TP9D09F3FD_MODEL_93.F00', int2str(i)];
%     G = readFfile(FFILE);
%     
% elseif (99<i) && (i<=999)
%     FFILE =['TP9D09F3FD_MODEL_93.F0', int2str(i)];
%     G = readFfile(FFILE);
% end
%     
% EGRID = input('name of Eclipse.EGRID?');
% UNRST = input('name of Eclipse.UNRST?');
% n = input('time of run');
% DX = input('DX');
% DY = input('DY');
% DZ = input('DZ');
%DX = 120; DY = 1; DZ = 1;

% F=readfile(EGRID,0);
% G = readfile(UNRST,0,{'PRESSURE','SWAT','DIMENS'});

EGRID ='TP9D09F3FD_MODEL_93.FEGRID';
[F,NNC]=readFfilezz(EGRID);

if n <= 9
    FFILE =['TP9D09F3FD_MODEL_93.F000', int2str(n)];
    G = readFfilezzz(FFILE);
    
elseif (9<n) && (n<=99)
    FFILE =['TP9D09F3FD_MODEL_93.F00', int2str(n)];
    G = readFfilezzz(FFILE);
    
elseif (99<n) && (n<=999)
    FFILE =['TP9D09F3FD_MODEL_93.F0', int2str(n)];
    G = readFfilezzz(FFILE);
end
    
SWAT = G.SWAT(:);
PRES = G.PRESSURE(:);    
ATN = F.ACTNUM;
% ATN = ones(size(SWAT));

jp = 1;
    for r = 1:length(ATN)        
        if ATN(r) == 0            
            saturation(r) = 1;            
        else saturation(r) = SWAT(jp);            
            jp = jp+1;        
        end        
    end
    
X = reshape(saturation, DX, DY, DZ);


    for ip = 1:DZ
        
        SWAT3D(:,:,ip) = X(:,:,ip)';
    end  
clear j
jp=1;

    for r = 1:length(ATN)
        if ATN(r) == 0
        pressure(r) = 0;
        else pressure(r) = PRES(jp);
        jp = jp+1;
        end  
    end
    
X = reshape(pressure, DX, DY, DZ);

    for ip = 1:DZ
        PRES3D(:,:,ip) = X(:,:,ip)';
    end

% x,y,z coordinates
    COORD = F.COORD;
    COORD1= reshape(COORD,6,[]);
    XCORD = COORD1(1,1:DX+1);
    ylength = size(COORD1);
    YCORD = COORD1(2,1:DX+1:ylength(2));  
    ZCORN = F.ZCORN;
    z1 = reshape(ZCORN,2*DX,2*DY,2*DZ);
    for ip = 1:2*DZ
        z2(:,:,ip) = z1(:,:,ip)';
     end;
    z3 = z2(:,[1:2:(2*DX-1),2*DX],:);
    z4 = z3([1:2:(2*DY-1),2*DY],:,:);
    z5 = z4(:,:,[1:2:(2*DZ-1),2*DZ]);
    ZCORD = z5(39,39,1:DZ+1);
    
    PRES3D = double(PRES3D);
    SWAT3D = double(SWAT3D);
    XCORD = double(XCORD);
    YCORD = double(YCORD);
    ZCORD = double(ZCORD);
    
    
% Saturation S, Pressure P, x, y and z coordinates
    S = SWAT3D;
    P = PRES3D*100000*0.0689;
    x = XCORD*0.3048;
    y = YCORD*0.3048;
    z = ZCORD*0.3048;
    boundarycon;
    
   
    DX3D = DX;
    DY3D = DY+5;
    DZ3D = DZ+20;
    S3D = ones(DY3D,DX3D,DZ3D);
    P3D = zeros(DY3D,DX3D,DZ3D);
    boundary3D = zeros(DY3D,DX3D);
    S3D(6:DY3D,1:DX3D,11:DZ3D-10)=S;
    P3D(6:DY3D,1:DX3D,11:DZ3D-10)=P;
    boundary3D(6:DY3D,1:DX3D) = boundary;
    boundary3D(1:5,95:96) = 1;
    boundary = boundary3D;
    
    y3D = ones(1,DY3D+1);
    z3D = ones(1,DZ3D+1);
         
    y3D(1:5) = [(y(1)-600) (y(1)-480) (y(1)-360) (y(1)-240) (y(1)-120)];
    y3D(6:DY3D+1)=y;
     
    z3D(1:11) = [(z(1)-500) (z(1)-450) (z(1)-400) (z(1)-350) (z(1)-300) (z(1)-250) (z(1)-200) (z(1)-150) (z(1)-100) (z(1)-50) z(1)];           
    z3D(11:DZ3D-9) = z;
    z3D(21:30) = [(z(end)+50) (z(end)+100) (z(end)+150) (z(end)+200) (z(end)+250) (z(end)+300) (z(end)+350) (z(end)+400) (z(end)+450) (z(end)+500)]; 
    
%     DX = DX3D;
%     DY = DY3D;
%     DZ = DZ3D;
    P = P3D;
    S = S3D;
    y = y3D;
    z = z3D;
    ZNNC;
    NNCI;
    
% for j=1:(max_x)
%     for i=1:(max_y)
%         for k = 1:(max_z)
%             if S(i,j,k) ==0
%             S(i,j,k)=1;
%             end
%         end
%     end
% end

    clearvars -except x y z S P F G DX DY DZ DX3D DY3D DZ3D boundary NNC NNCI1a NNCI2a NNCI1b NNCI2b
%EK SOLVER
    Swirr = 0;
    Sor = 0;
    SALT = 0.5 * ones(DY3D,DX3D,DZ3D) ;
    C_ek = (-1.36*(SALT.^-0.9123)).*10^-9;
    %C_ek = -(SALT.^-1.213).*10^-9;
    
    for r = 1:length(SALT(:))        
        if SALT(r) < 0.09            
            tna(r) = 0.39;            
        elseif SALT(r) > 0.09
            tna(r) = 0.366-(0.0212*(log(SALT(r))));      
        end        
    end

    %for r = 1:length(S(:))        
        %if S(r) <= 0.2003            
            %C_ecc(r) = ((-0.0861.*TEMP(r))./SALT(r)).*10^-3;            
        %else %C_ecc(r) = ((-0.0861.*TEMP(r))./SALT(r)).*10^-3;
             % C_ecc(r) = ((0.0861.*TEMP(r).*(2*tna(r)-1))./SALT(r)).*10^-3;
             %C_ecc(r) = ((-0.0207.*TEMP(r))./SALT(r)).*10^-3;                
       % end        
    %end
    
%C_ec = reshape(C_ecc, DY, DX, DZ);


    %j = 1;
    %for r = 1:length(S(:))        
        %if S(r) <= 0.2003            
           % C_tec(r) = ((-0.1984*(log(SALT(r))))+0.5953)*10^-3;            
        %else C_tec(r) = ((-0.1984*(2*tna(r)-1)*(log(SALT(r))))+(1.059*tna(r)-0.5673))*10^-3;
            %C_tec(r) = ((-0.1984*(log(SALT(r))))+0.5953)*10^-3;
            %j = j+1;        
       % end        
    %end
    
%C_te = reshape(C_tec, DY, DX, DZ);

    m_ek = 0.6;
    %m_ec = 3;
    Swn = (S-Swirr)./(1.0-Swirr-Sor);
    Swn(Swn<0)=0;
    %sigma1 = (log(SALT./0.4783))/0.1046;
    sigma1 = (SALT-0.0388)/0.1252;
    sigma = sigma_Swa(sigma1, S);
    Cr_ek = Cek_relative(Swn, m_ek);
    %Cr_ec = Cec_relative(Swn, m_ec);
    L_ek = sigma.*C_ek.*Cr_ek;
    %L_ec = sigma.*C_ec.*Cr_ec;
    %L_te = sigma.*C_te.*Cr_ec;
       
    ATN1 = reshape(F.ACTNUM, DX, DY, DZ);

    for ip = 1:DZ
        
        ATN3D(:,:,ip) = ATN1(:,:,ip)';
    end 
    
    ATN4D(6:53,1:DX3D,11:19) = ATN3D;
    
    for k=1:DZ
        for j=1:DX
            for i=1:DY
                n=(k-1)*(DY*DX)+(j-1)*DY+i;
                
                if ATN4D(n) == 0
                        L_ek(n) = 0;
                else L_ek(n) = L_ek(n);
                end          
                
            end
        end
    end
    
%=========================================    
%2D and 3D SP Solver
%=========================================
%     L_ek(1:9720) = 0;
%     L_ek(1:2,:,21:100)=0;
%     L_ek(8:9,:,21:100)=0;
%     L_ek(:,1:2,21:100)=0;
%     L_ek(:,53:54,21:100)=0;
    
    lek_pressurerc;
    %lec_concentration;
    %lte_temperature;
    sigma_swatc;
     NNCO;
     NNCOV;
     NNCwa2;
     NNCa2;
     NNCswa2;  
     NNCsa2;
     NNCwb2;
     NNCb2;
     NNCswb2;  
     NNCsb2;
    B(isnan(B)) = 0;
    A(isnan(A)) = 0;
    %C(isnan(C)) = 0;
    %T(isnan(T)) = 0;
        
% solve equation
    [U,flag,relres,iter,resvec] = bicg(A,B*P(:),1.0e-10,10000);
    %[U,flag,relres,iter,resvec] = bicg(A,C*SALT(:),1.0e-10,10000,[]);
    %[U,flag,relres,iter,resvec] = bicgstab(A,C*SALT(:),1.0e-10,100000,[]);
    %[U,flag,relres,iter,resvec] = bicgstab(A,T*TEMP(:),1.0e-10,100000,[]);
    
SP = reshape(U,DY3D,DX3D,DZ3D);
%DV = SP - SP(5,14,1);


%Create output file for petrel image
%    for ip = 1:DZ
%        U2(:,:,ip) = U1(:,:,ip)';
%    end
%U3 = U2(:);

%pressure difference with U3(1) being the reference
%U4 = U3 - U3(1);


%fout = fopen('SP.GRDECL','w+');
       %fprintf(fout,'SP\n')
               % fprintf(fout, '%f\n', U4);
%fprintf(fout,'/')
%fclose(fout);
return