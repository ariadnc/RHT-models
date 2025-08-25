%**********************************************************************************
%
%               A MATLAB FUNCTION WRITTEN BY HADI BORDBAR
%
%                     DEPARTMENT OF CIVIL ENGINEERING
%                       AALTO UNIVERSITY, FINLAND
%
% Discription: A Matlab function (suplimentary material) to calculate the 
% model coefficients and emissivity by 
% the pressure dependent WSGG model published in:
% Bordbar H., Coelho F.,Fraga G., Franca, F., Hostikka, S. 
% Pressure-dependent weighted-sum-of-gray-gases models for heterogeneous -
% mixtures at sub- and super-atmospheric pressure. 
% International Journal of Heat and Mass Transfer xxx (2021) 121207
%=========================================================================================
function [a_wsgg, kk_wsgg, emissivity]=Bordbar_PWSGG_MatlabFunction_IJHMT2021(Pt,T,x_h2o,x_co2, L);

% INPUTS:    
%           Total pressure (Pt) [atm] [0.1-80]
%           Temperature T [K] 300-3000
%           Mole fraction of H2O (xh2o) [0...1]
%           Mole fraction of CO2 (xco2) [0...1]
%           Path length (L) [m]; optional if emissivity of a path is needed
%
% OUTPUTS    
%           Gray gas absorption coefficients (kk_i) [atm ^-1 m^-1]
%           Weighting factors a_i
%           Emissivity of the mixture [0....1] if needed (path length is
%           given in the inputs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loading the model parameters corresponded to current Pt */
format long;
NPt=10; Ng=4; NPoly1=4; NPoly2=4; 

Pt_ref(:)=[0.1 0.5 1.0  2.0 5.0 10.0 20.0 30.0 50.0 80.0]; % total pressure of the mixture
 
if(x_co2 ~= 0)
    Mr_case=x_h2o/x_co2;
end

for i=1:NPt
    if(Pt==Pt_ref(i))
        iPt=i;
    end
end

if(iPt<10)
     Model_Address=['Final_PWSGG_Model_Coef_Pt00',num2str(iPt),'.dat'];
else
     Model_Address=['Final_PWSGG_Model_Coef_Pt0',num2str(iPt),'.dat'];
end
% model parameter
coef=dlmread(Model_Address,'',[27 0 (27+Ng*(NPoly1+2)-1) 3+NPoly2]);

count=1;
for iGas=1:Ng
    for iPoly1 =1:NPoly1+2
        poly_F=coef(count,(4:4+NPoly2));
        Coef_Data(iGas, iPoly1)=polyval(poly_F, Mr_case);
        count=count+1;
    end
end

A= zeros(Ng,NPoly1);
KK=zeros(Ng);
for i=1:Ng
    for j=1:NPoly1+1
        A(i,j)=Coef_Data(i,j);
    end
    KK(i)=Coef_Data(i,NPoly1+2);
end

Tr=T./1500;

sum_a=0;
for i=1:Ng
    a_wsgg(i)= 0;
    for j=1:NPoly1+1
        a_wsgg(i)=a_wsgg(i)+A(i,j)*Tr^(j-1);
    end
    sum_a=sum_a+a_wsgg(i);    
    kk_wsgg(i)=KK(i)*(x_co2+x_h2o)*Pt;
end

a_wsgg(Ng+1)=1-sum_a;
kk_wsgg(Ng+1)=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    Calculation of total emissivity(eq. 7)                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

emissivity=0;
if (nargin == 5) % only if emissivity of a path is needed
    PtL=(x_h2o+x_co2)*Pt*L;
    for i=1:4
        emissivity=emissivity+a_wsgg(i)*(1-exp(-kk_wsgg(i)*PtL));
    end
end
end
        