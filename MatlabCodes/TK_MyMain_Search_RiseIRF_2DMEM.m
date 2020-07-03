% TK_MyMain_Search_RiseIRF_2DMEM

%%
LoadFileName_IRF = 'C:\Users\Steph\Dropbox (MIT)\Data_MIT\2D_FLCS_data\Liposome\ForPaper\irf_liposome';
%LoadFileName_IRF = 'C:\Users\Julianne\Dropbox (Personal)\MIT\data\LHCSR3\FLC\FLC sample datga\SI_Fig3a_increasePhotonNum\IRF_norm' ;
%LoadFileName_IRF = 'irf_liposome' ;
%LHCSR
Hozzon_estimates = [0,...
                    1, 0.4, 0.3,...
                    1, 2.2, 0.3] ;
% LHCII
% Hozzon_estimates = [0,...
%                     1, 0.6, 0.2,...
%                     1, 1.5, 0.2,...
%                     1, 2.5, 0.2] ;
% 
%dT = [10^-4, 10^-3, 10^-2, 10^-1, 10^0, 10^1, 10^2] ;    % sec
dT = [1e-4] ;    % sec
ddT = 10^-5 ;

NumOfState = floor(length(Hozzon_estimates)/3) ;
Tau_Select = 0 ;

Tstart = 0 ;    % sec
Tend = 10^20 ;  % sec

testImax = 20;
RisePoint_IRF_start = 310 - floor(testImax/2) ;

RisePoint_FL = 300 ;
%%
load(LoadFileName_IRF)

TK_MyMain_Create2DFDC_cor_SeparateData_BootStrap_v02

ticAAA = tic ;
testI = 0 ;
while testI < testImax
    testI = testI + 1 ;
    
    RisePoint_IRF = (RisePoint_IRF_start-1) + testI
    
    Fix1orNot0orAbs2orAmp3_Mat_A = ones(NumOfState,1) * 2 ;%ones(26, 1) .* 3 ;
    Fix1orNot0orAbs2orSym3_Mat_G = ones(NumOfState, NumOfState) .* 3 ;%- (2 .* eye(NumOfState, NumOfState)) ;
    TK_MyMain_Fit_2DMEM_04
    
    Fix1orNot0orAbs2orAmp3_Mat_A  = ones(NumOfState,1) * 2 ;%ones(26, 1) .* 3;%
    Fix1orNot0orAbs2orSym3_Mat_G_dTfun = ones(NumOfState, NumOfState, length(dT)) .* 3 ;
    TK_MyMain_GFit_2DMEM
    
    
    if testI == 1
        Var = size(Gresult_Mat_A_Global) ;
        RisePoint_Result_2DFDC_Mat_A = zeros(Var(1), Var(2), testImax)  ;
        Var = size(Gresult_Mat_M_Model_lin_dTfun) ;
        RisePoint_Result_2DFDC_Mat_M_Model_lin_dT = zeros(Var(1), Var(2), testImax)  ;
        RisePoint_Result_2DFDC_EstQ = zeros(testImax, 2+3) ;
    end
    RisePoint_Result_2DFDC_Mat_A(:, :, testI) = Gresult_Mat_A_Global(:,:) ;
    % keep Model only at dT(1)
    RisePoint_Result_2DFDC_Mat_M_Model_lin_dT(:, :, testI) = Gresult_Mat_M_Model_lin_dTfun(:, :, 1) ;
    
    RisePoint_Result_2DFDC_EstQ(testI,1) = RisePoint_FL ;
    RisePoint_Result_2DFDC_EstQ(testI,2) = RisePoint_IRF ;
    Var = size(Gresult_EstQ_Global) ;
    RisePoint_Result_2DFDC_EstQ(testI,3:5) = Gresult_EstQ_Global(Var(1),2:Var(2)) ;
    
    close all ;
end

toc(ticAAA)

clear ticAAA Var
clear tt1 kin1

save('Search_RiseIRF_2DMEM_liposome')

