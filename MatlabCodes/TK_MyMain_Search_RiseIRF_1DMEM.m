% TK_MyMain_Search_RiseIRF_1DMEM

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%LoadFileName_IRF = 'C:\Users\Steph\Dropbox (MIT)\Data_MIT\2D_FLCS_Data\Disc\ForPaper\irf_disc';
%LoadFileName_IRF = '~/Dropbox (MIT)/2D_FLCS/SimulationData/PhotonNumber/IRF';
LoadFileName_IRF = 'irf_liposome';

testImax =15;
RisePoint_IRF_start =285;
RisePoint_FL =300;

Tstart = 0 ;    % sec
Tend = 10^20 ;  % sec

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dT = [10^-4] ; ddT = [10^-5] ;  % not important, just set.

Equi0orNonequi1 = 0;
%TK_MyMain_Create2DFDC_cor_02
TK_MyMain_Create2DFDC_cor_SeparateData_BootStrap_v02

load(LoadFileName_IRF)

ticAAA = tic;
testI = 0;
while testI < testImax
    testI = testI + 1 ;
    
    RisePoint_IRF = (RisePoint_IRF_start-1) + testI 

    TK_MyMain_Fit_1DMEM_02

    if testI == 1
       RisePoint_Result_1DFDC_Mat_A = zeros(length(Result_1DFDC_Mat_A), testImax)  ;
       RisePoint_Result_1DFDC_Mat_M_Model_lin = zeros(length(Result_1DFDC_Mat_M_Model_lin), testImax)  ;
       RisePoint_Result_1DFDC_EstQ = zeros(testImax, 3) ;
    end
    RisePoint_Result_1DFDC_Mat_A(:, testI) = Result_1DFDC_Mat_A(:) ;
    RisePoint_Result_1DFDC_Mat_M_Model_lin(:, testI) = Result_1DFDC_Mat_M_Model_lin(:) ;
    
    RisePoint_Result_1DFDC_EstQ(testI,1) = RisePoint_FL ;
    RisePoint_Result_1DFDC_EstQ(testI,2) = RisePoint_IRF ;
    Var = size(Result_1DFDC_EstQ) ;
    RisePoint_Result_1DFDC_EstQ(testI,3) = Result_1DFDC_EstQ(Var(1),Var(2)) ;
    
    close all ;
end

toc(ticAAA)

prompt = 'save file as:';
filename = input(prompt,'s');
save(filename)

clear ticAAA


