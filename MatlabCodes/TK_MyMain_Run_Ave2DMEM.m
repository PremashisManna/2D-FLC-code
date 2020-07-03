% TK_MyMain_Run_Ave2DMEM

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LoadFileName_IRF = 'IRF.mat' ; % load IRF
%LoadFileName_IRF = 'C:\Users\Julianne\Dropbox (Personal)\MIT\data\LHCSR3\FLC\FLC sample datga\3 states\IRF_norm' ;

Hozzon_estimates = [0,...
                    1, 1, 0.3,...
                    1, 3, 0.3] ;
% Hozzon_estimates = [0,...
%                     1, 0.6, 0.2,...
%                     1, 1.5, 0.2,...
%                     1, 2.5, 0.2] ;

                
Var_RisePoint_Imax =  5 ;
Var_RisePoint_Icenter = 300;

Tstart = 0 ;                       % sec
Tend = 10^20 ;                     % sec

%SaveFileName = 'Ave2D_BootStrap1_LHCSR1_V_7_5th' ;
SaveFileName = 'simulated_Ave2D' ;

Tau_Select = 0 ; 

if Tau_Select ~= 2  % If Tau_Select is not 2, Mat_A is free
    NumOfState = floor(length(Hozzon_estimates)/3) ;
    Fix1orNot0orAbs2orAmp3_Mat_A = ones(NumOfState,1) * 1; %changed from 2 to 1 on 190117
else
    LoadFileName_Stored_Mat_A = 'Gresult_Mat_A_Global' ;     % If Tau_Select = 2, Mat_A is fixed
    Var = size(LoadFileName_Stored_Mat_A) ;
    NumOfState = Var(2) ;
    Fix1orNot0orAbs2orAmp3_Mat_A = ones(NumOfState,1) * 1 ;
end

dT = [1e-4 3e-4 1e-3 3e-3 1e-2 3e-2 1e-1 3e-1 1 3];    % sec
ddT = 1e-5;

Equi0orNonequi1 = 0;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Var_RisePoint_I = 0 ;
while Var_RisePoint_I < Var_RisePoint_Imax
    Var_RisePoint_I = Var_RisePoint_I + 1 ;
    RisePoint_IRF = (Var_RisePoint_Icenter-1) - floor(Var_RisePoint_Imax/2) + Var_RisePoint_I ;
    
    if Tau_Select == 2
        load(LoadFileName_Stored_Mat_A) ;
    end
        
    load(LoadFileName_IRF)
    
    %%%%%%%%%%%%% fitting %%%%%%%%%%%%%%%%
    %TK_MyMain_Run
    if Var_RisePoint_I == 1
        TK_MyMain_Create2DFDC_cor_SeparateData_BootStrap_v02
    end
    
    if Equi0orNonequi1 == 0
        Fix1orNot0orAbs2orSym3_Mat_G = ones(NumOfState, NumOfState) .* 3 ;%- (2 .* eye(NumOfState, NumOfState)) ;
    elseif Equi0orNonequi1 == 1
        Fix1orNot0orAbs2orSym3_Mat_G = ones(NumOfState, NumOfState) .* 2 ;%- (2 .* eye(NumOfState, NumOfState)) ;
    end
    TK_MyMain_Fit_2DMEM_04
    
    if Equi0orNonequi1 == 0
        Fix1orNot0orAbs2orSym3_Mat_G_dTfun = ones(NumOfState, NumOfState, length(dT)) .* 3 ;
    elseif Equi0orNonequi1 == 1
        Fix1orNot0orAbs2orSym3_Mat_G_dTfun = ones(NumOfState, NumOfState, length(dT)) .* 2 ;
    end
    TK_MyMain_GFit_2DMEM
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if Var_RisePoint_I == 1
        clear kin1 tt1
        FileName = strcat(SaveFileName,'_IRFave') ;
        save(FileName)
    end
    
    clearvars -except...
        Gresult_estimates_dTfun Gresult_EstQ_Global Gresult_Mat_2DFDC_lin_t ...
        Gresult_Mat_2DFDC_log_t Gresult_Mat_A_Global Gresult_Mat_G_dTfun Gresult_Mat_M_2DFLC_dTfun ...
        Gresult_Mat_M_Model_lin_dTfun Gresult_Mat_M_Model_log_dTfun dT ddT Tau mi_TypeSelect UseCor1orNot0 ...
        Hozzon_estimates RisePoint_IRF Var_RisePoint_I Var_RisePoint_Imax Var_RisePoint_Icenter LoadFileName SaveFileName...
        Fix1orNot0orAbs2orAmp3_Mat_A Tau_Select LoadFileName_Stored_Mat_A NumOfState MeasurementTime_s...
        Tstart DeltaT Tend DataNameList LoadFileName_IRF;
    
    FileName = strcat(SaveFileName,'_IRF',num2str(RisePoint_IRF)) ;
    save(FileName)
    
    if Var_RisePoint_I == 1
        FileName = strcat(SaveFileName,'_IRFave.mat') ;
        load(FileName)
    elseif Var_RisePoint_I ~= 1
        Var_RisePoint_A1 = Gresult_estimates_dTfun ;
        Var_RisePoint_A2 = Gresult_EstQ_Global ;
        Var_RisePoint_A3 = Gresult_Mat_2DFDC_lin_t ;
        Var_RisePoint_A4 = Gresult_Mat_2DFDC_log_t ;
        Var_RisePoint_A5 = Gresult_Mat_A_Global ;
        Var_RisePoint_A6 = Gresult_Mat_G_dTfun ;
        Var_RisePoint_A7 = Gresult_Mat_M_2DFLC_dTfun ;
        Var_RisePoint_A8 = Gresult_Mat_M_Model_lin_dTfun ;
        Var_RisePoint_A9 = Gresult_Mat_M_Model_log_dTfun ;
        Var_RisePoint_A10 = Var_RisePoint_I ;
        
        clearvars -except ...
            Var_RisePoint_A1 Var_RisePoint_A2 Var_RisePoint_A3 Var_RisePoint_A4 Var_RisePoint_A5 Var_RisePoint_A6 ...
            Var_RisePoint_A7 Var_RisePoint_A8 Var_RisePoint_A9 Var_RisePoint_A10 ...
            SaveFileName Fix1orNot0orAbs2orAmp3_Mat_A Tau_Select LoadFileName_Stored_Mat_A NumOfState LoadFileName_IRF;
        
        FileName = strcat(SaveFileName,'_IRFave.mat') ;
        load(FileName)
        
        Gresult_estimates_dTfun = Var_RisePoint_A1 + Gresult_estimates_dTfun ;
        Gresult_EstQ_Global = Var_RisePoint_A2 + Gresult_EstQ_Global ;
        Gresult_Mat_2DFDC_lin_t = Var_RisePoint_A3 + Gresult_Mat_2DFDC_lin_t ;
        Gresult_Mat_2DFDC_log_t = Var_RisePoint_A4 + Gresult_Mat_2DFDC_log_t ;
        Gresult_Mat_A_Global = Var_RisePoint_A5 + Gresult_Mat_A_Global ;
        Gresult_Mat_G_dTfun = Var_RisePoint_A6 + Gresult_Mat_G_dTfun ;
        Gresult_Mat_M_2DFLC_dTfun = Var_RisePoint_A7 + Gresult_Mat_M_2DFLC_dTfun ;
        Gresult_Mat_M_Model_lin_dTfun = Var_RisePoint_A8 + Gresult_Mat_M_Model_lin_dTfun ;
        Gresult_Mat_M_Model_log_dTfun = Var_RisePoint_A9 + Gresult_Mat_M_Model_log_dTfun ;
        Var_RisePoint_I = Var_RisePoint_A10 ;   
        
        if Var_RisePoint_I == Var_RisePoint_Imax
            Gresult_estimates_dTfun = Gresult_estimates_dTfun / Var_RisePoint_Imax ;
            Gresult_EstQ_Global = Gresult_EstQ_Global / Var_RisePoint_Imax ;
            Gresult_Mat_2DFDC_lin_t = Gresult_Mat_2DFDC_lin_t / Var_RisePoint_Imax ;
            Gresult_Mat_2DFDC_log_t = Gresult_Mat_2DFDC_log_t / Var_RisePoint_Imax ;
            Gresult_Mat_A_Global = Gresult_Mat_A_Global / Var_RisePoint_Imax ;
            Gresult_Mat_G_dTfun = Gresult_Mat_G_dTfun / Var_RisePoint_Imax ;
            Gresult_Mat_M_2DFLC_dTfun = Gresult_Mat_M_2DFLC_dTfun / Var_RisePoint_Imax ;
            Gresult_Mat_M_Model_lin_dTfun = Gresult_Mat_M_Model_lin_dTfun / Var_RisePoint_Imax ;
            Gresult_Mat_M_Model_log_dTfun = Gresult_Mat_M_Model_log_dTfun / Var_RisePoint_Imax ;
        end
        
        clear Var_RisePoint_A1 Var_RisePoint_A2 Var_RisePoint_A3 Var_RisePoint_A4 Var_RisePoint_A5 Var_RisePoint_A6 ...
            Var_RisePoint_A7 Var_RisePoint_A8 Var_RisePoint_A9 Var_RisePoint_A10 ;
        
        FileName = strcat(SaveFileName,'_IRFave') ;
        save(FileName)
    end

    close all
end

clearvars -except ...
    Hozzon_estimates RisePoint_IRF Var_RisePoint_I Var_RisePoint_Imax Var_RisePoint_Icenter LoadFileName SaveFileName ...
    Fix1orNot0orAbs2orAmp3_Mat_A Tau_Select LoadFileName_Stored_Mat_A NumOfState...
    Tstart DeltaT Tend DataNameList LoadFileName_IRF;
