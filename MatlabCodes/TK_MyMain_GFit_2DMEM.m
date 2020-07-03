% TK_MyMain_GFit_2DMEM

%%
%%%%%%%%%%%%%% Fix1orNot0orAbs2orAmp3_Mat_A = [2,2,2,2] ;
%%%%%%%%%%%%%% Fix1orNot0orAbs2_Mat_G_dTfun = ones(NumOfState, NumOfState, length(dT)) .* 3 ;

% For fitting parameters
%%%%%%%%%%%%%% mi_TypeSelect = 0 ; % Select mi-Type used in MEM: 0 = average constant, 1 = Initial Tau distribution

DisplayFig = 1 ;    % Display figure during calculation (1) or not (0)

% For fitting parameters to calculate initial Mat_A and Mat_U
BreakFactor = 100 ;
RegulatorConst = 10^(-1+0) ;
RegulatorFactor = 1.4 ;
TrialNumFor_MinimizeEstQ = 1 ;  % Trial number for minimizing Q
TrialNumFor_RegulatorConst = 10 ;    % Trial number for increasing RegulatorConst

% For Global fitting parameters
G_BreakFactor = 100 ;
G_RegulatorConst = 10^(-1+0) ;
G_RegulatorFactor = 1.4;
G_TrialNumFor_MinimizeEstQ = 1 ;    % Global Trial number for minimizing Q
G_TrialNumFor_RegulatorConst = 100 ; % Global Trial number for increasing RegulatorConst

y0_dTfun = ones(length(dT), 1) ;
y0_dTfun = y0_dTfun .* Result_Short_y0 ;
Fix1orNot0orAbs2_y0_dTfun = ones(length(dT), 1) ;
Fix1orNot0orAbs2_y0_dTfun = Fix1orNot0orAbs2_y0_dTfun .* Fix1orNot0orAbs2_y0 ;

%% Same parameters as Fit_2DMEM_04
% Tau_min = 0.1 ;
% Tau_max = 10.1 ;
% Tau_step = 0.1 ;

% For fitting parameters
% UseCor1orNot0 = 0 ; % Use or not cor-matrix (after subtraction of longest-matrix to remove uncorrelated contribution)
% Linear0orLog1 = 1 ; % Fit linear data: 0, or log data: 1
% FitStartI = 30 ;    % Start point for fitting

% convoluted exp curve
% RisePoint_FL = 50 ;
% RisePoint_IRF = 50 ;
% RangePoint_IRF_min = 50 ;
% RangePoint_IRF_max = 1500 ;

%% Set parameters
Var_StartTime = tic ;

% % model 2D map
% if Tau_Select == 0
%     NumOfState = size(Tau_Initial) ;
%     NumOfState = NumOfState(1) ;
% elseif Tau_Select == 1
%     NumOfState = size(Tau_SplitPos) ;
%     NumOfState = NumOfState(1) ;
% elseif Tau_Select == 2
%     NumOfState = length(Tau_StateSelect) ;
% end
% 
% % Set Initial Tau distribution = Initial Mat_A
% Tau = [Tau_min : Tau_step : Tau_max]' ;  % array of tau  /ns
% NumOfComponent = length(Tau) ;
% Tau_Initial_distribution = zeros(NumOfComponent, NumOfState) ;
% 
% if Tau_Select == 0
%     K = 0 ;
%     while K < NumOfState
%         K = K + 1 ;
%         I = 0 ;
%         while I < NumOfComponent
%             I = I + 1 ;
%             VarAmp = Tau_Initial(K, 1) ;
%             VarTau0 = Tau_Initial(K, 2) ;
%             VarGaussWidth = Tau_Initial(K, 3) ;
%             Tau_Initial_distribution(I, K) = VarAmp * exp(-1*((Tau(I)-VarTau0)/VarGaussWidth)^2) ;
%         end
%     end
% elseif Tau_Select == 1        
%     K = 0 ;
%     while K < NumOfState
%         K = K + 1 ;
%         Vmin = Tau_SplitPos(K, 1) ;
%         Vmax = Tau_SplitPos(K, 2) ;
%         if (Vmin == 0) && (Vmax == 0)
%             Tau_Initial_distribution(:, K) = 1 ;
%         else
%             I = 0 ;
%             while I < NumOfComponent
%                 I = I + 1 ;
%                 if (Vmin <= Tau(I)) && (Tau(I) < Vmax)
%                     Tau_Initial_distribution(I, K) = sum(Tau_Initial(I,:)) ;
%                 elseif Tau(I) >= Vmax
%                     break
%                 end
%             end
%         end
%     end
% elseif Tau_Select == 2     
%     K = 0 ;
%     while K < NumOfState
%         K = K + 1 ;
%         Var = Tau_StateSelect(K) ;
%         Tau_Initial_distribution(:,K) = Tau_Initial(:,Var) ;
%     end
% end
% 
% clear VarAmp VarTau0 VarGaussWidth Vmin Vmax
% 
% Set ExpCurve
% [ExpCurve, ExpCurve_Binned_lin, ExpCurve_Binned_log] = TK_CreateExpCurve(...
%     Tau, xdata, IRF, NumOfComponent,...
%     tMin, tMax, tStep, lint_BinFactor, Mat_2DFDC_log_t,...
%     RisePoint_FL, RisePoint_IRF, RangePoint_IRF_min, RangePoint_IRF_max) ;


%% Get Mat_U by using Mat_A obtained above (To estimate Initial_Mat_U)
display('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@')
display('@@@@ Calculating Mat_U for each dT with Mat_A estimated at dT=shortest @@@@')
display('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@')

% fit Mat_U with fixed Mat_A
Var_Fix1orNot0orAbs2orAmp3_Mat_A = ones(NumOfState, 1) ; % Mat_A is fixed
Tau_Initial_distribution = Result_Short_Mat_A ;
Var_y0 = y0_dTfun ./ (Mat_2DFDC_lin_t(2) - Mat_2DFDC_lin_t(1)) ;

Tmax = length(dT) ;
Initial_Mat_G_dTfun = zeros(NumOfState, NumOfState, Tmax) ;
Initial_y0_dTfun = zeros(Tmax) ;

T = 0 ;
while T < Tmax
    T = T + 1 ;
    
    Fix1orNot0orAbs2_Mat_G = Fix1orNot0orAbs2orSym3_Mat_G_dTfun(:,:,T) ;
    
    InMat_2DFDC_lin = Mat_2DFDC_lin(:, :, T) ;
    InMat_2DFDC_log = Mat_2DFDC_log(:, :, T) ;
    if UseCor1orNot0 == 1
        InMat_2DFDC_cor_lin = Mat_2DFDC_cor_lin(:, :, T) ;
        InMat_2DFDC_cor_log = Mat_2DFDC_cor_log(:, :, T) ;
    elseif UseCor1orNot0 == 0
        InMat_2DFDC_cor_lin = Mat_2DFDC_lin(:, :, T) ;
        InMat_2DFDC_cor_log = Mat_2DFDC_log(:, :, T) ;    
    end
    
    [Result_estimates, Vart, Result_Mat_M_2DFLC,...
        Result_Mat_M_Model_lin, Result_Mat_M_Model_log,...
        Result_Mat_A, Result_Mat_G, Result_y0, Result_EstQ] =...
        TK_FitF_MinimizeQ_09(...
        Mat_2DFDC_lin_t, InMat_2DFDC_lin, InMat_2DFDC_cor_lin,...
        Mat_2DFDC_log_t, InMat_2DFDC_log, InMat_2DFDC_cor_log,...
        Tau, Tau_Initial_distribution, Var_Fix1orNot0orAbs2orAmp3_Mat_A, Fix1orNot0orAbs2_Mat_G, Var_y0(T), Fix1orNot0orAbs2_y0_dTfun(T),...
        ExpCurve_Binned_lin, ExpCurve_Binned_log,...
        NumOfState, RegulatorConst, RegulatorFactor,...
        Linear0orLog1, FitStartI, TrialNumFor_MinimizeEstQ, TrialNumFor_RegulatorConst, mi_TypeSelect, BreakFactor,...
        tMin, tMax, 0) ;
    
    Initial_Mat_G_dTfun(:, :, T) = Result_Mat_G(:, :) ;
    Initial_y0_dTfun(T) = Result_y0 ;
end

% Tmax = length(dT) ;
% Initial_Mat_G_dTfun = zeros(NumOfState, NumOfState, Tmax) ;
% T = 0 ;
% while T < Tmax
%     T = T + 1 ;
%     Initial_Mat_G_dTfun(:, :, T) = eye(NumOfState, NumOfState) ;
% end

%% Global fitting to get Result_Mat_A and Result_Mat_U_dTfun
display('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@')
display('@@@@ Calculating Mat_A (global variation) and Mat_U (local variation) for all dT data sets @@@@')
display('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@')

% Fix1orNot0orAbs2orSym3_Mat_G_dTfun = zeros(length(dT), 1) ;
% Fix1orNot0orAbs2orSym3_Mat_G_dTfun(:) = 2 ;

Initial_Mat_A_Global = Result_Short_Mat_A ;
% Initial_Mat_G_dTfun = Initial_Mat_G ;
% Initial_Mat_A_Global = Gresult_Mat_A_Global ;
% Initial_Mat_G_dTfun = Gresult_Mat_G_dTfun ;
Var_y0 = y0 / (Mat_2DFDC_lin_t(2) - Mat_2DFDC_lin_t(1)) ;

InMat_2DFDC_lin     = Mat_2DFDC_lin(:, :, :) ;
InMat_2DFDC_log     = Mat_2DFDC_log(:, :, :) ;
if UseCor1orNot0 == 1
    InMat_2DFDC_cor_lin = Mat_2DFDC_cor_lin ;
    InMat_2DFDC_cor_log = Mat_2DFDC_cor_log ;
elseif UseCor1orNot0 == 0
    InMat_2DFDC_cor_lin = Mat_2DFDC_lin ;
    InMat_2DFDC_cor_log = Mat_2DFDC_log ; 
end

[Gresult_estimates_dTfun, Gresult_Mat_M_2DFLC_dTfun, Gresult_Mat_2DFDC_lin_t, Gresult_Mat_2DFDC_log_t,...
    Gresult_Mat_M_Model_lin_dTfun, Gresult_Mat_M_Model_log_dTfun,...
    Gresult_Mat_A_Global, Gresult_Mat_G_dTfun, GResult_y0_dTfun, Gresult_EstQ_Global] =...
    TK_GFitF_MinimizeQ_04(...
    Mat_2DFDC_lin_t, InMat_2DFDC_lin, InMat_2DFDC_cor_lin,...
    Mat_2DFDC_log_t, InMat_2DFDC_log, InMat_2DFDC_cor_log,...
    Tau, Initial_Mat_A_Global, Initial_Mat_G_dTfun, Fix1orNot0orAbs2orAmp3_Mat_A, Fix1orNot0orAbs2orSym3_Mat_G_dTfun, Initial_y0_dTfun, Fix1orNot0orAbs2_y0_dTfun, ...
    ExpCurve_Binned_lin, ExpCurve_Binned_log,...
    NumOfState, G_RegulatorConst, G_RegulatorFactor,...
    Linear0orLog1, FitStartI, G_TrialNumFor_MinimizeEstQ, G_TrialNumFor_RegulatorConst, mi_TypeSelect, G_BreakFactor,...
    tMin, tMax, DisplayFig) ;

GResult_y0_dTfun = GResult_y0_dTfun .* (Mat_2DFDC_lin_t(2) - Mat_2DFDC_lin_t(1)) ;

%%

display(' ')
display(strcat('Total elasped time is...___', num2str(toc(Var_StartTime)), '   seconds'))
display(' ')

clear I Initial_Mat_A_Global InMat_2DFDC_cor_dTfun InMat_2DFDC_cor_lin InMat_2DFDC_cor_log
clear InMat_2DFDC_dTfun InMat_2DFDC_lin InMat_2DFDC_log InMat_2DFDC_t K T Var1 Var2 Var3 Var_StartTime Vart
