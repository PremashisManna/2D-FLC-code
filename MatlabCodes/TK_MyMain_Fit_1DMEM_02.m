% TK_MyMain_Fit_1DMEM_02

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tau_min = 0.05 ;
Tau_max = 7.05 ;
Tau_step = 0.05 ;

Tau_Select = 0 ;
if Tau_Select == 0
    Tau_Initial = [10000, 3, 10^10*1 ] ; % raw = NumOfState,   column = [Amp, Tau0, GaussWidth]
elseif Tau_Select == 1
    Tau_SplitPos = [] ;
    Tau_Initial = Stored_Result_Short_Mat_A ;
end

% fit Mat_A
Fix1orNot0orAbs2 = [2] ; % column 1 for Mat_A

% initial y0   y0 should be baseline intensity of 1D-FD at linear scale 
y0 = 5000;
% fit y0
Fix1orNot0orAbs2_y0 = [2];

% For fitting parameters
UseCor1orNot0 = 0 ; % Use or not cor-matrix (after subtraction of longest-matrix to remove uncorrelated contribution)
mi_TypeSelect = 0 ; % Select mi-Type used in MEM: 0 = average constant, 1 = Initial Tau distribution
Linear0orLog1 = 1 ; % Fit linear data: 0, or log data: 1
FitStartI = 1 ;    % Start point for fitting

DisplayFig = 1 ;    % Display figure during calculation (1) or not (0)

% For fitting parameters to calculate initial Mat_A and Mat_U
BreakFactor = 100 ;
RegulatorConst = 10^-1 ;
RegulatorFactor = 1.40 ;
TrialNumFor_MinimizeEstQ = 1 ;  % Trial number for minimizing Q
%for Ntrial=1:trialLength
TrialNumFor_RegulatorConst = 200;    % Trial number for increasing RegulatorConst

Tau = [Tau_min : Tau_step : Tau_max]' ;  % array of tau  /ns

trialLength=100;
displayMatA=zeros(length(Tau),trialLength);
displayResultsDecay = zeros(732,trialLength);


% convoluted exp curve
%%%%%%%%%%%%%%%%%%%%%%%%%%% RisePoint_FL = 300 ;
%%%%%%%%%%%%%%%%%%%%%%%%%%% RisePoint_IRF = 300-20 ;
RangePoint_IRF_min = 50 ;
RangePoint_IRF_max = length(xdata)-70 ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set parameters
Var_StartTime = tic ;

% model 2D map
if Tau_Select == 0
    NumOfState = size(Tau_Initial) ;
    NumOfState = NumOfState(1,1) ;
elseif Tau_Select == 1
    NumOfState = length(Tau_SplitPos) + 1 ;
end

% Set Initial Tau distribution = Initial Mat_A

NumOfComponent = length(Tau) ;
Tau_Initial_distribution = zeros(NumOfComponent, 1) ;

if Tau_Select == 0
    K = 0 ;
    while K < NumOfState
        K = K + 1 ;
        I = 0 ;
        while I < NumOfComponent
            I = I + 1 ;
            VarAmp = Tau_Initial(K, 1) ;
            VarTau0 = Tau_Initial(K, 2) ;
            VarGaussWidth = Tau_Initial(K, 3) ;
            Tau_Initial_distribution(I) = Tau_Initial_distribution(I) + VarAmp * exp(-1*((Tau(I)-VarTau0)/VarGaussWidth)^2) ;
        end
    end
end

if Tau_Select == 1        
    K = 1 ;
    Vmin = min(Tau) ;
    if NumOfState == 1
        Vmax = max(Tau)+abs(max(Tau)) ; 
    else
        Vmax = Tau_SplitPos(K) ;
    end
    I = 0 ;
    while I < NumOfComponent
        I = I + 1 ;
        if (Vmin <= Tau(I)) && (Tau(I) < Vmax)
            Tau_Initial_distribution(I) = Tau_Initial_distribution(I) + Tau_Initial(I) ;
        else
            I = I - 1 ;
            K = K + 1 ;
            Vmin = Vmax ;
            if K <= length(Tau_SplitPos)
                Vmax = Tau_SplitPos(K) ;
            else
                Vmax = max(Tau) * 2 ;
            end
        end
    end
end

clear VarAmp VarTau0 VarGaussWidth Vmin Vmax

% Set ExpCurve
[ExpCurve, ExpCurve_Binned_lin, ExpCurve_Binned_log] = TK_CreateExpCurve(...
    Tau, xdata, IRF, NumOfComponent,...
    tMin, tMax, tStep, lint_BinFactor, Mat_2DFDC_log_t,...
    RisePoint_FL, RisePoint_IRF, RangePoint_IRF_min, RangePoint_IRF_max) ;

% Adjust amplitude of Mat_A
Var1 = sum(sum(Tau_Initial_distribution))  ;
if Linear0orLog1 == 0
    Var2 = max(max(ExpCurve_Binned_lin)) ;
    if UseCor1orNot0 == 0
        Var3 = max(max(max(Mat_1DFDC_lin))) ;
    elseif UseCor1orNot0 == 1
        Var4 = max(max(max(Mat_1DFDC_cor_lin))) ;
    end
elseif Linear0orLog1 == 1
    Var2 = max(max(ExpCurve_Binned_log)) ;
    if UseCor1orNot0 == 0
        Var3 = max(max(max(Mat_1DFDC_log))) ;
    elseif UseCor1orNot0 == 1
        Var4 = max(max(max(Mat_1DFDC_cor_log))) ;
    end
end
Var = (Var1^1) * (Var2^1) ;
if UseCor1orNot0 == 0 
    Var = sqrt(Var3 / Var) ;
elseif UseCor1orNot0 == 1
    Var = sqrt(Var4 / Var) ;
end
Tau_Initial_distribution = Tau_Initial_distribution .* Var ;

Var_y0 = y0 / (Mat_2DFDC_lin_t(2) - Mat_2DFDC_lin_t(1)) ;

%%
%% Get Mat_A by fitting dT = shortest data (To estimate Initial_Mat_A)
display('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@')
display('@@@@ Calculating Mat_A for dT=shortest @@@@')
display('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@')

% Get Result_Short_Mat_A
InMat_1DFDC_lin = Mat_1DFDC_lin ;
InMat_1DFDC_log = Mat_1DFDC_log ;
if UseCor1orNot0 == 1
    InMat_1DFDC_cor_lin = Mat_1DFDC_cor_lin ;
    InMat_1DFDC_cor_log = Mat_1DFDC_cor_log ;
elseif UseCor1orNot0 == 0
    InMat_1DFDC_cor_lin = Mat_1DFDC_lin ;
    InMat_1DFDC_cor_log = Mat_1DFDC_log ;
end

[Result_1DFDC_estimates, ~, Result_1DFDC_Mat_M_Model_lin, Result_1DFDC_Mat_M_Model_log,...
    Result_1DFDC_Mat_A, Result_y0, Result_1DFDC_EstQ] =...
    TK_FitF_1DMEM_MinimizeQ_02(...
    Mat_2DFDC_lin_t, InMat_1DFDC_lin, InMat_1DFDC_cor_lin,...
    Mat_2DFDC_log_t, InMat_1DFDC_log, InMat_1DFDC_cor_log,...
    Tau, Tau_Initial_distribution, Fix1orNot0orAbs2, Var_y0, Fix1orNot0orAbs2_y0,...
    ExpCurve_Binned_lin, ExpCurve_Binned_log,...
    RegulatorConst, RegulatorFactor,...
    Linear0orLog1, FitStartI, TrialNumFor_MinimizeEstQ, TrialNumFor_RegulatorConst, mi_TypeSelect, BreakFactor,...
    tMin, tMax, DisplayFig) ;

% if Fix1orNot0orAbs2_y0 ~= 1
%     Var = length(Result_1DFDC_estimates) ;
%     Result_1DFDC_estimates(Var) = Result_1DFDC_estimates(Var) * (Mat_2DFDC_lin_t(2) - Mat_2DFDC_lin_t(1)) ;
% end
Result_y0 = Result_y0 * (Mat_2DFDC_lin_t(2) - Mat_2DFDC_lin_t(1)) ;

%%

%displayMatA(:,Ntrial)=Result_1DFDC_Mat_A;
%displayResultsDecay(:,Ntrial)=Result_1DFDC_Mat_M_Model_lin;

display(' ')
display(strcat('Total elasped time is...___', num2str(toc(Var_StartTime)), '   seconds'))
display(' ')

clear I Initial_Mat_A_Global Initial_Mat_U InMat_2DFDC_cor_dTfun InMat_2DFDC_cor_lin InMat_2DFDC_cor_log
clear InMat_2DFDC_dTfun InMat_2DFDC_lin InMat_2DFDC_log InMat_2DFDC_t K T Var1 Var2 Var3 Var_StartTime Vart
%end 