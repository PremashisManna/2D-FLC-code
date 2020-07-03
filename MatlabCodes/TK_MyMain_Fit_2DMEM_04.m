% TK_MyMain_Fit_2DMEM_04

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tau_min = 0.05 ;
Tau_max = 5.05 ;

Tau_Lin0orLog1 = 0 ;

if Tau_Lin0orLog1 == 0
    Tau_step = 0.05 ;
elseif Tau_Lin0orLog1 == 1
    Tau_step = 100 ; % split points
end
%%%%%%%%%%%%%%%%%%%% Tau_Select = 0 ;

if Tau_Select == 0
    VarImax = length(Hozzon_estimates) ;
    Tau_Initial = zeros(VarImax, 3) ;
    VarCount = 0 ;
    VarI = -1 ;
    while VarI+2 < VarImax
        VarI = VarI + 3 ;
        VarCount = VarCount + 1 ;
        Tau_Initial(VarCount, 1) = Hozzon_estimates(VarI) ;
        Tau_Initial(VarCount, 2) = Hozzon_estimates(VarI+1) ;
        Tau_Initial(VarCount, 3) = Hozzon_estimates(VarI+2) ;
    end
    Tau_Initial = Tau_Initial(1:VarCount, :) ;
    clear VarImax VarCount VarI
%     Tau_Initial = [...%min([Hozzon_estimates(2), Hozzon_estimates(5),Hozzon_estimates(8)]), 0, 0.1 ; ...
%                    Hozzon_estimates(2), Hozzon_estimates(3), Hozzon_estimates(4); ...
%                    Hozzon_estimates(5), Hozzon_estimates(6), Hozzon_estimates(7); ...
%                    Hozzon_estimates(8), Hozzon_estimates(9), Hozzon_estimates(10)] ;    % raw = NumOfState,   column = [Amp, Tau0, GaussWidth]
%     
elseif Tau_Select == 1  % split Tau into N state at given position 
    Tau_SplitPos = [0.0, 0.2 ; ...
                    0.2, 1.1 ; ...
                    1.1, 2.3 ; ...
                    2.3, 4.2] ;
    Tau_Initial = Stored_Result_Short_Mat_A ;
    
elseif Tau_Select == 2
    Tau_StateSelect = [1,2] ;
     Tau_Initial = Gresult_Mat_A_Global ;
    %Tau_Initial = Stored_Result_Short_Mat_A ;
end

%%%%%%%%%%%%%% Fix1orNot0orAbs2orAmp3_Mat_A = [2,2,2,2] ;
%%%%%%%%%%%%%% Fix1orNot0orAbs2_Mat_G = ones(4) ;

% initial y0   y0 should be baseline intensity of 1D-FD at linear scale 
y0 = 5000;
% fit y0
Fix1orNot0orAbs2_y0 = [2] ; %was2

% For fitting parameters
UseCor1orNot0 = 0 ; % Use or not cor-matrix (after subtraction of longest-matrix to remove uncorrelated contribution)
mi_TypeSelect = 0 ; % Select mi-Type used in MEM: 0 = average constant, 1 = Initial Tau distribution
Linear0orLog1 = 1 ; % Fit linear data: 0, or log data: 1
FitStartI = 30 ;    % Start point for fitting

DisplayFig = 1 ;    % Display figure during calculation (1) or not (0)

% For fitting parameters to calculate initial Mat_A and Mat_U
BreakFactor = 100 ;
RegulatorConst = 10^(-1+0) ;
RegulatorFactor = 1.4 ;
TrialNumFor_MinimizeEstQ = 1 ;  % Trial number for minimizing Q
TrialNumFor_RegulatorConst = 100 ;    % Trial number for increasing RegulatorConst

% convoluted exp curve
RisePoint_FL = 300;
%%%%%%%%%%%%%% RisePoint_IRF = 50 ;
RangePoint_IRF_min = 50 ;
RangePoint_IRF_max = length(xdata)-70-20 ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set parameters
Var_StartTime = tic ;

% model 2D map
if Tau_Select == 0
    NumOfState = size(Tau_Initial) ;
    NumOfState = NumOfState(1) ;
elseif Tau_Select == 1
    NumOfState = size(Tau_SplitPos) ;
    NumOfState = NumOfState(1) ;
elseif Tau_Select == 2
    NumOfState = length(Tau_StateSelect) ;
end

% Set Initial Tau distribution = Initial Mat_A
if Tau_Lin0orLog1 == 0
    Tau = [Tau_min : Tau_step : Tau_max]' ;  % array of tau  /ns
elseif Tau_Lin0orLog1 == 1
    Var = (log10(Tau_max) - log10(Tau_min)) / Tau_step ;
    Tau = [log10(Tau_min) : Var : log10(Tau_max)]' ;
    Tau = 10.^Tau ;  % array of tau  /ns
end

NumOfComponent = length(Tau) ;
Tau_Initial_distribution = zeros(NumOfComponent, NumOfState) ;

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
            Tau_Initial_distribution(I, K) = VarAmp * exp(-1*((Tau(I)-VarTau0)/VarGaussWidth)^2) ; %****
        end
    end
elseif Tau_Select == 1        
    K = 0 ;
    while K < NumOfState
        K = K + 1 ;
        Vmin = Tau_SplitPos(K, 1) ;
        Vmax = Tau_SplitPos(K, 2) ;
        if (Vmin == 0) && (Vmax == 0)
            Tau_Initial_distribution(:, K) = 1 ;
        else
            I = 0 ;
            while I < NumOfComponent
                I = I + 1 ;
                if (Vmin <= Tau(I)) && (Tau(I) < Vmax)
                    Tau_Initial_distribution(I, K) = sum(Tau_Initial(I,:)) ;
                elseif Tau(I) >= Vmax
                    break
                end
            end
        end
    end
elseif Tau_Select == 2     
    K = 0 ;
    while K < NumOfState
        K = K + 1 ;
        Var = Tau_StateSelect(K) ;
        Tau_Initial_distribution(:,K) = Tau_Initial(:,Var) ;
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
    Var3 = max(max(max(Mat_2DFDC_lin))) ;
    Var4 = max(max(max(Mat_2DFDC_cor_lin))) ;
elseif Linear0orLog1 == 1
    Var2 = max(max(ExpCurve_Binned_log)) ;
    Var3 = max(max(max(Mat_2DFDC_log))) ;
    Var4 = max(max(max(Mat_2DFDC_cor_log))) ;
end
Var = (Var1^2) * (Var2^2) ;
if UseCor1orNot0 == 0 
    Var = sqrt(Var3 / Var) ;
elseif UseCor1orNot0 == 1
    Var = sqrt(Var4 / Var) ;
end
Tau_Initial_distribution = Tau_Initial_distribution .* Var ;

Var_y0 = y0 / (Mat_2DFDC_lin_t(2) - Mat_2DFDC_lin_t(1)) ;

%% Get Mat_A by fitting dT = shortest data (To estimate Initial_Mat_A)
display('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@')
display('@@@@ Calculating Mat_A for dT=shortest @@@@')
display('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@')

% fit Mat_A with fixed Mat_G = Unit matrix
%Fix1orNot0orAbs2_Mat_A = [2, 2] ;
%Fix1orNot0orAbs2_Mat_G = ones(NumOfState) ;    % Mat_G is fixed

% Get Result_Short_Mat_A
InMat_2DFDC_lin     = Mat_2DFDC_Short_lin ;
InMat_2DFDC_log     = Mat_2DFDC_Short_log ;
if UseCor1orNot0 == 1
    InMat_2DFDC_cor_lin = Mat_2DFDC_Short_cor_lin ;
    InMat_2DFDC_cor_log = Mat_2DFDC_Short_cor_log ;
elseif UseCor1orNot0 == 0
    InMat_2DFDC_cor_lin = Mat_2DFDC_Short_lin ;
    InMat_2DFDC_cor_log = Mat_2DFDC_Short_log ;
end
    
[Result_Short_estimates, Vart, Result_Short_Mat_M_2DFLC,...
    Result_Short_Mat_M_Model_lin, Result_Short_Mat_M_Model_log,...
    Result_Short_Mat_A, Result_Short_Mat_G, Result_Short_y0, Result_Short_EstQ] =...
    TK_FitF_MinimizeQ_09(...
    Mat_2DFDC_lin_t, InMat_2DFDC_lin, InMat_2DFDC_cor_lin,...
    Mat_2DFDC_log_t, InMat_2DFDC_log, InMat_2DFDC_cor_log,...
    Tau, Tau_Initial_distribution, Fix1orNot0orAbs2orAmp3_Mat_A, Fix1orNot0orAbs2orSym3_Mat_G, Var_y0, Fix1orNot0orAbs2_y0,...
    ExpCurve_Binned_lin, ExpCurve_Binned_log,...
    NumOfState, RegulatorConst, RegulatorFactor,...
    Linear0orLog1, FitStartI, TrialNumFor_MinimizeEstQ, TrialNumFor_RegulatorConst, mi_TypeSelect, BreakFactor,...
    tMin, tMax, DisplayFig) ;

% if Fix1orNot0orAbs2_y0 ~= 1
%     Var = length(Result_Short_estimates) ;
%     Result_Short_estimates(Var) = Result_Short_estimates(Var) * (Mat_2DFDC_lin_t(2) - Mat_2DFDC_lin_t(1)) ;
% end
Result_Short_y0 = Result_Short_y0 * (Mat_2DFDC_lin_t(2) - Mat_2DFDC_lin_t(1)) ;

%%

display(' ')
display(strcat('Total elasped time is...___', num2str(toc(Var_StartTime)), '   seconds'))
display(' ')

clear I Initial_Mat_A_Global InMat_2DFDC_cor_dTfun InMat_2DFDC_cor_lin InMat_2DFDC_cor_log
clear InMat_2DFDC_dTfun InMat_2DFDC_lin InMat_2DFDC_log InMat_2DFDC_t K T Var1 Var2 Var3 Var_StartTime Vart
