%TK_FitF_MinimizeQ_09

function [estimates, Result_Mat_2DFDC_t, Result_Mat_M_2DFLC, Result_Mat_M_Model_lin, Result_Mat_M_Model_log,...
    Result_Mat_A, Result_Mat_G, Result_y0, Result_EstQ] =...
    TK_FitF_MinimizeQ_09(...
    Mat_2DFDC_lint, Mat_2DFDC_lin, Mat_2DFDC_cor_lin,...
    Mat_2DFDC_logt, Mat_2DFDC_log, Mat_2DFDC_cor_log,...
    Tau, Tau_Initial_distribution, Fix1orNot0orAbs2orAmp3_Mat_A, Fix1orNot0orAbs2orSym3_Mat_G, Initial_y0, Fix1orNot0orAbs2_y0,...
    ExpCurve_Binned_lin, ExpCurve_Binned_log,...
    NumOfState, RegulatorConst, RegulatorFactor,...
    Linear0orLog1, FitStartI, TrialNumFor_MinimizeEstQ, TrialNumFor_RegulatorConst, mi_TypeSelect, BreakFactor,...
    tMin, tMax, DisplayFig)

%%

% tMin = 1;%1.9    ;   % ns
% tMax = 12 ;   % ns
% tStep = 0.004 ; % ns
% 
% lint_BinFactor = 4 ;   % data points in linear scale = t_Imax / lint_BinFactor
% 
% % convoluted exp curve
% RisePoint_FL = 50 ;
% RisePoint_IRF = 50 ;
% RangePoint_IRF_min = 100 ;
% RangePoint_IRF_max = 1500 ;
% 
% Tau_Initial = [1.5 ; 3.5]  ; % raw = NumOfState,  column = Tau distribution
% Tau_min = 0.1 ;
% Tau_max = 6.1 ;
% Tau_step = 0.1 ;
% 
% Fix1orNot0orAbs2 = [2, 1] ; % column 1 for Mat_A, column 2 for Mat_U
% 
% % model 2D map
% NumOfState = 2 ;
% RegulatorConst = 0.1 ;
% 
% % Fitting
% Linear0orLog1 = 1 ;
% FitStartI = 30 ;
% TrialNumFor_MinimizeEstQ = 30 ;
% TrialNumFor_RegulatorConst = 100 ;


%% Set Tau distribution = Mat_A
Initial_Mat_A = Tau_Initial_distribution ;
Initial_Mat_G = eye(NumOfState, NumOfState) ;

%% Set Matrix for fitting
% adjust lin data length
if Linear0orLog1 == 0   
    Imax = length(Mat_2DFDC_lint) ;
    Result_Mat_2DFDC_t = Mat_2DFDC_lint(FitStartI:Imax) ;
    In_Mat_2DFDC = Mat_2DFDC_lin(FitStartI:Imax, FitStartI:Imax) ;
    In_Mat_2DFDC_cor = Mat_2DFDC_cor_lin(FitStartI:Imax, FitStartI:Imax) ;
    In_ExpCurve = ExpCurve_Binned_lin(FitStartI:Imax, :) ;
else
    Imax = length(Mat_2DFDC_logt) ;
    Result_Mat_2DFDC_t = Mat_2DFDC_logt(FitStartI:Imax) ;
    In_Mat_2DFDC = Mat_2DFDC_log(FitStartI:Imax, FitStartI:Imax) ;
    In_Mat_2DFDC_cor = Mat_2DFDC_cor_log(FitStartI:Imax, FitStartI:Imax) ;
    In_ExpCurve = ExpCurve_Binned_log(FitStartI:Imax, :) ;
end

%% Minimize Q value by increasing RegulatorConst

tic

Var1 = Initial_Mat_A ;
Var2 = Initial_Mat_G ;
Var5 = Initial_y0 ;
Result_EstQ = zeros(TrialNumFor_RegulatorConst*TrialNumFor_MinimizeEstQ, 2) ;
Count = 0 ;
I = 0 ;
while I < TrialNumFor_RegulatorConst
    I = I + 1 ;
    Var3 = RegulatorConst * (RegulatorFactor^(I-1)) ;
    
    % mi_func calculation
    [mi] = TK_mi_ModelFunction(Var1, Tau, tMax, tMin, NumOfState, mi_TypeSelect) ;

    Check = 0 ;
    K = 0 ;
    while K < TrialNumFor_MinimizeEstQ
        K = K + 1 ;
        Count = Count + 1 ;
       
        [estimates, model] = TK_FitF_2DMEM_07(...
            Var1, Var2, Var5, Fix1orNot0orAbs2orAmp3_Mat_A, Fix1orNot0orAbs2orSym3_Mat_G, Fix1orNot0orAbs2_y0, Var3, ...
            Result_Mat_2DFDC_t, In_Mat_2DFDC, In_Mat_2DFDC_cor, ...
            NumOfState, Tau, In_ExpCurve, mi) ;
        
        [Result_Mat_M_2DFLC, Result_Mat_M_Model, Result_Mat_A, Result_Mat_G, Result_y0, Var_EstQ, Var_Kai2, Var_EntropyS] = TK_FitF_Reproduct2DFDCand2DFLC_03(...
            estimates, Var1, Var2, Var5, Fix1orNot0orAbs2orAmp3_Mat_A, Fix1orNot0orAbs2orSym3_Mat_G, Fix1orNot0orAbs2_y0, Var3, ...
            Result_Mat_2DFDC_t, In_Mat_2DFDC, In_Mat_2DFDC_cor, ...
            NumOfState, Tau, In_ExpCurve, mi) ;
             
        Var1 = Result_Mat_A ;
        Var2 = Result_Mat_G ;
        Var5 = Result_y0 ;
        
        Var = Check - Var_EstQ ;
        Check = Var_EstQ ;
        Result_EstQ(Count, :) = [Var3, Var_EstQ] ;
        
        toc
        %Var4 = abs(Var) <= (abs(Var_EstQ) * 0.01)
        if (abs(Var) <= (abs(Var_EstQ) * 10^-10))
            break
        end
    end

    display(strcat('TrialNumFor_RegulatorConst =', ' ', num2str(I), ' /', num2str(TrialNumFor_RegulatorConst)))
    display(strcat('Kai2 / (2*EntropyS/RegulatorConst) =', num2str(abs(Var_Kai2/(2*Var_EntropyS/Var3))) ) )
    display(' ')
    
    
    if DisplayFig == 1
        [Plot_Var1, Plot_Var2] = max(max(Result_Mat_M_Model(:,:,1))) ;
        Var = size(Result_Mat_A) ;
        if I == 1
            figure; Plot_hp1 = plot(Tau, Result_Mat_A) ;
            T = 0 ;
            while T < Var(2)
                T = T + 1 ;
                set(Plot_hp1(T),'YDataSource', strcat('Result_Mat_A','(:,',num2str(T),')')) ;
            end 
            
            Plot_X = Result_Mat_2DFDC_t ;
            Plot_Y1 = In_Mat_2DFDC_cor(Plot_Var2,:,1) ;
            Plot_Y2 = Result_Mat_M_Model(Plot_Var2,:,1) ;
            figure; Plot_hp2 = plot(Plot_X,Plot_Y1, Plot_X,Plot_Y2) ;
            %set(Plot_hp2(1),'YDataSource', 'Plot_Y1') ;
            set(Plot_hp2(2),'YDataSource', 'Plot_Y2') ;
        end
        T = 0 ;
        while T < Var(2)
            T = T + 1 ;
            refreshdata(Plot_hp1(T), 'caller')
        end
        
        %refreshdata(Plot_hp2(1), 'caller')
        Plot_Y2 = Result_Mat_M_Model(Plot_Var2,:,1) ;
        refreshdata(Plot_hp2(2), 'caller')
        drawnow
    end
    
    if (abs(Var_Kai2/(2*Var_EntropyS/Var3)) > 10^BreakFactor)
        display('break:  Kai2 / (2*EntropyS/RegulatorConst) > 10^BreakFactor')
        break
    end
end

Result_EstQ = Result_EstQ(1:Count, :) ;

%% Create Mat_M_Model without binning for both linear and log scale.
[Result_Mat_M_2DFLC, Result_Mat_M_Model_lin, Result_Mat_A, Result_Mat_G, Result_y0, Var_EstQ, Var_Kai2, Var_EntropyS] = TK_FitF_Reproduct2DFDCand2DFLC_03(...
    estimates, Var1, Var2, Var5, Fix1orNot0orAbs2orAmp3_Mat_A, Fix1orNot0orAbs2orSym3_Mat_G, Fix1orNot0orAbs2_y0, Var3, ...
    Mat_2DFDC_lint, Mat_2DFDC_lin, Mat_2DFDC_cor_lin, ...
    NumOfState, Tau, ExpCurve_Binned_lin, mi) ;

[Result_Mat_M_2DFLC, Result_Mat_M_Model_log, Result_Mat_A, Result_Mat_G, Result_y0, Var_EstQ, Var_Kai2, Var_EntropyS] = TK_FitF_Reproduct2DFDCand2DFLC_03(...
    estimates, Var1, Var2, Var5, Fix1orNot0orAbs2orAmp3_Mat_A, Fix1orNot0orAbs2orSym3_Mat_G, Fix1orNot0orAbs2_y0, Var3, ...
    Mat_2DFDC_logt, Mat_2DFDC_log, Mat_2DFDC_cor_log, ...
    NumOfState, Tau, ExpCurve_Binned_log, mi) ;

end
