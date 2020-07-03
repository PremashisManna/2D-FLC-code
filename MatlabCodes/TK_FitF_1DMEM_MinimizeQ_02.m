%TK_FitF_1DMEM_MinimizeQ_02

function [estimates, Result_Mat_1DFDC_t, Result_Mat_M_Model_lin, Result_Mat_M_Model_log,...
    Result_Mat_A, Result_y0, Result_EstQ] =...
    TK_FitF_1DMEM_MinimizeQ_02(...
    Mat_1DFDC_lint, Mat_1DFDC_lin, Mat_1DFDC_cor_lin,...
    Mat_1DFDC_logt, Mat_1DFDC_log, Mat_1DFDC_cor_log,...
    Tau, Tau_Initial_distribution, Fix1orNot0orAbs2, Initial_y0, Fix1orNot0orAbs2_y0,...
    ExpCurve_Binned_lin, ExpCurve_Binned_log,...
    RegulatorConst, RegulatorFactor,...
    Linear0orLog1, FitStartI, TrialNumFor_MinimizeEstQ, TrialNumFor_RegulatorConst, mi_TypeSelect, BreakFactor,...
    tMin, tMax, DisplayFig)

%%
NumOfState = 1 ;
%% Set Tau distribution = Mat_A
Initial_Mat_A = Tau_Initial_distribution ;

%% Set Matrix for fitting
% adjust lin data length
if Linear0orLog1 == 0   
    Imax = length(Mat_1DFDC_lint) ;
    Result_Mat_1DFDC_t = Mat_1DFDC_lint(FitStartI:Imax) ;
    In_Mat_1DFDC = Mat_1DFDC_lin(FitStartI:Imax) ;
    In_Mat_1DFDC_cor = Mat_1DFDC_cor_lin(FitStartI:Imax) ;
    In_ExpCurve = ExpCurve_Binned_lin(FitStartI:Imax, :) ;
else
    Imax = length(Mat_1DFDC_logt) ;
    Result_Mat_1DFDC_t = Mat_1DFDC_logt(FitStartI:Imax) ;
    In_Mat_1DFDC = Mat_1DFDC_log(FitStartI:Imax) ;
    In_Mat_1DFDC_cor = Mat_1DFDC_cor_log(FitStartI:Imax) ;
    In_ExpCurve = ExpCurve_Binned_log(FitStartI:Imax, :) ;
end

%% Minimize Q value by increasing RegulatorConst

tic

Var1 = Initial_Mat_A ;
Var2 = Initial_y0 ;

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

        [estimates, model] = TK_FitF_1DMEM_02(...
            Var1, Fix1orNot0orAbs2, Var2, Fix1orNot0orAbs2_y0, Var3, ...
            Result_Mat_1DFDC_t, In_Mat_1DFDC, In_Mat_1DFDC_cor, ...
            Tau, In_ExpCurve, mi) ;
        
        [Result_Mat_M_Model, Result_Mat_A, Result_y0, Var_EstQ, Var_Kai2, Var_EntropyS] = TK_FitF_Reproduct1DFDC_02(...
            estimates, Var1, Fix1orNot0orAbs2, Var2, Fix1orNot0orAbs2_y0, Var3, ...
            Result_Mat_1DFDC_t, In_Mat_1DFDC, In_Mat_1DFDC_cor, ...
            Tau, In_ExpCurve, mi) ;

        Var1 = Result_Mat_A ;
        
%         Var1 = Result_Mat_A .*( (1 - exp(-1*(tMax-tMin) ./ Tau)) ).* 1;
%        Var1 = Result_Mat_A .* Initial_Mat_A ./ max(max(Initial_Mat_A)) ;      
        Var2 = Result_y0 ;

        Var = Check - Var_EstQ ;
        Check = Var_EstQ ;
        Result_EstQ(Count, :) = [Var3, Var_EstQ] ;
        
        toc
        
        if (abs(Var) <= (abs(Var_EstQ) * 10^-5))
            break
        end
    end
    
    display(strcat('TrialNumFor_RegulatorConst =', ' ', num2str(I), ' /', num2str(TrialNumFor_RegulatorConst)))
    display(strcat('Kai2 / (2*EntropyS/RegulatorConst) =', num2str(abs(Var_Kai2/(2*Var_EntropyS/Var3))) ) )
    display(' ')
    
%     if DisplayFig == 1
%         if I == 1
%             figure; plot(Tau, Result_Mat_A) ;
%             drawnow
%         else
%             plot(Tau, Result_Mat_A) ;
%             drawnow
%         end
%     end
    if DisplayFig == 1
        %[Plot_Var1, Plot_Var2] = max(max(Result_Mat_M_Model) ;
        Var = size(Result_Mat_A) ;
        
        
        
        if I == 1
            figure; Plot_hp1 = plot(Tau, Result_Mat_A) ;
            set(gcf,'color','w')
            xlabel('\tau (ns)')
            ylabel('a(\tau)')
            set(gca,'fontSize',20)
            T = 0 ;
            while T < Var(2)
                T = T + 1 ;
                set(Plot_hp1(T),'YDataSource', strcat('Result_Mat_A','(:,',num2str(T),')')) ;
            end 
            
            
            Plot_X = Result_Mat_1DFDC_t ;
            Plot_Y1 = In_Mat_1DFDC_cor ;
            Plot_Y2 = Result_Mat_M_Model ;
            figure; Plot_hp2 = plot(Plot_X,Plot_Y1, Plot_X,Plot_Y2) ;
            %set(Plot_hp2(1),'YDataSource', 'Plot_Y1') ;
            set(Plot_hp2(2),'YDataSource', 'Plot_Y2') ;
            set(gcf,'color','w')
            xlabel('t (ns)')
            ylabel('counts')
            set(gca,'fontSize',20)
            
        end
        T = 0 ;
        while T < Var(2)
            T = T + 1 ;
            refreshdata(Plot_hp1(T), 'caller')
        end
        
        %refreshdata(Plot_hp2(1), 'caller')
        Plot_Y2 = Result_Mat_M_Model ;
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

[Result_Mat_M_Model_lin, Result_Mat_A, Result_y0, Var_EstQ, Var_Kai2, Var_EntropyS] = TK_FitF_Reproduct1DFDC_02(...
    estimates, Var1, Fix1orNot0orAbs2, Var2, Fix1orNot0orAbs2_y0, Var3, ...
    Mat_1DFDC_lint, Mat_1DFDC_lin, Mat_1DFDC_cor_lin, ...
    Tau, ExpCurve_Binned_lin, mi) ;

[Result_Mat_M_Model_log, Result_Mat_A, Result_y0, Var_EstQ, Var_Kai2, Var_EntropyS] = TK_FitF_Reproduct1DFDC_02(...
    estimates, Var1, Fix1orNot0orAbs2, Var2, Fix1orNot0orAbs2_y0, Var3, ...
    Mat_1DFDC_logt, Mat_1DFDC_log, Mat_1DFDC_cor_log, ...
    Tau, ExpCurve_Binned_log, mi) ;

end
