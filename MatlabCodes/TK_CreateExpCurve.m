% TK_CreateExpCurve

function [ExpCurve, ExpCurve_Binned_lin, ExpCurve_Binned_log] = TK_CreateExpCurve(...
    Tau, xdata, IRF, NumOfComponent,...
    tMin, tMax, tStep, lint_BinFactor, Mat_2DFDC_log_t,...
RisePoint_FL, RisePoint_IRF, RangePoint_IRF_min, RangePoint_IRF_max)

%% convoluted exp curve

[ExpCurve] = TK_ExpMultiDeco_For2DFLC(Tau, xdata,...
    IRF, RisePoint_FL, RisePoint_IRF, RangePoint_IRF_min, RangePoint_IRF_max) ;

Var1 = round(tMin/tStep) ;
Var2 = round(tMax/tStep) ;
ExpCurve = ExpCurve(Var1:Var2, :) ;
Var_xdata = xdata(Var1:Var2) ;
Var_xdata = Var_xdata - Var_xdata(1) ;

% linearly binning ExpCurve
Imax = ceil(Var2-Var1) ;
Var3 = ceil( (Imax)/lint_BinFactor ) ;
ExpCurve_Binned_lin = ExpCurve(1:Var3, :) ;
I = 0 ;
while I < Var3
    I = I + 1 ;
    Var1 = 1+lint_BinFactor*(I-1) ;
    Var2 = lint_BinFactor*(I) ;
    if Var2 <= Imax
        ExpCurve_Binned_lin(I, :) = sum(ExpCurve(Var1:Var2, :)) ;
    else
        ExpCurve_Binned_lin(I, :) = sum(ExpCurve(Var1:Imax, :)) + sum(ExpCurve(Imax-(Var2-Imax)+1:Imax, :)) ;
    end
end

% logarithmically binning ExpCurve
Var = size(ExpCurve) ;
Imax = Var(1) ;
Kmax = length(Mat_2DFDC_log_t) ;
ExpCurve_Binned_log = zeros(Kmax, NumOfComponent) ;
I = 0 ;
K = 2 ;
while I < Imax 
    I = I + 1 ;
    if K <= Kmax
        if Var_xdata(I) <= Mat_2DFDC_log_t(K)
            ExpCurve_Binned_log(K-1, :) = ExpCurve_Binned_log(K-1, :) + ExpCurve(I, :) ;
        else
            K = K + 1 ;
            I = I - 1 ;
        end
    else
        ExpCurve_Binned_log(K-1, :) = ExpCurve_Binned_log(K-1, :) + ExpCurve(I, :) ;
    end
end