% TK_ExpMultiDeco_v02

function [TotalCurve] = TK_ExpMultiDeco_For2DFLC(Para_Tau, xdata,...
    IRF, RisePoint_FL, RisePoint_IRF, RangePoint_IRF_min, RangePoint_IRF_max)

%%
% ExpCurvePara = [0.0005, 0.1399, 0.0700, ...] ;      
% y0, A1, tau1, A2, tau2, ....
%%

I_IRF_RisePoint = RisePoint_IRF ;
I_xdata_RisePoint = RisePoint_FL ; %165;
Imax_xdata = length(xdata) ;
Idev = I_IRF_RisePoint - I_xdata_RisePoint ; % IRF rise - Signal rise

Imax = length(xdata) ;
NumOfComponent = length(Para_Tau) ;
TotalCurve = zeros(Imax, NumOfComponent) ;
Curve = zeros(Imax, 1) ;
ComNum = 0 ;
while ComNum < NumOfComponent
    ComNum = ComNum + 1 ;
    %I = 1 + 2*(ComNum-1) ;
    
    %A = ExpCurvePara(I+1) ;
    %tau = ExpCurvePara(I+2) ;
    
    Imax = RangePoint_IRF_max ; %length(IRF) - Idev ;
    I = RangePoint_IRF_min - 1;
    
    Var2 = xdata(1) ;
    Var = exp(-1 / Para_Tau(ComNum) .* (xdata-Var2)) ;
    
    Curve(:) = 0 ;
    while I < Imax
        I = I + 1 ;
        I_IRF = I ;
        
        Var1 = IRF(I_IRF) ;

        I_xdata = I - Idev ;
        Var3 = zeros(Imax_xdata, 1) ;

        Var3(I_xdata:Imax_xdata) = Var(1:Imax_xdata-I_xdata+1) ;
        Curve = Curve + Var1 .* Var3 ;
    end
    TotalCurve(:,ComNum) = Curve ;
end

Var = max(max(TotalCurve)) ;
TotalCurve = TotalCurve / Var ;


