% TK_MyMain_Fit_GaussianMulti

%%
TrialNum = 20 ;
SSE1D = zeros(1, TrialNum) ;

Var_xdata = Tau ;
Var_ydata = A ;

NumOfComponent = 4 ;
switch NumOfComponent
    case 1
        InitialPara =   [0, ...
                        1, 1, 1] ;      % Initial parameters  % y0, A, x0, dx
        Fix1orNot0 =    [1, ...
                        0, 0, 0] ;        % Fix para => 1
    case 2
        InitialPara =   [0, ...
                        1, 2, 1, ...
                        0.2, 8, 3] ;     % Initial parameters  % y0, A, x0, dx
        Fix1orNot0 =    [1, ...
                        0, 0, 0, ...
                        0, 0, 0] ;       % Fix para => 1
    case 3
        InitialPara =   [0, ...
                        1, 0.4, 0.2, ...
                        1, 1.8, 0.2, ...
                        1, 2.8, 0.2] ;     % Initial parameters  % y0, A, x0, dx
        InitialPara =   [Hozzon_estimates(1), ...
                        Hozzon_estimates(2), Hozzon_estimates(3), Hozzon_estimates(4), ...
                        Hozzon_estimates(5), Hozzon_estimates(6), Hozzon_estimates(7), ...
                        Hozzon_estimates(8), Hozzon_estimates(9), Hozzon_estimates(10)] ;     % Initial parameters  % y0, A, x0, dx
                    
        Fix1orNot0 =    [1, ...
                        0, 1, 0, ...
                        0, 0, 0, ...
                        0, 0, 0] ;       % Fix para => 1
    case 4
        InitialPara =   [0, ...
                        1, 0.4, 0.2, ...
                        1, 1.0, 0.2, ...
                        1, 1.8, 0.2, ...
                        1, 2.8, 0.2] ;     % Initial parameters  % y0, A, x0, dx
        InitialPara =   [Hozzon_estimates(1), ...
                        Hozzon_estimates(2), Hozzon_estimates(3), Hozzon_estimates(4), ...
                        Hozzon_estimates(5), Hozzon_estimates(6), Hozzon_estimates(7), ...
                        Hozzon_estimates(8), Hozzon_estimates(9), Hozzon_estimates(10), ...
                        Hozzon_estimates(11), Hozzon_estimates(12), Hozzon_estimates(13)] ;     % Initial parameters  % y0, A, x0, dx
                    
        Fix1orNot0 =    [1, ...
                        1, 1, 1, ...
                        1, 1, 1, ...
                        1, 1, 1, ...
                        1, 1, 1] ;       % Fix para => 1
end
%FitRange = [-inf inf] ;        % fittng range / points     [-inf, inf] -> whole range
FitRange = [1 inf] ;


%%

if isinf(abs(FitRange(1,1))) == 1
    Imin = 1 ;
else
    Imin = FitRange(1,1) ;
end
if isinf(abs(FitRange(1,2))) == 1
    Imax = length(Var_xdata) ;
else
    Imax = FitRange(1,2) ;
end
Var_xdata = Var_xdata(Imin:Imax) ;
Var_ydata = Var_ydata(Imin:Imax) ;
clear Imin Imax ;


I = 0 ;
while I < TrialNum
    I = I + 1 ;
    
    [estimates, model] = TK_FitF_GaussianMulti(Var_xdata, Var_ydata, NumOfComponent, InitialPara, Fix1orNot0) ;
    
    [sse, FittedCurve] = model(estimates);
    
    SSE1D(1, I) = sse ;
    InitialPara = estimates ;
end

figure ; plot(Var_xdata, Var_ydata, Var_xdata, FittedCurve, 'r') ;
estimates ;

switch NumOfComponent
    case 1
        disp(['A = ', num2str(estimates(1,2))])
        disp(['x0 = ', num2str(estimates(1,3))])
        disp(['dx = ', num2str(estimates(1,4))])
    case 2
        disp(['A1 = ', num2str(estimates(1,2))])
        disp(['x0_1 = ', num2str(estimates(1,3))])
        disp(['dx_1 = ', num2str(estimates(1,4))])
        disp(' ')
        disp(['A2 = ', num2str(estimates(1,5))])
        disp(['x0_2 = ', num2str(estimates(1,6))])
        disp(['dx_2 = ', num2str(estimates(1,7))])
    case 3
        disp(['A1 = ', num2str(estimates(1,2))])
        disp(['x0_1 = ', num2str(estimates(1,3))])
        disp(['dx_1 = ', num2str(estimates(1,4))])
        disp(['A2 = ', num2str(estimates(1,5))])
        disp(['x0_2 = ', num2str(estimates(1,6))])
        disp(['dx_2 = ', num2str(estimates(1,7))])
        disp(['A3 = ', num2str(estimates(1,8))])
        disp(['x0_3 = ', num2str(estimates(1,9))])
        disp(['dx_3 = ', num2str(estimates(1,10))])
    case 4
        disp(['A1 = ', num2str(estimates(1,2))])
        disp(['x0_1 = ', num2str(estimates(1,3))])
        disp(['dx_1 = ', num2str(estimates(1,4))])
        disp(['A2 = ', num2str(estimates(1,5))])
        disp(['x0_2 = ', num2str(estimates(1,6))])
        disp(['dx_2 = ', num2str(estimates(1,7))])
        disp(['A3 = ', num2str(estimates(1,8))])
        disp(['x0_3 = ', num2str(estimates(1,9))])
        disp(['dx_3 = ', num2str(estimates(1,10))])
        disp(['A4 = ', num2str(estimates(1,11))])
        disp(['x0_4 = ', num2str(estimates(1,12))])
        disp(['dx_4 = ', num2str(estimates(1,13))])
end

clear TrialNum I NumOfComponent Var_ydata
%xlabel('xdata')
%ylabel('f(estimates,xdata)')
%title(['Fitting to function ', func2str(model)]);
%legend('data', ['fit using ', func2str(model)])
%hold off