function [estimates, model] = TK_FitF_GaussianMulti(xdata, ydata, NumOfComponent, InitialPara, Fix1orNot0)
% InitialPara: y0, A, x0, dx .....
% Num of vars: 1 + 3*NumOfComponent = y0, (A, x0, dx)*NumOfComponent

% Call fminsearch with a random starting point.
model = @expfun ;
estimates = fminsearch(model, InitialPara) ;
% expfun accepts curve parameters as inputs, and outputs sse,
% the sum of squares error for A*exp(-lambda*xdata)-ydata,
% and the FittedCurve. FMINSEARCH only needs sse, but we want
% to plot the FittedCurve at the end.
    function [sse, FittedCurve] = expfun(params)
        if (Fix1orNot0(1) == 1)
            y0 = InitialPara(1) ;
        else
            y0 = params(1) ;
        end
        
        FittedCurve = y0 ;
        
        ComNum = 0 ;
        while ComNum < NumOfComponent
            ComNum = ComNum + 1 ;
            I = 1 + 3*(ComNum-1) ;
            if (Fix1orNot0(I+1) == 1)
                A = InitialPara(I+1) ;
            else
                A = params(I+1) ;
            end
            if (Fix1orNot0(I+2) == 1)
                x0 = InitialPara(I+2) ;
            else
                x0 = params(I+2) ;
            end
            if (Fix1orNot0(I+3) == 1)
                dx = InitialPara(I+3) ;
            else
                dx = params(I+3) ;
            end
            
            Var = (xdata - x0)./ dx ;
            FittedCurve = FittedCurve + A .* exp(-1 .* Var.^2) ;
        end
        ErrorVector = FittedCurve - ydata;
        sse = sum(ErrorVector .^ 2/length(ydata));
    end

Imax = 1 + ComNum*3 ;
I = 0 ;
while I < Imax
    I = I + 1 ;
    if (Fix1orNot0(I) == 1)
        estimates(I) = InitialPara(I) ;
    end
end

end


