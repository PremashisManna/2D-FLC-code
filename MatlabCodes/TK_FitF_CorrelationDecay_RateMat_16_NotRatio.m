% TK_FitF_CorrelationDecay_RateMat_16_NotRatio

function [Ana_estimates, Ana_Free0orFix1orAbs2, Ana_estimates_A, Ana_estimates_y0,...
    Ana_estimates_E, Ana_estimates_Q, Ana_estimates_RateM, Ana_FittedX, Ana_FittedY, model] = ...
    TK_FitF_CorrelationDecay_RateMat_16_NotRatio(...
    xdata, ydata, ComponentNumber, StateAssign,...
    FitA, Fity0, FitE, FitQ, FitRateM,...
    Free0orFix1orAbs2_A, Free0orFix1orAbs2_y0, Free0orFix1orAbs2_E, Free0orFix1orAbs2_Q, Free0orFix1orAbs2_RateM,...
    FittedPoints, FitNormalizedCor_y1n0, FitWeightFactor)

ComNumMax = ComponentNumber ;
Kmax = zeros(ComNumMax, 1) ;
VarN = zeros(ComNumMax, 1) ;
Var1 = size(FitRateM) ;
Var1 = Var1(1) ;
Var2 = max(StateAssign) ;
Var2 = Var2*(Var2-1)/2 ;
Imax = 1 + Var2*3 + Var1 + Var1 + Var1*Var1 ;
InitialPara = zeros(Imax*ComNumMax, 1) ;

Global_Check_E = Free0orFix1orAbs2_E .* 0 ;
Global_Para_E = zeros(max(StateAssign)*ComNumMax, 1) ;
Global_Check_Q = Free0orFix1orAbs2_Q .* 0 ;
Global_Para_Q = zeros(max(StateAssign)*ComNumMax, 1) ;
Global_Check_RateM = Free0orFix1orAbs2_RateM .* 0 ;
Global_Para_RateM = zeros(Imax*ComNumMax, 1) ;
Global_Check_y0 = Free0orFix1orAbs2_y0 .* 0 ;
Global_Para_y0 = zeros(max(StateAssign)*max(StateAssign), 1) ;

Count = 0 ;

if (Free0orFix1orAbs2_A ~= 1)
    Count = Count + 1 ;
    if (Free0orFix1orAbs2_A == 0)
        InitialPara(Count) = FitA ;
    elseif (Free0orFix1orAbs2_A == 2)
        InitialPara(Count) = abs(FitA) ;
    end
end

ComNum = 0 ;
while ComNum < ComNumMax
    ComNum = ComNum + 1 ;
    
    Kmax(ComNum) = length(StateAssign) ;
    VarN(ComNum) = max(StateAssign) ;
       
    K = 0 ;
    while K < Kmax(ComNum)
        K = K + 1 ;
        
        if (Free0orFix1orAbs2_E(ComNum, K) ~= 1)
            if Free0orFix1orAbs2_E(ComNum, K) >= 3 && Global_Check_E(ComNum, K) == 0
                Count = Count + 1 ;
                InitialPara(Count) = abs(FitE(ComNum, K)) ;
                Global_Check_E(:,:) = Global_Check_E(:,:) + (Free0orFix1orAbs2_E==Free0orFix1orAbs2_E(ComNum, K))...
                    *Free0orFix1orAbs2_E(ComNum, K) ;
            elseif Free0orFix1orAbs2_E(ComNum, K) >= 3 && Global_Check_E(ComNum, K) ~= 0
                
            else
                Count = Count + 1 ;
                if (Free0orFix1orAbs2_E(ComNum, K) == 0)
                    InitialPara(Count) = FitE(ComNum, K) ;
                elseif (Free0orFix1orAbs2_E(ComNum, K) == 2)
                    InitialPara(Count) = abs(FitE(ComNum, K)) ;
                end
            end
        end
        
        
        if (Free0orFix1orAbs2_Q(ComNum, K) ~= 1)
            if Free0orFix1orAbs2_Q(ComNum, K) >= 3 && Global_Check_Q(ComNum, K) == 0
                Count = Count + 1 ;
                InitialPara(Count) = abs(FitQ(ComNum, K)) ;
                Global_Check_Q(:,:) = Global_Check_Q(:,:) + (Free0orFix1orAbs2_Q==Free0orFix1orAbs2_Q(ComNum, K))...
                    *Free0orFix1orAbs2_Q(ComNum, K) ;
            elseif Free0orFix1orAbs2_Q(ComNum, K) >= 3 && Global_Check_Q(ComNum, K) ~= 0
                
            else
                Count = Count + 1 ;
                if (Free0orFix1orAbs2_Q(ComNum, K) == 0)
                    InitialPara(Count) = FitQ(ComNum, K) ;
                elseif (Free0orFix1orAbs2_Q(ComNum, K) == 2)
                    InitialPara(Count) = abs(FitQ(ComNum, K)) ;
                end
            end
        end
        
        
        I = 0 ;
        while I < Kmax(ComNum)
            I = I + 1 ;
            if (Free0orFix1orAbs2_RateM(I, K, ComNum) ~= 1)
                if Free0orFix1orAbs2_RateM(I, K, ComNum) >= 3 && Global_Check_RateM(I, K, ComNum) == 0
                    Count = Count + 1 ;
                    InitialPara(Count) = abs(FitRateM(I, K, ComNum)) ;
                    Global_Check_RateM(:,:,:) = Global_Check_RateM(:,:,:) + (Free0orFix1orAbs2_RateM==Free0orFix1orAbs2_RateM(I, K, ComNum))...
                                                                                  *Free0orFix1orAbs2_RateM(I, K, ComNum) ;
                elseif Free0orFix1orAbs2_RateM(I, K, ComNum) >= 3 && Global_Check_RateM(I, K, ComNum) ~= 0
                    
                else
                    Count = Count + 1 ;
                    if (Free0orFix1orAbs2_RateM(I, K, ComNum) == 0)
                        InitialPara(Count) = FitRateM(I, K, ComNum) ;
                    elseif (Free0orFix1orAbs2_RateM(I, K, ComNum) == 2)
                        InitialPara(Count) = abs(FitRateM(I, K, ComNum)) ;
                    end
                end
            end
        end
    end
    
end

Kmaxmax = max(Kmax) ;
VarNmax = max(VarN) ;

K = 0 ;
while K < VarNmax
    K = K + 1 ;
    I = 0 ;
    while I < VarNmax
        I = I + 1 ;

        if (Free0orFix1orAbs2_y0(I, K) ~= 1)
            if Free0orFix1orAbs2_y0(I, K) >= 3 && Global_Check_y0(I, K) == 0
                Count = Count + 1 ;
                InitialPara(Count) = abs(Fity0(I, K)) ;
                Global_Check_y0(:,:) = Global_Check_y0(:,:) + (Free0orFix1orAbs2_y0==Free0orFix1orAbs2_y0(I, K))...
                    *Free0orFix1orAbs2_y0(I, K) ;
            elseif Free0orFix1orAbs2_y0(I, K) >= 3 && Global_Check_y0(I, K) ~= 0
                
            else
                Count = Count + 1 ;
                if (Free0orFix1orAbs2_y0(I, K) == 0)
                    InitialPara(Count) = Fity0(I, K) ;
                elseif (Free0orFix1orAbs2_y0(I, K) == 2)
                    InitialPara(Count) = abs(Fity0(I, K)) ;
                end
            end
        end
    end
end

InitialPara = InitialPara(1:Count) ;

model = @expfun ;
options = optimset('MaxFunEvals',10^5,'MaxIter',10^5);
estimates = fminsearch(model, InitialPara, options) ;

    function [sse, FittedCurve] = expfun(params)
        %Amp = zeros(1) ;
        y0 = zeros(VarNmax, VarNmax) ;
        Eps = zeros(ComNumMax, Kmaxmax) ;
        Qua = zeros(ComNumMax, Kmaxmax) ;
        RateM = zeros(Kmaxmax, Kmaxmax, ComNumMax) ;
        %Pop_V = zeros(length(xdata), Kmax) ;
            
        Count = 0 ;
        
        if (Free0orFix1orAbs2_A == 1)
            Amp = FitA ;
        elseif (Free0orFix1orAbs2_A == 0)
            Count = Count + 1 ;
            Amp = params(Count) ;
        elseif (Free0orFix1orAbs2_A == 2)
            Count = Count + 1 ;
            Amp = abs(params(Count)) ;
        end
        
        Global_Check_E(:,:) = 0 ;
        Global_Check_Q(:,:) = 0 ;
        Global_Check_RateM(:,:,:) = 0 ;
        Global_Check_y0(:,:) = 0 ;
        
        ComNum = 0 ;
        while ComNum < ComNumMax
            ComNum = ComNum + 1 ;
                     
            K = 0 ;
            while K < Kmax(ComNum)
                K = K + 1 ;
               
                if Free0orFix1orAbs2_E(ComNum, K) >= 3 && Global_Check_E(ComNum, K) == 0
                    Count = Count + 1 ;
                    Eps(ComNum, K) = abs(params(Count)) ;
                    Var = Global_Check_E(:,:) + (Free0orFix1orAbs2_E==Free0orFix1orAbs2_E(ComNum, K))...
                        *Free0orFix1orAbs2_E(ComNum, K) ;
                    Global_Check_E = Var ;
                    Global_Para_E(Free0orFix1orAbs2_E(ComNum, K)) = Eps(ComNum, K) ;
                elseif Free0orFix1orAbs2_E(ComNum, K) >= 3 && Global_Check_E(ComNum, K) ~= 0
                    Eps(ComNum, K) = Global_Para_E(Free0orFix1orAbs2_E(ComNum, K)) ;
                else
                    if (Free0orFix1orAbs2_E(ComNum, K) == 1)
                        Eps(ComNum, K) = FitE(ComNum, K) ;
                    elseif (Free0orFix1orAbs2_E(ComNum, K) == 0)
                        Count = Count + 1 ;
                        Eps(ComNum, K) = params(Count) ;
                    elseif (Free0orFix1orAbs2_E(ComNum, K) == 2)
                        Count = Count + 1 ;
                        Eps(ComNum, K) = abs(params(Count)) ;
                    end
                end            
                
                
                if Free0orFix1orAbs2_Q(ComNum, K) >= 3 && Global_Check_Q(ComNum, K) == 0
                    Count = Count + 1 ;
                    Qua(ComNum, K) = abs(params(Count)) ;
                    Var = Global_Check_Q(:,:) + (Free0orFix1orAbs2_Q==Free0orFix1orAbs2_Q(ComNum, K))...
                        *Free0orFix1orAbs2_Q(ComNum, K) ;
                    Global_Check_Q = Var ;
                    Global_Para_Q(Free0orFix1orAbs2_Q(ComNum, K)) = Qua(ComNum, K) ;
                elseif Free0orFix1orAbs2_Q(ComNum, K) >= 3 && Global_Check_Q(ComNum, K) ~= 0
                    Qua(ComNum, K) = Global_Para_Q(Free0orFix1orAbs2_Q(ComNum, K)) ;
                else
                    if (Free0orFix1orAbs2_Q(ComNum, K) == 1)
                        Qua(ComNum, K) = FitQ(ComNum, K) ;
                    elseif (Free0orFix1orAbs2_Q(ComNum, K) == 0)
                        Count = Count + 1 ;
                        Qua(ComNum, K) = params(Count) ;
                    elseif (Free0orFix1orAbs2_Q(ComNum, K) == 2)
                        Count = Count + 1 ;
                        Qua(ComNum, K) = abs(params(Count)) ;
                    end
                end
                
                 
                I = 0 ;
                while I < Kmax(ComNum)
                    I = I + 1 ;
                    
                    if Free0orFix1orAbs2_RateM(I, K, ComNum) >= 3 && Global_Check_RateM(I, K, ComNum) == 0
                        Count = Count + 1 ;
                        RateM(I, K, ComNum) = abs(params(Count)) ;
                        Var = Global_Check_RateM(:,:,:) + (Free0orFix1orAbs2_RateM==Free0orFix1orAbs2_RateM(I, K, ComNum))...
                                                                        *Free0orFix1orAbs2_RateM(I, K, ComNum) ;
                        Global_Check_RateM = Var ;                                           
                        Global_Para_RateM(Free0orFix1orAbs2_RateM(I, K, ComNum)) = RateM(I, K, ComNum) ;
                    elseif Free0orFix1orAbs2_RateM(I, K, ComNum) >= 3 && Global_Check_RateM(I, K, ComNum) ~= 0
                        RateM(I, K, ComNum) = Global_Para_RateM(Free0orFix1orAbs2_RateM(I, K, ComNum)) ;
                    else
                        if (Free0orFix1orAbs2_RateM(I, K, ComNum) == 1)
                            RateM(I, K, ComNum) = FitRateM(I, K, ComNum) ;
                        elseif (Free0orFix1orAbs2_RateM(I, K, ComNum) == 0)
                            Count = Count + 1 ;
                            RateM(I, K, ComNum) = params(Count) ;
                        elseif (Free0orFix1orAbs2_RateM(I, K, ComNum) == 2)
                            Count = Count + 1 ;
                            RateM(I, K, ComNum) = abs(params(Count)) ;
                        end
                    end
                end
            end
            
        end
        
        
        K = 0 ;
        while K < VarNmax
            K = K + 1 ;
            I = 0 ;
            while I < VarNmax
                I = I + 1 ;
                
                if Free0orFix1orAbs2_y0(I, K) >= 3 && Global_Check_y0(I, K) == 0
                    Count = Count + 1 ;
                    y0(I, K) = abs(params(Count)) ;
                    Var = Global_Check_y0(:,:) + (Free0orFix1orAbs2_y0==Free0orFix1orAbs2_y0(I, K))...
                        *Free0orFix1orAbs2_y0(I, K) ;
                    Global_Check_y0 = Var ;
                    Global_Para_y0(Free0orFix1orAbs2_y0(I, K)) = y0(I, K) ;
                elseif Free0orFix1orAbs2_y0(I, K) >= 3 && Global_Check_y0(I, K) ~= 0
                    y0(I, K) = Global_Para_y0(Free0orFix1orAbs2_y0(I, K)) ;
                else
                    if (Free0orFix1orAbs2_y0(I, K) == 1)
                        y0(I, K) = Fity0(I, K) ;
                    elseif (Free0orFix1orAbs2_y0(I, K) == 0)
                        Count = Count + 1 ;
                        y0(I, K) = params(Count) ;
                    elseif (Free0orFix1orAbs2_y0(I, K) == 2)
                        Count = Count + 1 ;
                        y0(I, K) = abs(params(Count)) ;
                    end
                end
            end
        end
        
        Var = length(xdata) ;
        Var_FittedCurvePinf = zeros(Kmaxmax, Kmaxmax, Var, ComNumMax) ;
        Var_FittedCurveRinf = zeros(Kmaxmax, Kmaxmax, ComNumMax) ;
        Var_FittedCurveDelta = zeros(Kmaxmax, Kmaxmax, Var, ComNumMax) ;
        
        ComNum = 0 ;
        while ComNum < ComNumMax
            ComNum = ComNum + 1 ;
            
            VarRateM = zeros(Kmax(ComNum), Kmax(ComNum)) ;
            VarRateM(:,:) = RateM(:,:,ComNum) ;
            
            %%% calcurate Ci (concentration of each state i)
            Var = ones(Kmax(ComNum)) - eye(Kmax(ComNum)) ;
            Var = Var .* VarRateM ;    % avoid diffution/photobleach effects
            [EigenValues, ExpMatrix, EigenVector, InvEigenVector] = TK_RateEq_MakeExpMatrix(Var) ;
            
            Var = ones(Kmax(ComNum), 1) ;
            InitialPop_V = Var ./ Kmax(ComNum) ;
            Var = VarRateM + (VarRateM==0)*max(max(VarRateM)) ;
            Var = min(min(Var)) ;%+ 10^-10 ;
            Var = ExpMatrix .^ (1/Var*(10^10)) ;
            RinfM = EigenVector * Var * InvEigenVector ;
            EquiribriumPop_V = RinfM * InitialPop_V ;
            %%%

            %%% Change into Matrix form
            Var = Eps(ComNum, 1:Kmax(ComNum)) ;
            EpsM = diag(Var) ;
            Var = Qua(ComNum, 1:Kmax(ComNum)) ;
            QuaM = diag(Var) ;
            PinfM = diag(EquiribriumPop_V) ;
            %%%
            
            Var_FittedCurveRinf(:,:, ComNum) = EpsM * QuaM * RinfM ;
            
            %%%
            [EigenValues, ExpMatrix, EigenVector, InvEigenVector] = TK_RateEq_MakeExpMatrix(VarRateM) ;
            
            I = 0 ;
            while I < length(xdata)
                I = I + 1 ;
                if I == 1
                    Var = ExpMatrix .^ xdata(I) ;
                else
                    Var = ExpMatrix .^ xdata(I) ;
                end
                RdeltaM = EigenVector * Var * InvEigenVector ;
                
                Var_FittedCurveDelta(:,:,I, ComNum) = EpsM * QuaM * RdeltaM ;
                Var_FittedCurvePinf(:,:,I, ComNum) = PinfM * QuaM * EpsM ; 
            end
        end
        
        Var = length(xdata) ;
        CorTotal = zeros(Kmaxmax, Kmaxmax, Var) ;
        Var1 = zeros(Kmaxmax, Kmaxmax) ;
        Var2 = zeros(Kmaxmax, Kmaxmax) ;
        I = 0 ;
        while I < length(xdata)
            I = I + 1 ;
            ComNum1 = 0 ;
            while ComNum1 < ComNumMax
                ComNum1 = ComNum1 + 1 ;
                CorDelta = zeros(Kmaxmax, Kmaxmax, Var) ;
                ComNum2 = 0 ;
                while ComNum2 < ComNumMax
                    ComNum2 = ComNum2 + 1 ;
                    if ComNum1 ~= ComNum2
                        CorDelta(:,:,I) = CorDelta(:,:,I) + Var_FittedCurveRinf(:,:, ComNum2) ;
                    else
                        CorDelta(:,:,I) = CorDelta(:,:,I) + Var_FittedCurveDelta(:,:,I, ComNum1) ;
                    end
                end
                CorTotal(:,:,I) = CorTotal(:,:,I) + CorDelta(:,:,I) * Var_FittedCurvePinf(:,:,I, ComNum1) ;
            end
        end
        CorTotal(:,:,I) = CorTotal(:,:,I) .* Amp ;
        
        Var = length(xdata) ;
        Var_FittedCurve2 = zeros(VarNmax, VarNmax, Var) ;
        I = 0 ;
        while I < Kmaxmax
            I = I + 1 ;
            StateI = StateAssign(I) ;
            K = 0 ;
            while K < Kmaxmax
                K = K + 1 ;
                StateK = StateAssign(K) ;
                Var_FittedCurve2(StateI,StateK,:) = Var_FittedCurve2(StateI,StateK,:) + CorTotal(I,K,:) ;
            end
        end

        Var = length(xdata) ;
        I = 0 ;
        while I < Var
            I = I + 1 ;
            Var_FittedCurve2(:,:,I) = (Var_FittedCurve2(:,:,I) + Var_FittedCurve2(:,:,I).') /2 ;
        end
        
        FittedCurve = Var_FittedCurve2 ;
        I = 0 ;
        while I < VarNmax
            I = I + 1 ;
            K = 0 ;
            while K < VarNmax
                K = K + 1 ;
                FittedCurve(I,K,:) = y0(I, K) + Var_FittedCurve2(I,K,:) ;
            end
        end
     
        ErrorVector = FittedCurve - ydata ;
        
        if FitNormalizedCor_y1n0 == 1
            I = 0 ;
            while I < max(StateAssign)
                I = I + 1 ;
                K = 0 ;
                while K < max(StateAssign)
                    K = K + 1 ;
                    Var = max(ydata(I,K,:)) ;
                    ErrorVector(I,K,:) = ErrorVector(I,K,:) ./ Var ;
                end
            end
        elseif FitNormalizedCor_y1n0 ~= 0
            I = 0 ;
            while I < max(StateAssign)
                I = I + 1 ;
                K = 0 ;
                while K < max(StateAssign)
                    K = K + 1 ;
                    Var = max(ydata(I,K,:)) ;
                    ErrorVector(I,K,:) = ErrorVector(I,K,:) ./ Var ;
                end
            end
            I = floor(FitNormalizedCor_y1n0/10) ;
            K = FitNormalizedCor_y1n0 - I*10 ;
            ErrorVector(I,K,:) = ErrorVector(I,K,:) .* FitWeightFactor ;
%             Var = ErrorVector(I,K,:) ;
%             ErrorVector = ErrorVector .* 0 ;
%             ErrorVector(I,K,:) = Var ;
        end
        
        Var = size(ydata) ;
        sse = sum(sum(sum(ErrorVector .^ 2)))/(Var(1)*Var(2)*Var(3)) ; 
        
    end

Ana_estimates_A = FitA ;
Ana_estimates_y0 = Fity0 ;
Ana_estimates_y0(:,:) = 0 ;
Ana_estimates_E = FitE ;
Ana_estimates_Q = FitQ ;
Ana_estimates_RateM = FitRateM ;

Count = 0 ;

if (Free0orFix1orAbs2_A == 1)
    Ana_estimates_A = FitA ;
elseif (Free0orFix1orAbs2_A == 0)
    Count = Count + 1 ;
    Ana_estimates_A = estimates(Count) ;
elseif (Free0orFix1orAbs2_A == 2)
    Count = Count + 1 ;
    Ana_estimates_A = abs(estimates(Count)) ;
end

Global_Check_E(:,:) = 0 ;
Global_Check_Q(:,:) = 0 ;
Global_Check_RateM(:,:,:) = 0 ;
Global_Check_y0(:,:) = 0 ;
ComNum = 0 ;
while ComNum < ComNumMax
    ComNum = ComNum + 1 ;

    K = 0 ;
    while K < Kmax(ComNum)
        K = K + 1 ;       
        
        if Free0orFix1orAbs2_E(ComNum, K) >= 3 && Global_Check_E(ComNum, K) == 0
            Count = Count + 1 ;
            Ana_estimates_E(ComNum, K) = abs(estimates(Count)) ;
            Var = Global_Check_E(:,:) + (Free0orFix1orAbs2_E==Free0orFix1orAbs2_E(ComNum, K))...
                *Free0orFix1orAbs2_E(ComNum, K) ;
            Global_Check_E = Var ;
            Global_Para_E(Free0orFix1orAbs2_E(ComNum, K)) = Ana_estimates_E(ComNum, K) ;
        elseif Free0orFix1orAbs2_E(ComNum, K) >= 3 && Global_Check_E(ComNum, K) ~= 0
            Ana_estimates_E(ComNum, K) = Global_Para_E(Free0orFix1orAbs2_E(ComNum, K)) ;
        else
            if (Free0orFix1orAbs2_E(ComNum, K) == 1)
                Ana_estimates_E(ComNum, K) = FitE(ComNum, K) ;
            elseif (Free0orFix1orAbs2_E(ComNum, K) == 0)
                Count = Count + 1 ;
                Ana_estimates_E(ComNum, K) = estimates(Count) ;
            elseif (Free0orFix1orAbs2_E(ComNum, K) == 2)
                Count = Count + 1 ;
                Ana_estimates_E(ComNum, K) = abs(estimates(Count)) ;
            end
        end
        
        
        if Free0orFix1orAbs2_Q(ComNum, K) >= 3 && Global_Check_Q(ComNum, K) == 0
            Count = Count + 1 ;
            Ana_estimates_Q(ComNum, K) = abs(estimates(Count)) ;
            Var = Global_Check_Q(:,:) + (Free0orFix1orAbs2_Q==Free0orFix1orAbs2_Q(ComNum, K))...
                *Free0orFix1orAbs2_Q(ComNum, K) ;
            Global_Check_Q = Var ;
            Global_Para_Q(Free0orFix1orAbs2_Q(ComNum, K)) = Ana_estimates_Q(ComNum, K) ;
        elseif Free0orFix1orAbs2_Q(ComNum, K) >= 3 && Global_Check_Q(ComNum, K) ~= 0
            Ana_estimates_Q(ComNum, K) = Global_Para_Q(Free0orFix1orAbs2_Q(ComNum, K)) ;
        else
            if (Free0orFix1orAbs2_Q(ComNum, K) == 1)
                Ana_estimates_Q(ComNum, K) = FitQ(ComNum, K) ;
            elseif (Free0orFix1orAbs2_Q(ComNum, K) == 0)
                Count = Count + 1 ;
                Ana_estimates_Q(ComNum, K) = estimates(Count) ;
            elseif (Free0orFix1orAbs2_Q(ComNum, K) == 2)
                Count = Count + 1 ;
                Ana_estimates_Q(ComNum, K) = abs(estimates(Count)) ;
            end
        end
        
        
        I = 0 ;
        while I < Kmax(ComNum)
            I = I + 1 ;
            
            if Free0orFix1orAbs2_RateM(I, K, ComNum) >= 3 && Global_Check_RateM(I, K, ComNum) == 0
                Count = Count + 1 ;
                Ana_estimates_RateM(I, K, ComNum) = abs(estimates(Count)) ;
                Var = Global_Check_RateM(:,:,:) + (Free0orFix1orAbs2_RateM==Free0orFix1orAbs2_RateM(I, K, ComNum))...
                                                                *Free0orFix1orAbs2_RateM(I, K, ComNum) ;
                Global_Check_RateM = Var ;
                Global_Para_RateM(Free0orFix1orAbs2_RateM(I, K, ComNum)) = Ana_estimates_RateM(I, K, ComNum) ;
            elseif Free0orFix1orAbs2_RateM(I, K, ComNum) >= 3 && Global_Check_RateM(I, K, ComNum) ~= 0
                Ana_estimates_RateM(I, K, ComNum) = Global_Para_RateM(Free0orFix1orAbs2_RateM(I, K, ComNum)) ;
            else
                if (Free0orFix1orAbs2_RateM(I, K, ComNum) == 1)
                    Ana_estimates_RateM(I, K, ComNum) = FitRateM(I, K, ComNum) ;
                elseif (Free0orFix1orAbs2_RateM(I, K, ComNum) == 0)
                    Count = Count + 1 ;
                    Ana_estimates_RateM(I, K, ComNum) = estimates(Count) ;
                elseif (Free0orFix1orAbs2_RateM(I, K, ComNum) == 2)
                    Count = Count + 1 ;
                    Ana_estimates_RateM(I, K, ComNum) = abs(estimates(Count)) ;
                end
            end
        end
    end
end

K = 0 ;
while K < VarN
    K = K + 1 ;
    I = 0 ;
    while I < VarN
        I = I + 1 ;

        if Free0orFix1orAbs2_y0(I, K) >= 3 && Global_Check_y0(I, K) == 0
            Count = Count + 1 ;
            Ana_estimates_y0(I, K) = abs(estimates(Count)) ;
            Var = Global_Check_y0(:,:) + (Free0orFix1orAbs2_y0==Free0orFix1orAbs2_y0(I, K))...
                *Free0orFix1orAbs2_y0(I, K) ;
            Global_Check_y0 = Var ; 
            Global_Para_y0(Free0orFix1orAbs2_y0(I, K)) = Ana_estimates_y0(I, K) ;
        elseif Free0orFix1orAbs2_y0(I, K) >= 3 && Global_Check_y0(I, K) ~= 0
            Ana_estimates_y0(I, K) = Global_Para_y0(Free0orFix1orAbs2_y0(I, K)) ;
        else
            if (Free0orFix1orAbs2_y0(I, K) == 1)
                Ana_estimates_y0(I, K) = Fity0(I, K) ;
            elseif (Free0orFix1orAbs2_y0(I, K) == 0)
                Count = Count + 1 ;
                Ana_estimates_y0(I, K) = estimates(Count) ;
            elseif (Free0orFix1orAbs2_y0(I, K) == 2)
                Count = Count + 1 ;
                Ana_estimates_y0(I, K) = abs(estimates(Count)) ;
            end
        end
    end
end

Step = (max(xdata)-min(xdata)) / (FittedPoints-1) ;
Ana_FittedX = [min(xdata) : Step : max(xdata)]' ;

Var = length(Ana_FittedX) ;
Var_FittedCurveRinf = zeros(Kmaxmax, Kmaxmax, ComNumMax) ;
Var_FittedCurvePinf = zeros(Kmaxmax, Kmaxmax, Var, ComNumMax) ;
Var_FittedCurveDelta = zeros(Kmaxmax, Kmaxmax, Var, ComNumMax) ;
        
ComNum = 0 ;
while ComNum < ComNumMax
    ComNum = ComNum + 1 ;
    
    VarRateM = zeros(Kmax(ComNum), Kmax(ComNum)) ;
    VarRateM(:,:) = Ana_estimates_RateM(:,:,ComNum) ;
    
    %%% calcurate Ci (concentration of each state i)
    Var = ones(Kmax(ComNum)) - eye(Kmax(ComNum)) ;
    Var = Var .* VarRateM ;    % avoid diffution/photobleach effects
    [EigenValues, ExpMatrix, EigenVector, InvEigenVector] = TK_RateEq_MakeExpMatrix(Var) ;
    
    Var = ones(Kmax(ComNum), 1) ;
    InitialPop_V = Var ./ Kmax(ComNum) ;
    Var = VarRateM + (VarRateM==0)*max(max(VarRateM)) ;
    Var = min(min(Var)) ;%+ 10^-10 ;
    Var = ExpMatrix .^ (1/Var*(10^10)) ;
    RinfM = EigenVector * Var * InvEigenVector ;
    EquiribriumPop_V = RinfM * InitialPop_V ;
    %%%

    %%% Change into Matrix form
    Var = Ana_estimates_E(ComNum, 1:Kmax(ComNum)) ;
    EpsM = diag(Var) ;
    Var = Ana_estimates_Q(ComNum, 1:Kmax(ComNum)) ;
    QuaM = diag(Var) ;
    PinfM = diag(EquiribriumPop_V) ;
    %%%
    
    Var_FittedCurveRinf(:,:, ComNum) = EpsM * QuaM * RinfM ;
    
    %%%
    [EigenValues, ExpMatrix, EigenVector, InvEigenVector] = TK_RateEq_MakeExpMatrix(VarRateM) ;
    
    I = 0 ;
    while I < length(Ana_FittedX)
        I = I + 1 ;
        if I == 1
            Var = ExpMatrix .^ Ana_FittedX(I) ;
        else
            Var = ExpMatrix .^ Ana_FittedX(I) ;
        end
        RdeltaM = EigenVector * Var * InvEigenVector ;
        
        Var_FittedCurveDelta(:,:,I, ComNum) = EpsM * QuaM * RdeltaM ;
        Var_FittedCurvePinf(:,:,I, ComNum) = PinfM * QuaM * EpsM ;
    end
end

Var = length(Ana_FittedX) ;
CorTotal = zeros(Kmaxmax, Kmaxmax, Var) ;
Var1 = zeros(Kmaxmax, Kmaxmax) ;
Var2 = zeros(Kmaxmax, Kmaxmax) ;
I = 0 ;
while I < length(Ana_FittedX)
    I = I + 1 ;
    ComNum1 = 0 ;
    while ComNum1 < ComNumMax
        ComNum1 = ComNum1 + 1 ;
        CorDelta = zeros(Kmaxmax, Kmaxmax, Var) ;
        ComNum2 = 0 ;
        while ComNum2 < ComNumMax
            ComNum2 = ComNum2 + 1 ;
            if ComNum1 ~= ComNum2
                CorDelta(:,:,I) = CorDelta(:,:,I) + Var_FittedCurveRinf(:,:, ComNum2) ;
            else
                CorDelta(:,:,I) = CorDelta(:,:,I) + Var_FittedCurveDelta(:,:,I, ComNum1) ;
            end
        end
        CorTotal(:,:,I) = CorTotal(:,:,I) + CorDelta(:,:,I) * Var_FittedCurvePinf(:,:,I, ComNum1) ;
    end
end
CorTotal(:,:,:) = CorTotal(:,:,:) .* Ana_estimates_A ;

Var = length(Ana_FittedX) ;
Var_FittedCurve2 = zeros(VarNmax, VarNmax, Var) ;
I = 0 ;
while I < Kmaxmax
    I = I + 1 ;
    StateI = StateAssign(I) ;
    K = 0 ;
    while K < Kmaxmax
        K = K + 1 ;
        StateK = StateAssign(K) ;
        Var_FittedCurve2(StateI,StateK,:) = Var_FittedCurve2(StateI,StateK,:) + CorTotal(I,K,:) ;
    end
end

Var = length(Ana_FittedX) ;
I = 0 ;
while I < Var
    I = I + 1 ;
    Var_FittedCurve2(:,:,I) = (Var_FittedCurve2(:,:,I) + Var_FittedCurve2(:,:,I).') /2 ;
end

Ana_FittedY = Var_FittedCurve2 ;
I = 0 ;
while I < VarNmax
    I = I + 1 ;
    K = 0 ;
    while K < VarNmax
        K = K + 1 ;
        Ana_FittedY(I,K,:) = Ana_estimates_y0(I, K) + Var_FittedCurve2(I,K,:) ;
    end
end


Var = length(StateAssign) ;
Ana_estimates = zeros(5+1*Var, Var*ComNumMax) ;
Ana_Free0orFix1orAbs2 = zeros(5+1*Var, Var*ComNumMax) ;

ComNum = 0 ;
while ComNum < ComNumMax
    ComNum = ComNum + 1 ;
    
    Var = length(StateAssign) ;
    Var1 = Var * (ComNum-1) + 1 ;
    Var2 = Var * ComNum ;
    
    Ana_estimates(1,Var1:Var2) = ComNum ;
    Ana_estimates(2,Var1:Var2) = StateAssign(:) ;
    Ana_estimates(3,Var1:Var2) = Ana_estimates_A ;
    Ana_estimates(4,Var1:Var2) = Ana_estimates_E(ComNum,:) ;
    Ana_estimates(5,Var1:Var2) = Ana_estimates_Q(ComNum,:) ;
    Ana_estimates(6:5+Var,Var1:Var2) = Ana_estimates_RateM(:,:,ComNum) ;
    
    Ana_Free0orFix1orAbs2(1,Var1:Var2) = ComNum ;
    Ana_Free0orFix1orAbs2(2,Var1:Var2) = StateAssign(:) ;
    Ana_Free0orFix1orAbs2(3,Var1:Var2) = Free0orFix1orAbs2_A ;
    Ana_Free0orFix1orAbs2(4,Var1:Var2) = Free0orFix1orAbs2_E(ComNum,:) ;
    Ana_Free0orFix1orAbs2(5,Var1:Var2) = Free0orFix1orAbs2_Q(ComNum,:) ;
    Ana_Free0orFix1orAbs2(6:5+Var,Var1:Var2) = Free0orFix1orAbs2_RateM(:,:,ComNum) ;
end

end


