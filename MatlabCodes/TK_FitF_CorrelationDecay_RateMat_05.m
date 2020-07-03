% TK_FitF_CorrelationDecay_RateMat_05

function [Ana_estimates, Ana_Free0orFix1orAbs2, Ana_estimates_Free0orFix1orAbs2_y0, Ana_estimates_A, Ana_estimates_y0, Ana_estimates_E, Ana_estimates_Q, Ana_estimates_RateM, Ana_FittedX, Ana_FittedY, Ana_FittedY_Ratio, model] = ...
    TK_FitF_CorrelationDecay_RateMat_05(...
    xdata, ydata, StateAssign, ...
    FitA, Fity0, FitE, FitQ, FitRateM,...
    Free0orFix1orAbs2_A, Free0orFix1orAbs2_y0, Free0orFix1orAbs2_E, Free0orFix1orAbs2_Q, Free0orFix1orAbs2_RateM,...
    FittedPoints)

Kmax = length(FitE) ;
VarN = max(StateAssign) ;
N1N2max = VarN*(VarN-1)/2 ;
Imax = 1 + N1N2max*3 + Kmax + Kmax + Kmax*Kmax ;
InitialPara = zeros(Imax, 1) ;

Count = 0 ;
if (Free0orFix1orAbs2_A ~= 1)
    Count = Count + 1 ;
    if (Free0orFix1orAbs2_A == 0)
        InitialPara(Count) = FitA ;
    elseif (Free0orFix1orAbs2_A == 2)
        InitialPara(Count) = abs(FitA) ;
    end
end

K = 0 ;
while K < Kmax
    K = K + 1 ;
    
    if (Free0orFix1orAbs2_E(K) ~= 1)
        Count = Count + 1 ;
        if (Free0orFix1orAbs2_E(K) == 0)
            InitialPara(Count) = FitE(K) ;
        elseif (Free0orFix1orAbs2_E(K) == 2)
            InitialPara(Count) = abs(FitE(K)) ;
        end
    end
    if (Free0orFix1orAbs2_Q(K) ~= 1)
        Count = Count + 1 ;
        if (Free0orFix1orAbs2_Q(K) == 0)
            InitialPara(Count) = FitQ(K) ;
        elseif (Free0orFix1orAbs2_Q(K) == 2)
            InitialPara(Count) = abs(FitQ(K)) ;
        end
    end
    
    I = 0 ;
    while I < Kmax
        I = I + 1 ;
        if (Free0orFix1orAbs2_RateM(I, K) ~= 1)
            Count = Count + 1 ;
            if (Free0orFix1orAbs2_RateM(I, K) == 0)
                InitialPara(Count) = FitRateM(I, K) ;
            elseif (Free0orFix1orAbs2_RateM(I, K) == 2)
                InitialPara(Count) = abs(FitRateM(I, K)) ;
            end
        end
    end
end

K = 0 ;
while K < VarN
    K = K + 1 ;
    I = K ;
    while I < VarN
        I = I + 1 ;
        if (Free0orFix1orAbs2_y0(I, K, 1) ~= 1)
            Count = Count + 1 ;
            if (Free0orFix1orAbs2_y0(I, K, 1) == 0)
                InitialPara(Count) = Fity0(I, K, 1) ;
            elseif (Free0orFix1orAbs2_y0(I, K, 1) == 2)
                InitialPara(Count) = abs(Fity0(I, K, 1)) ;
            end
        end
        if (Free0orFix1orAbs2_y0(I, K, 2) ~= 1)
            Count = Count + 1 ;
            if (Free0orFix1orAbs2_y0(I, K, 2) == 0)
                InitialPara(Count) = Fity0(I, K, 2) ;
            elseif (Free0orFix1orAbs2_y0(I, K, 2) == 2)
                InitialPara(Count) = abs(Fity0(I, K, 2)) ;
            end
        end
        if (Free0orFix1orAbs2_y0(I, K, 3) ~= 1)
            Count = Count + 1 ;
            if (Free0orFix1orAbs2_y0(I, K, 3) == 0)
                InitialPara(Count) = Fity0(I, K, 3) ;
            elseif (Free0orFix1orAbs2_y0(I, K, 3) == 2)
                InitialPara(Count) = abs(Fity0(I, K, 3)) ;
            end
        end
    end
end
InitialPara = InitialPara(1:Count) ;

ErrorScaleFactor = zeros(3,N1N2max) ;
Var = size(ydata) ;
Var_ydata_G1 = zeros(Var(1),Var(2),Var(3)) ;
Var_ydata_G2 = zeros(Var(1),Var(2),Var(3)) ;
Var_ydata_G3 = zeros(Var(1),Var(2),Var(3)) ;
Count = 0 ;
I = 0 ;
while I < VarN
    I = I + 1 ;
    K = I ;
    while K < VarN
        K = K + 1 ;
        Var_ydata_G1(I,K,:) = ydata(K,I,:) ./ ydata(K,K,:) ;
        Var_ydata_G2(I,K,:) = (ydata(I,K,:) + ydata(K,I,:)) ./ (ydata(I,I,:) + ydata(K,K,:)) ;
        Var_ydata_G3(I,K,:) = (ydata(I,K,:) .* ydata(K,I,:)) ./ (ydata(I,I,:) .* ydata(K,K,:)) ;
        Var_ydata_G1_mean = mean(Var_ydata_G1(I,K,:)) ;
        Var_ydata_G2_mean = mean(Var_ydata_G2(I,K,:)) ;
        Var_ydata_G3_mean = mean(Var_ydata_G3(I,K,:)) ;
        Var_mean_II = mean(ydata(I,I,:)) ;
        Var_mean_KK = mean(ydata(K,K,:)) ;
        Var_mean_AveIIKK = (Var_mean_II+Var_mean_KK)/2 ;
        VarVar1 = Var_mean_KK*Var_mean_KK/(Var_mean_AveIIKK^2) ;
        VarVar2 = 1 ;
        VarVar3 = Var_mean_II*Var_mean_KK/(Var_mean_AveIIKK^2) ;
        VarVar = max([VarVar1, VarVar2, VarVar3])  ;
        Count = Count + 1 ;
        ErrorScaleFactor(1,Count) = (1)./Var_ydata_G1_mean ;%.* VarVar1 ./ VarVar ;
        ErrorScaleFactor(2,Count) = (1)./Var_ydata_G2_mean ;%.* VarVar2 ./ VarVar ;
        ErrorScaleFactor(3,Count) = (1)./Var_ydata_G3_mean ;%.* VarVar3 ./ VarVar ;
    end
end

model = @expfun ;
options = optimset('MaxFunEvals',10^20,'MaxIter',10^20);
estimates = fminsearch(model, InitialPara, options) ;

    function [sse, FittedCurve] = expfun(params)
        Amp = zeros(1) ;
        y0 = zeros(VarN, VarN, 3) ;
        Eps = zeros(Kmax, 1) ;
        Qua = zeros(Kmax, 1) ;
        RateM = zeros(Kmax, Kmax) ;
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
        
        K = 0 ;
        while K < Kmax
            K = K + 1 ;
            
            if (Free0orFix1orAbs2_E(K) == 1)
                Eps(K) = FitE(K) ;
            elseif (Free0orFix1orAbs2_E(K) == 0)
                Count = Count + 1 ;
                Eps(K) = params(Count) ;
            elseif (Free0orFix1orAbs2_E(K) == 2)
                Count = Count + 1 ;
                Eps(K) = abs(params(Count)) ;
            end
            
            if (Free0orFix1orAbs2_Q(K) == 1)
                Qua(K) = FitQ(K) ;
            elseif (Free0orFix1orAbs2_Q(K) == 0)
                Count = Count + 1 ;
                Qua(K) = params(Count) ;
            elseif (Free0orFix1orAbs2_Q(K) == 2)
                Count = Count + 1 ;
                Qua(K) = abs(params(Count)) ;
            end
            
            I = 0 ;
            while I < Kmax
                I = I + 1 ;
                if (Free0orFix1orAbs2_RateM(I, K) == 1)
                    RateM(I, K) = FitRateM(I, K) ;
                elseif (Free0orFix1orAbs2_RateM(I, K) == 0)
                    Count = Count + 1 ;
                    RateM(I, K) = params(Count) ;
                elseif (Free0orFix1orAbs2_RateM(I, K) == 2)
                    Count = Count + 1 ;
                    RateM(I, K) = abs(params(Count)) ;
                end
            end
        end
        K = 0 ;
        while K < VarN
            K = K + 1 ;
            I = K ;
            while I < VarN
                I = I + 1 ;
                if (Free0orFix1orAbs2_y0(I, K, 1) == 1)
                    y0(I, K, 1) = Fity0(I, K, 1) ;
                elseif (Free0orFix1orAbs2_y0(I, K, 1) == 0)
                    Count = Count + 1 ;
                    y0(I, K, 1) = params(Count) ;
                elseif (Free0orFix1orAbs2_y0(I, K, 1) == 2)
                    Count = Count + 1 ;
                    y0(I, K, 1) = abs(params(Count)) ;
                end
                if (Free0orFix1orAbs2_y0(I, K, 2) == 1)
                    y0(I, K, 2) = Fity0(I, K, 2) ;
                elseif (Free0orFix1orAbs2_y0(I, K, 2) == 0)
                    Count = Count + 1 ;
                    y0(I, K, 2) = params(Count) ;
                elseif (Free0orFix1orAbs2_y0(I, K, 2) == 2)
                    Count = Count + 1 ;
                    y0(I, K, 2) = abs(params(Count)) ;
                end
                if (Free0orFix1orAbs2_y0(I, K, 3) == 1)
                    y0(I, K, 3) = Fity0(I, K, 3) ;
                elseif (Free0orFix1orAbs2_y0(I, K, 3) == 0)
                    Count = Count + 1 ;
                    y0(I, K, 3) = params(Count) ;
                elseif (Free0orFix1orAbs2_y0(I, K, 3) == 2)
                    Count = Count + 1 ;
                    y0(I, K, 3) = abs(params(Count)) ;
                end
            end
        end
       % [Fity0(1, 2, 1), Fity0(1, 2, 2), Fity0(1, 2, 3) ; y0(1, 2, 1), y0(1, 2, 2), y0(1, 2, 3)]
        [EigenValues, ExpMatrix, EigenVector, InvEigenVector] = TK_RateEq_MakeExpMatrix(RateM) ;
        
       % Var = ones(Kmax, 1) ;
       % InitialPop_V = Var ./ Kmax ;
        
        Var = length(xdata) ;
        Var_FittedCurve1 = zeros(Kmax, Kmax, Var) ;
        
        I = 0 ;
        while I < length(xdata)
            I = I + 1 ;
            if I == 1
                Var = ExpMatrix .^ xdata(I) ;
            else
                Var = ExpMatrix .^ xdata(I) ;
            end
            Var_FittedCurve1(:,:,I) = EigenVector * Var * InvEigenVector ;
        end
        
%         I = 0 ;
%         while I < Kmax
%             I = I + 1 ;
%             K = 0 ;
%             while K < Kmax
%                 K = K + 1 ;
%                 Var_FittedCurve(I,K,:) = Var_FittedCurve(I,K,:) * Amp * Eps(I)*Eps(K) * Qua(I)*Qua(K) ;
%             end
%         end
        
        Var = length(xdata) ;
        %Var1 = max(StateAssign) ;
        Var_FittedCurve2 = zeros(VarN, VarN, Var) ;
        I = 0 ;
        while I < Kmax
            I = I + 1 ;
            StateI = StateAssign(I) ;
            K = 0 ;
            while K < Kmax
                K = K + 1 ;
                StateK = StateAssign(K) ;
                %Var_FittedCurve2(StateI,StateK,:) = Var_FittedCurve2(StateI,StateK,:) + Var_FittedCurve1(I,K,:) * Amp * Eps(I)*Eps(K) * Qua(I)*Qua(K) ;
                Var_FittedCurve2(StateI,StateK,:) = Var_FittedCurve2(StateI,StateK,:) + Var_FittedCurve1(I,K,:) * 1 * Eps(I)*Eps(K) * Qua(I)*Qua(K) ;
            end
        end
        
        Var = length(xdata) ;
        I = 0 ;
        while I < Var
            I = I + 1 ;
            Var_FittedCurve2(:,:,I) = (Var_FittedCurve2(:,:,I) + Var_FittedCurve2(:,:,I).') /2 ;
        end
        
        Var2 = VarN*(VarN-1)/2 ;
        ErrorVector = zeros(Var, 3, Var2) ;
        FittedCurve = zeros(VarN, VarN, Var, 3) ;
        Count = 0 ;
        I = 0 ;
        while I < VarN
            I = I + 1 ;
            K = I ;
            while K < VarN
                K = K + 1 ;
                Var_G1 = y0(K, I, 1) + Amp .* (Var_FittedCurve2(K,I,:) ./ Var_FittedCurve2(K,K,:)) ;
                Var_G2 = y0(K, I, 2) + Amp .* ((Var_FittedCurve2(I,K,:) + Var_FittedCurve2(K,I,:)) ./ (Var_FittedCurve2(I,I,:) + Var_FittedCurve2(K,K,:))) ;
                Var_G3 = y0(K, I, 3) + Amp .* ((Var_FittedCurve2(I,K,:) .* Var_FittedCurve2(K,I,:)) ./ (Var_FittedCurve2(I,I,:) .* Var_FittedCurve2(K,K,:))) ;
              %  [y0(K, I, 1), y0(K, I, 2), y0(K, I, 3)]
                FittedCurve(I, K, :, 1) = Var_G1 ;
                FittedCurve(I, K, :, 2) = Var_G2 ;
                FittedCurve(I, K, :, 3) = Var_G3 ;
%                 Var_ydata_G1 = ydata(K,I,:) ./ ydata(K,K,:) ;
%                 Var_ydata_G2 = (ydata(I,K,:) + ydata(K,I,:)) ./ (ydata(I,I,:) + ydata(K,K,:)) ;
%                 Var_ydata_G3 = (ydata(I,K,:) .* ydata(K,I,:)) ./ (ydata(I,I,:) .* ydata(K,K,:)) ;
%                 
%                 Var_ydata_G1_mean = mean(Var_ydata_G1) ;
%                 Var_ydata_G2_mean = mean(Var_ydata_G2) ;
%                 Var_ydata_G3_mean = mean(Var_ydata_G3) ;
%                 Var_mean_II = mean(ydata(I,I,:)) ;
%                 Var_mean_KK = mean(ydata(K,K,:)) ;
%                 Var_mean_AveIIKK = (Var_mean_II+Var_mean_KK)/2 ;
%                 VarVar1 = Var_mean_KK*Var_mean_KK/(Var_mean_AveIIKK^2) ;
%                 VarVar2 = 1 ;
%                 VarVar3 = Var_mean_II*Var_mean_KK/(Var_mean_AveIIKK^2) ;
%                 VarVar = max([VarVar1, VarVar2, VarVar3])  ; 
                Count = Count + 1 ;
                ErrorVector(:,1,Count) = (Var_G1-Var_ydata_G1(I,K,:)).* ErrorScaleFactor(1,Count) ;
                ErrorVector(:,2,Count) = (Var_G2-Var_ydata_G2(I,K,:)).* ErrorScaleFactor(2,Count) ;
                ErrorVector(:,3,Count) = (Var_G3-Var_ydata_G3(I,K,:)).* ErrorScaleFactor(3,Count) ;
                
                ErrorVector(:,:,Count) = ErrorVector(:,:,Count) .* mean(ydata(K,I,:)) ;
            end
        end
        sse = sum(sum(sum(ErrorVector .^ 2)))/(Var*3*Var2) ;
        
%         ErrorVector = FittedCurve - ydata ;
%         I = 0 ;
%         while I < Var1
%             I = I + 1 ;
%             K = 0 ;
%             while K < Var1
%                 K = K + 1 ;
%                 ErrorVector(I,K,:) = ErrorVector(I,K,:) ./ (max(ydata(I,K,:))-min(ydata(I,K,:))) ;
%             end
%         end
       
%         Var = size(ydata) ;
%         sse = sum(sum(sum(ErrorVector .^ 2)))/(Var(1)*Var(2)*Var(3)) ; 
        
    end

%Ana_estimates_A = FitA ;
Ana_estimates_y0 = Fity0 ;
Ana_estimates_y0(:,:,:) = 0 ;
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

K = 0 ;
while K < Kmax
    K = K + 1 ;
    
    if (Free0orFix1orAbs2_E(K) == 1)
        Ana_estimates_E(K) = FitE(K) ;
    elseif (Free0orFix1orAbs2_E(K) == 0)
        Count = Count + 1 ;
        Ana_estimates_E(K) = estimates(Count) ;
    elseif (Free0orFix1orAbs2_E(K) == 2)
        Count = Count + 1 ;
        Ana_estimates_E(K) = abs(estimates(Count)) ;
    end
    
    if (Free0orFix1orAbs2_Q(K) == 1)
        Ana_estimates_Q(K) = FitQ(K) ;
    elseif (Free0orFix1orAbs2_Q(K) == 0)
        Count = Count + 1 ;
        Ana_estimates_Q(K) = estimates(Count) ;
    elseif (Free0orFix1orAbs2_Q(K) == 2)
        Count = Count + 1 ;
        Ana_estimates_Q(K) = abs(estimates(Count)) ;
    end
    
    I = 0 ;
    while I < Kmax
        I = I + 1 ;
        if (Free0orFix1orAbs2_RateM(I, K) == 1)
            Ana_estimates_RateM(I, K) = FitRateM(I, K) ;
        elseif (Free0orFix1orAbs2_RateM(I, K) == 0)
            Count = Count + 1 ;
            Ana_estimates_RateM(I, K) = estimates(Count) ;
        elseif (Free0orFix1orAbs2_RateM(I, K) == 2)
            Count = Count + 1 ;
            Ana_estimates_RateM(I, K) = abs(estimates(Count)) ;
        end
    end
end
K = 0 ;
while K < VarN
    K = K + 1 ;
    I = K ;
    while I < VarN
        I = I + 1 ;
        if (Free0orFix1orAbs2_y0(I, K, 1) == 1)
            Ana_estimates_y0(I, K, 1) = Fity0(I, K, 1) ;
        elseif (Free0orFix1orAbs2_y0(I, K, 1) == 0)
            Count = Count + 1 ;
            Ana_estimates_y0(I, K, 1) = estimates(Count) ;
        elseif (Free0orFix1orAbs2_y0(I, K, 1) == 2)
            Count = Count + 1 ;
            Ana_estimates_y0(I, K, 1) = abs(estimates(Count)) ;
        end
        if (Free0orFix1orAbs2_y0(I, K, 2) == 1)
            Ana_estimates_y0(I, K, 2) = Fity0(I, K, 2) ;
        elseif (Free0orFix1orAbs2_y0(I, K, 2) == 0)
            Count = Count + 1 ;
            Ana_estimates_y0(I, K, 2) = estimates(Count) ;
        elseif (Free0orFix1orAbs2_y0(I, K, 2) == 2)
            Count = Count + 1 ;
            Ana_estimates_y0(I, K, 2) = abs(estimates(Count)) ;
        end
        if (Free0orFix1orAbs2_y0(I, K, 3) == 1)
            Ana_estimates_y0(I, K, 3) = Fity0(I, K, 3) ;
        elseif (Free0orFix1orAbs2_y0(I, K, 3) == 0)
            Count = Count + 1 ;
            Ana_estimates_y0(I, K, 3) = estimates(Count) ;
        elseif (Free0orFix1orAbs2_y0(I, K, 3) == 2)
            Count = Count + 1 ;
            Ana_estimates_y0(I, K, 3) = abs(estimates(Count)) ;
        end
    end
end

Step = (max(xdata)-min(xdata)) / (FittedPoints-1) ;
Ana_FittedX = [min(xdata) : Step : max(xdata)]' ;
Var_Ana_FittedY = zeros(Kmax, Kmax, FittedPoints) ;
        
[EigenValues, ExpMatrix, EigenVector, InvEigenVector] = TK_RateEq_MakeExpMatrix(Ana_estimates_RateM) ;

I = 0 ;
while I < length(Ana_FittedX)
    I = I + 1 ;
    if I == 1
        Var = ExpMatrix .^ Ana_FittedX(I) ;
    else
        Var = ExpMatrix .^ Ana_FittedX(I) ;
    end
    Var_Ana_FittedY(:,:,I) = EigenVector * Var * InvEigenVector ;
end

I = 0 ;
while I < Kmax
    I = I + 1 ;
    K = 0 ;
    while K < Kmax
        K = K + 1 ;
        Var_Ana_FittedY(I,K,:) = Var_Ana_FittedY(I,K,:) * Ana_estimates_A ...
                             * Ana_estimates_E(I)*Ana_estimates_E(K) * Ana_estimates_Q(I)*Ana_estimates_Q(K) ;
    end
end

Var = length(Ana_FittedX) ;
%Var1 = max(StateAssign) ;
Ana_FittedY = zeros(VarN, VarN, Var) ;
I = 0 ;
while I < Kmax
    I = I + 1 ;
    StateI = StateAssign(I) ;
    K = 0 ;
    while K < Kmax
        K = K + 1 ;
        StateK = StateAssign(K) ;
        Ana_FittedY(StateI,StateK,:) = Ana_FittedY(StateI,StateK,:) + Var_Ana_FittedY(I,K,:) ;
    end
end

Ana_FittedY_Ratio = zeros(VarN, VarN, Var, 3) ;
I = 0 ;
while I < VarN
    I = I + 1 ;
    K = I ;
    while K < VarN
        K = K + 1 ;
        Ana_FittedY_Ratio(I,K,:,1) = Ana_estimates_y0(K, I, 1) + Ana_estimates_A .* (Ana_FittedY(K,I,:) ./ Ana_FittedY(K,K,:)) ;
        Ana_FittedY_Ratio(I,K,:,2) = Ana_estimates_y0(K, I, 2) + Ana_estimates_A .* ((Ana_FittedY(I,K,:) + Ana_FittedY(K,I,:)) ./ (Ana_FittedY(I,I,:) + Ana_FittedY(K,K,:))) ;
        Ana_FittedY_Ratio(I,K,:,3) = Ana_estimates_y0(K, I, 3) + Ana_estimates_A .* ((Ana_FittedY(I,K,:) .* Ana_FittedY(K,I,:)) ./ (Ana_FittedY(I,I,:) .* Ana_FittedY(K,K,:))) ;
        Ana_FittedY_Ratio(K,I,:,1) =  Ana_FittedY_Ratio(I,K,:,1) ;
        Ana_FittedY_Ratio(K,I,:,2) =  Ana_FittedY_Ratio(I,K,:,2) ;
        Ana_FittedY_Ratio(K,I,:,3) =  Ana_FittedY_Ratio(I,K,:,3) ;
    end
end

Var = length(StateAssign) ;
Ana_estimates = zeros(4+1*Var, Var) ;
Ana_estimates(1,:) = StateAssign(:) ;
Ana_estimates(2,:) = Ana_estimates_A(:) ;
Ana_estimates(3,:) = Ana_estimates_E(:) ;
Ana_estimates(4,:) = Ana_estimates_Q(:) ;
Ana_estimates(5:4+Var,:) = Ana_estimates_RateM(:,:) ;
%Ana_estimates(5+Var:4+2*Var,:) = Ana_estimates_y0(:,:) ;

Ana_Free0orFix1orAbs2 = zeros(4+1*Var, Var) ;
Ana_Free0orFix1orAbs2(1,:) = StateAssign(:) ;
Ana_Free0orFix1orAbs2(2,:) = Free0orFix1orAbs2_A(:) ;
Ana_Free0orFix1orAbs2(3,:) = Free0orFix1orAbs2_E(:) ;
Ana_Free0orFix1orAbs2(4,:) = Free0orFix1orAbs2_Q(:) ;
Ana_Free0orFix1orAbs2(5:4+Var,:) = Free0orFix1orAbs2_RateM(:,:) ;
%Ana_Free0orFix1orAbs2(5+Var:4+2*Var,:) = Free0orFix1orAbs2_y0(:,:) ;

Ana_estimates_Free0orFix1orAbs2_y0 = zeros(2+3+3, N1N2max) ;
Count = 0 ;
K = 0 ;
while K < VarN
    K = K + 1 ;
    I = K ;
    while I < VarN
        I = I + 1 ;
        Count = Count + 1 ;
        Ana_estimates_Free0orFix1orAbs2_y0(1, Count) = I ;
        Ana_estimates_Free0orFix1orAbs2_y0(2, Count) = K ;
        Ana_estimates_Free0orFix1orAbs2_y0(3:5, Count) = Ana_estimates_y0(I, K, 1:3) ;
        Ana_estimates_Free0orFix1orAbs2_y0(6:8, Count) = Free0orFix1orAbs2_y0(I, K, 1:3) ;
    end
end

end


