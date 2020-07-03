% TK_MyMain_Analyze

%% Input parameters

Ana_TauRange = 0.4 ;%0.62 ;
PeakPosEachState = [0.05, 0.6, 1.85, 3.35] ; % LHCSR_V_7
PeakPosEachState = [0.45, 1.6, 2.8] ; %0.75

%%%% Ana_TargetN1N2 = [1, 2] ;
FuncType = 2 ;  % 1: gij/gii  2: (gij+gji)/(gii+gjj)

FitPara = [0, 0.03, 10, 1]' ;
Free0orFix1orAbs2 = [0, 2, 2, 1]' ;
% FitPara = [0, 0.02, 1000, 1 ; ...
%            0, 0.02, 1, 1]' ;
% Free0orFix1orAbs2 = [1, 1, 1, 1 ; ...
%                      0, 1, 2, 1]' ;
FittedPoints = 100000 ;

%%
TauRange = zeros(NumOfState, 2) ;
N = 0 ;
while N < NumOfState
    N = N + 1 ;
%    [VVV, Var] = max(Gresult_Mat_A_Global(:,N)) ;
%    Var = Var*Tau_step ;
    Var = PeakPosEachState(N) ;

    Var1 = Var - Ana_TauRange/2 ;
    Tau_I_Min = floor(Var1/Tau_step) ;
    if (Tau_I_Min < 1)
    	Tau_I_Min = 1 ;
    elseif (Tau_I_Min > length(Tau))
        Tau_I_Min = length(Tau) ;
    end
    
    Var1 = Var + Ana_TauRange/2 ;
    Tau_I_Max = ceil(Var1/Tau_step) ;
    if (Tau_I_Max < 1)
    	Tau_I_Max = 1 ;
    elseif (Tau_I_Max > length(Tau))
        Tau_I_Max = length(Tau) ;
    end
    
    TauRange(N, 1) = Tau_I_Min ;
    TauRange(N, 2) = Tau_I_Max ;
end

Tmax = length(dT) ;
Correlation = zeros(NumOfState, NumOfState, Tmax) ;
CorrelationRatio = zeros(NumOfState, NumOfState, Tmax) ;
CorrelationRatio_N1N2 = zeros(Tmax, 1) ;
T = 0 ;
while T < Tmax
    T = T + 1 ;
    TargetMat = Gresult_Mat_M_2DFLC_dTfun(:,:,T) ;
    
    %% estimate Correlation at each dT
    N1 = 0 ;
    while N1 < NumOfState
        N1 = N1 + 1 ;
        Extract_Range_I = [TauRange(N1, 1) TauRange(N1, 2)] ;
        
        N2 = 0 ;
        while N2 < NumOfState
            N2 = N2 + 1 ;
            Extract_Range_K = [TauRange(N2, 1) TauRange(N2, 2)] ;
            
            Sum_ExtractMat = 0 ;
            I = Extract_Range_I(1) - 1 ;
            while I < Extract_Range_I(2)
                I = I + 1 ;
                K = Extract_Range_K(1) - 1 ;
                while K < Extract_Range_K(2)
                    K = K + 1 ;
                    Sum_ExtractMat = Sum_ExtractMat + TargetMat(I, K) ;
                end
            end
            Correlation(N1, N2, T) = Sum_ExtractMat ;
        end
    end
    
    %% estimate CorrelationRatio at each dT
    N1 = 0 ;
    while N1 < NumOfState
        N1 = N1 + 1 ;      
        N2 = 0 ;
        while N2 < NumOfState
            N2 = N2 + 1 ;
            if FuncType == 1
                CorrelationRatio(N1, N2, T) = Correlation(N1, N2, T) / Correlation(N1, N1, T) ;
            elseif FuncType == 2
                CorrelationRatio(N1, N2, T) = (Correlation(N1, N2, T)+Correlation(N2, N1, T)) / (Correlation(N1, N1, T)+Correlation(N2, N2, T)) ;
            end
        end
    end
    
    CorrelationRatio_N1N2(T) = CorrelationRatio(Ana_TargetN1N2(1), Ana_TargetN1N2(2), T) ;
    %%
end

Ana_TauRange = TauRange ;
Ana_Correlation = Correlation ;
Ana_CorrelationRatio = CorrelationRatio ;
Ana_CorrelationRatio_N1N2 = CorrelationRatio_N1N2 ;

%% Fitting
[estimates, Ana_FittedX, Ana_FittedY, model] = TK_FitF_CorrelationDecay_Multi(dT, Ana_CorrelationRatio_N1N2, FitPara, Free0orFix1orAbs2, FuncType, FittedPoints) ;
% Ana_estimates(1) = y0
% Ana_estimates(2) = A
% Ana_estimates(3) = K
% Ana_estimates(4) = RateN
%%

%clear Var Var1 Tau_I_Min Tau_I_Max N1 N2 T Tmax K I TargetMat Extract_Range_I
%clear Extract_Range_K Sum_ExtractMat TauRange Correlation CorrelationRatio CorrelationRatio_N1N2



