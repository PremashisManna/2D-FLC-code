function [estimates, model] = TK_GFitF_2DMEM_05(...
    Initial_Mat_A_Global, Initial_Mat_G_dTfun, Initial_y0_dTfun, Fix1orNot0orAbs2orAmp3_Mat_A, Fix1orNot0orAbs2_Mat_G_dTfun, Fix1orNot0orAbs2_y0_dTfun, RegulatorConst, ...
    Mat_2DFDC_It, Mat_2DFDC_dTfun, Mat_2DFDC_cor_dTfun, ...
    NumOfState, Tau, ExpCurve, mi)
%%
NumOfdT = size(Initial_Mat_G_dTfun) ;

if length(NumOfdT) == 2
    NumOfdT = 1 ;
else
    NumOfdT = NumOfdT(3) ;
end
NumOfComponent = length(Tau) ;
Var1 = NumOfComponent*NumOfState ;  % number of Mat_A elements
Var2 = (NumOfState+1)*NumOfState/2 ;    % number of triangular Mat_U elements
Fit_InitialPara = zeros((Var1+Var2+1)*NumOfdT, 1) ;

Count = 0 ;
K = 0 ;
while K < NumOfState
    K = K + 1 ;
    if Fix1orNot0orAbs2orAmp3_Mat_A(K) == 1
    elseif Fix1orNot0orAbs2orAmp3_Mat_A(K) == 3
        Count = Count + 1 ;
        Fit_InitialPara(Count) = 1 ;
    else
        I = 0 ;
        while I < NumOfComponent
            I = I + 1 ;
            Count = Count + 1 ;
            % InitialPara = [Tau_distribution at n=1], [Tau_dis. at n=2], [Tau_dis. at n=3],...
            if Fix1orNot0orAbs2orAmp3_Mat_A(K) == 2
                Fit_InitialPara(Count) = abs(Initial_Mat_A_Global(I, K)) ;
            else
                Fit_InitialPara(Count) = Initial_Mat_A_Global(I, K) ;
            end
        end
    end
end

T = 0 ;
while T < NumOfdT
    T = T + 1 ;
    K = 0 ;
    while K < NumOfState
        K = K + 1 ;
        I = 0 ;
        while I < NumOfState
            I = I + 1 ;
            if Fix1orNot0orAbs2_Mat_G_dTfun(I, K, T) == 1
                % none
            elseif Fix1orNot0orAbs2_Mat_G_dTfun(I, K, T) == 2
                Count = Count + 1 ;
                Fit_InitialPara(Count) = abs(Initial_Mat_G_dTfun(I, K, T)) ;
            elseif Fix1orNot0orAbs2_Mat_G_dTfun(I, K, T) == 0
                Count = Count + 1 ;
                Fit_InitialPara(Count) = Initial_Mat_G_dTfun(I, K, T) ;
            elseif Fix1orNot0orAbs2_Mat_G_dTfun(I, K, T) == 3
                if I <= K
                    Count = Count + 1 ;
                    Fit_InitialPara(Count) = abs(Initial_Mat_G_dTfun(I, K, T)) ;
                end
            end
        end
    end
    
    if Fix1orNot0orAbs2_y0_dTfun(T) ~= 1
        Count = Count + 1 ;
        % InitialPara = [Tau_distribution at n=1], [Tau_dis. at n=2], [Tau_dis. at n=3],...
        if Fix1orNot0orAbs2_y0_dTfun(T) == 2
            Fit_InitialPara(Count) = abs(Initial_y0_dTfun(T)) ;
        else
            Fit_InitialPara(Count) = Initial_y0_dTfun(T) ;
        end
    end
end

Fit_InitialPara = Fit_InitialPara(1:Count) ;

Dif_Mat_2DFDC_It = Mat_2DFDC_It ;
Var = length(Mat_2DFDC_It) ;
Dif_Mat_2DFDC_It(2:Var) = Mat_2DFDC_It(2:Var) - Mat_2DFDC_It(1:Var-1) ;
Dif_Mat_2DFDC_It(1) = Dif_Mat_2DFDC_It(2) ;
Dif_Mat_2DFDC_It = kron(Dif_Mat_2DFDC_It, Dif_Mat_2DFDC_It') ;

%%
% Call fminsearch with a random starting point.
%start_point = rand(1, 2) ;
model = @expfun ;
options = optimset('MaxFunEvals',10^4);
estimates = fminsearch(model, Fit_InitialPara, options) ;
% expfun accepts curve parameters as inputs, and outputs Variance,
% the sum of squares error for A*exp(-lambda*xdata)-ydata,
% and the FittedCurve. FMINSEARCH only needs Variance, but we want
% to plot the FittedCurve at the end.

    function [EstQ, Mat_M_Model_dTfun] = expfun(params)
        %% Set fitting parameters   
        Count = 0 ;
        
        Mat_A = zeros(NumOfComponent, NumOfState) ;
        K = 0 ;
        while K < NumOfState
            K = K + 1 ;
            if Fix1orNot0orAbs2orAmp3_Mat_A(K) == 1
                Mat_A(:, K) = Initial_Mat_A_Global(:, K) ;
            elseif Fix1orNot0orAbs2orAmp3_Mat_A(K) == 3
                Count = Count + 1 ;
                Var_Amp = abs(params(Count)) ;
                Mat_A(:, K) = Initial_Mat_A_Global(:, K) .* Var_Amp ;
            else
                I = 0 ;
                while I < NumOfComponent
                    I = I + 1 ;
                    Count = Count + 1 ;
                    % InitialPara = [Tau_distribution at n=1], [Tau_dis. at n=2], [Tau_dis. at n=3],...
                    if Fix1orNot0orAbs2orAmp3_Mat_A(K) == 2
                        Mat_A(I, K) = abs(params(Count)) ;
                    else
                        Mat_A(I, K) = params(Count) ;
                    end
                end
            end
        end

        % just for avoiding error of log(0) = inf
        Var = max(max(Mat_A)) - min(min(Mat_A)) ;
        Mat_A = Mat_A + (Mat_A==0).* (Var*0.0000001) ;
        
        Kai2_dTfun = zeros(NumOfdT, 1) ;
        T = 0 ;
        while T < NumOfdT
            T = T + 1 ;
             
            Mat_G = zeros(NumOfState, NumOfState) ;
            K = 0 ;
            while K < NumOfState
                K = K + 1 ;
                I = 0 ;
                while I < NumOfState
                    I = I + 1 ;
                    if Fix1orNot0orAbs2_Mat_G_dTfun(I, K, T) == 1
                        Mat_G(I, K) = Initial_Mat_G_dTfun(I, K, T) ;
                    elseif Fix1orNot0orAbs2_Mat_G_dTfun(I, K, T) == 2
                        Count = Count + 1 ;
                        Mat_G(I, K) = abs(params(Count)) ;
                    elseif Fix1orNot0orAbs2_Mat_G_dTfun(I, K, T) == 0
                        Count = Count + 1 ;
                        Mat_G(I, K) = params(Count) ;
                    elseif Fix1orNot0orAbs2_Mat_G_dTfun(I, K, T) == 3
                        if I <= K
                            Count = Count + 1 ;
                            Mat_G(I, K) = abs(params(Count)) ;
                            Mat_G(K, I) = abs(params(Count)) ;
                        end
                    end
                end
            end

            if Fix1orNot0orAbs2_y0_dTfun(T) ~= 1
                Count = Count + 1 ;
                % InitialPara = [Tau_distribution at n=1], [Tau_dis. at n=2], [Tau_dis. at n=3],...
                if Fix1orNot0orAbs2_y0_dTfun(T) == 2
                    y0 = abs(params(Count)) ;
                else
                    y0 = params(Count) ;
                end
            else
                y0 = Initial_y0_dTfun(T) ;
            end
            
            %% model 2D map at dT = shortest
            
            Mat_M_2DFLC = Mat_A * Mat_G * Mat_A' ;
            Mat_M_Model = ExpCurve * Mat_M_2DFLC * ExpCurve' + y0*Dif_Mat_2DFDC_It ;
            if T == 1
                Var = size(Mat_M_Model) ;
                Mat_M_Model_dTfun = zeros(Var(1), Var(2), NumOfdT) ;
            end
            Mat_M_Model_dTfun(:, :, T) = Mat_M_Model(:,:) ;
            
            %% Kai^2 estimation
            %MatCorIn = Mat0_2DFDC_cor ; % linear Mat is used
            %MatIn = Mat0_2DFDC ;        % linear Mat is used
            Mat_2DFDC = Mat_2DFDC_dTfun(:, :, T) ;
            Mat_2DFDC_cor = Mat_2DFDC_cor_dTfun(:, :, T) ;
            
%            Var1 = mean(mean(mean(Mat_2DFDC_dTfun))) ;
            Var1 = mean(mean(Mat_2DFDC)) ;
%            Var1 = max(max(max(Mat_2DFDC_dTfun))) ;
%             Var = Mat_2DFDC + (Mat_2DFDC==0) * max(max(Mat_2DFDC)) ;
%             Var1 = min(min(Var)) / 1 ;
            Var = ((Mat_2DFDC_cor - Mat_M_Model).^2) ./ (Mat_2DFDC+Var1) ;  % Var1 is just for avoiding error of Num/0 = inf
%             Var = Area .* Var ;
            Imax = size(Mat_2DFDC) ;
            Imax = Imax(1,1) ;
            Kai2_dTfun(T) = sum(sum(Var)) / (Imax*Imax) ;           
        end
        Kai2 = sum(Kai2_dTfun) / NumOfdT ;
        
        %% entropy S estimation
        EntropyS = 0 ;
        
        K = 0 ;
        while K < NumOfState
            K = K + 1 ;
            Var = Mat_A(:, K) ;
            Var1 = sum(Var) ;
            Var = (Mat_A(:, K) + max(Var)*10^-10) ./ (mi(:, K) + max(mi(:, K))*10^-10) ;
            if sum(isinf(Var)) > 0
                display('inf in Mat_A / mi')
            end
            Var = Mat_A(:, K) .* log(Var) ;
            Var2 = sum(Var) ;
            Var = mi(:, K) ;
            Var3 = sum(Var) ;
            EntropyS = EntropyS + (Var1 - Var3 - Var2) ;
        end
        
        %% Estimator Q
        %RegulatorConst = 0.01 ;
        EstQ = Kai2 - 2*EntropyS/RegulatorConst ;
         
    end

end

