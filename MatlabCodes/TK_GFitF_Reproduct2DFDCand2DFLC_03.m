% TK_GFitF_Reproduct2DFDCand2DFLC_03
% Reproduce 2DFLC and 2DFDC map from estimates.

function [Mat_M_2DFLC_dTfun, Mat_M_Model_dTfun, Mat_A, Mat_G_dTfun, y0_dTfun, EstQ, Kai2, EntropyS] = TK_GFitF_Reproduct2DFDCand2DFLC_03(...
    estimates, Initial_Mat_A_Global, Initial_Mat_G_dTfun, Initial_y0_dTfun, Fix1orNot0orAbs2orAmp3_Mat_A, Fix1orNot0orAbs2_Mat_G_dTfun, Fix1orNot0orAbs2_y0_dTfun, RegulatorConst,...
    Mat_2DFDC_It, Mat_2DFDC_dTfun, Mat_2DFDC_cor_dTfun,...
    NumOfState, Tau, ExpCurve, mi)

NumOfComponent = length(Tau) ;
NumOfdT = size(Mat_2DFDC_dTfun) ;
if length(NumOfdT) == 2
    NumOfdT = 1 ;
else
    NumOfdT = NumOfdT(3) ;
end

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
        Var_Amp = abs(estimates(Count)) ;
        Mat_A(:, K) = Initial_Mat_A_Global(:, K) .* Var_Amp ;
    else
        I = 0 ;
        while I < NumOfComponent
            I = I + 1 ;
            Count = Count + 1 ;
            % InitialPara = [Tau_distribution at n=1], [Tau_dis. at n=2], [Tau_dis. at n=3],...
            if Fix1orNot0orAbs2orAmp3_Mat_A(K) == 2
                Mat_A(I, K) = abs(estimates(Count)) ;
            else
                Mat_A(I, K) = estimates(Count) ;
            end
        end
    end
end

% just for avoiding error of log(0) = inf
Var = max(max(Mat_A)) - min(min(Mat_A)) ;
Mat_A = Mat_A + (Mat_A==0).* (Var*0.0000001) ;


%%

Kai2_dTfun = zeros(NumOfdT, 1) ;
Mat_G_dTfun = zeros(NumOfState, NumOfState, NumOfdT) ;
y0_dTfun = zeros(NumOfdT, 1) ;
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
                Mat_G_dTfun(I, K, T) = Initial_Mat_G_dTfun(I, K, T) ;    
            elseif Fix1orNot0orAbs2_Mat_G_dTfun(I, K, T) == 2
                Count = Count + 1 ;
                Mat_G_dTfun(I, K, T) = abs(estimates(Count)) ;
            elseif Fix1orNot0orAbs2_Mat_G_dTfun(I, K, T) == 0
                Count = Count + 1 ;
                Mat_G_dTfun(I, K, T) = estimates(Count) ;
            elseif Fix1orNot0orAbs2_Mat_G_dTfun(I, K, T) == 3
                if I <= K
                    Count = Count + 1 ;
                    Mat_G_dTfun(I, K, T) = abs(estimates(Count)) ;
                    Mat_G_dTfun(K, I, T) = abs(estimates(Count)) ;
                end
            end
        end
    end
    Mat_G = Mat_G_dTfun(:, :, T) ;
    
    if Fix1orNot0orAbs2_y0_dTfun(T) ~= 1
        Count = Count + 1 ;
        % InitialPara = [Tau_distribution at n=1], [Tau_dis. at n=2], [Tau_dis. at n=3],...
        if Fix1orNot0orAbs2_y0_dTfun(T) == 2
            y0_dTfun(T) = abs(estimates(Count)) ;
        else
            y0_dTfun(T) = estimates(Count) ;
        end
    else
        y0_dTfun(T) = Initial_y0_dTfun(T) ;
    end

    %% model 2D map at dT = shortest
    
    Dif_Mat_2DFDC_It = Mat_2DFDC_It ;
    Var = length(Mat_2DFDC_It) ;
    Dif_Mat_2DFDC_It(2:Var) = Mat_2DFDC_It(2:Var) - Mat_2DFDC_It(1:Var-1) ;
    Dif_Mat_2DFDC_It(1) = Dif_Mat_2DFDC_It(2) ;
    Dif_Mat_2DFDC_It = kron(Dif_Mat_2DFDC_It, Dif_Mat_2DFDC_It') ;

    Mat_M_2DFLC = Mat_A * Mat_G * Mat_A' ;
    Mat_M_Model = ExpCurve * Mat_M_2DFLC * ExpCurve' + y0_dTfun(T)*Dif_Mat_2DFDC_It ;
    
    if T == 1
        Var = size(Mat_M_Model) ;
        Mat_M_Model_dTfun = zeros(Var(1), Var(2), NumOfdT) ;
        Var = size(Mat_M_2DFLC) ;
        Mat_M_2DFLC_dTfun = zeros(Var(1), Var(2), NumOfdT) ;
    end
    
    Mat_M_Model_dTfun(:, :, T) = Mat_M_Model(:,:) ;
    Mat_M_2DFLC_dTfun(:, :, T) = Mat_M_2DFLC(:,:) ;
    
    %% Kai^2 estimation
    %MatCorIn = Mat0_2DFDC_cor ; % linear Mat is used
    %MatIn = Mat0_2DFDC ;        % linear Mat is used
    Mat_2DFDC = Mat_2DFDC_dTfun(:, :, T) ;
    Mat_2DFDC_cor = Mat_2DFDC_cor_dTfun(:, :, T) ;

%    Var1 = mean(mean(mean(Mat_2DFDC_dTfun))) ;
    Var1 = mean(mean(Mat_2DFDC)) ;
%    Var1 = max(max(max(Mat_2DFDC_dTfun))) ;
%     Var = Mat_2DFDC + (Mat_2DFDC==0) * max(max(Mat_2DFDC)) ;
%     Var1 = min(min(Var)) / 1 ;
    Var = ((Mat_2DFDC_cor - Mat_M_Model).^2) ./ (Mat_2DFDC+Var1) ;  % Var1 is just for avoiding error of Num/0 = inf
%     Var = Area .* Var ;
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

EstQ = Kai2 - 2*EntropyS/RegulatorConst ;

end