% TK_FitF_Reproduct1DFDC
% Reproduce 1DFDC from estimates.

function [Mat_M_Model, Mat_A, EstQ, Kai2, EntropyS] = TK_FitF_Reproduct1DFDC(...
    estimates, Initial_Mat_A, Fix1orNot0orAbs2, RegulatorConst,...
    Mat_1DFDC_It, Mat_1DFDC, Mat_1DFDC_cor,...
    Tau, ExpCurve, mi)

NumOfState = 1 ;
NumOfComponent = length(Tau) ;

%% Set fitting parameters
Count = 0 ;

if Fix1orNot0orAbs2(1) ~= 1
    Mat_A = zeros(NumOfComponent, NumOfState) ;
    K = 0 ;
    while K < NumOfState
        K = K + 1 ;
        I = 0 ;
        while I < NumOfComponent
            I = I + 1 ;
            Count = Count + 1 ;
            % InitialPara = [Tau_distribution at n=1], [Tau_dis. at n=2], [Tau_dis. at n=3],...
            if Fix1orNot0orAbs2(1) == 2
                Mat_A(I, K) = abs(estimates(Count)) ;
            else
                Mat_A(I, K) = estimates(Count) ;
            end
        end
    end
else
    Mat_A = Initial_Mat_A ;
end

% just for avoiding error
Var = max(max(Mat_A)) - min(min(Mat_A)) ;
Mat_A = Mat_A + (Mat_A==0).* (Var*0.0000001) ;

Imax = length(Mat_1DFDC) ;
Area = Mat_1DFDC ;
Area(1) = (Mat_1DFDC_It(2)-Mat_1DFDC_It(1)) ;
I = 1 ;
while I < Imax
    I = I + 1 ;
    Area(I) = (Mat_1DFDC_It(I)-Mat_1DFDC_It(I-1)) ;
end

%% model 2D map at dT = shortest

Mat_M_Model = (Mat_A' * ExpCurve')' ;

%% Kai^2 estimation
%MatCorIn = Mat0_2DFDC_cor ; % linear Mat is used
%MatIn = Mat0_2DFDC ;        % linear Mat is used

Var1 = mean(Mat_1DFDC) ;
Var = ((Mat_1DFDC_cor - Mat_M_Model).^2) ./ (Mat_1DFDC+Var1) ;  % Var1 is just for avoiding error
Var = Area .* Var ;
Imax = length(Mat_1DFDC) ;
Kai2 = sum(Var) / (Imax) ;

%% entropy S estimation
EntropyS = 0 ;

K = 0 ;
while K < NumOfState
    K = K + 1 ;
    Var = Mat_A(:, K) ;
    Var1 = sum(Var) ;
    
    Var = Mat_A(:, K) ./ (mi(:, K) + max(mi(:, K))*10^-10) ;
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