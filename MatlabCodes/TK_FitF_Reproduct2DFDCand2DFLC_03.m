% TK_FitF_Reproduct2DFDCand2DFLC_03
% Reproduce 2DFLC and 2DFDC map from estimates.

function [Mat_M_2DFLC, Mat_M_Model, Mat_A, Mat_G, y0, EstQ, Kai2, EntropyS] = TK_FitF_Reproduct2DFDCand2DFLC_03(...
    estimates, Initial_Mat_A, Initial_Mat_G, Initial_y0, Fix1orNot0orAbs2orAmp3_Mat_A, Fix1orNot0orAbs2_Mat_G, Fix1orNot0orAbs2_y0, RegulatorConst,...
    Mat_2DFDC_It, Mat_2DFDC, Mat_2DFDC_cor,...
    NumOfState, Tau, ExpCurve, mi)

NumOfComponent = length(Tau) ;

%% Set fitting parameters
Count = 0 ;

% if Fix1orNot0orAbs2(1) ~= 1
    Mat_A = zeros(NumOfComponent, NumOfState) ;
    K = 0 ;
    while K < NumOfState
        K = K + 1 ;
        if Fix1orNot0orAbs2orAmp3_Mat_A(K) == 1
            Mat_A(:, K) = Initial_Mat_A(:, K) ;
        elseif Fix1orNot0orAbs2orAmp3_Mat_A(K) == 3
            Count = Count + 1 ;
            Var_Amp = abs(estimates(Count)) ;
            Mat_A(:, K) = Initial_Mat_A(:, K) .* Var_Amp ;
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
% else
%     Mat_A = Initial_Mat_A ;
% end

% just for avoiding error
Var = max(max(Mat_A)) - min(min(Mat_A)) ;
Mat_A = Mat_A + (Mat_A==0).* (Var*0.0000001) ;

Mat_G = zeros(NumOfState, NumOfState) ;
K = 0 ;
while K < NumOfState
    K = K + 1 ;
    I = 0 ;
    while I < NumOfState
        I = I + 1 ;
        if Fix1orNot0orAbs2_Mat_G(I, K) == 1
            Mat_G(I, K) = Initial_Mat_G(I, K) ;
        elseif Fix1orNot0orAbs2_Mat_G(I, K) == 2
            Count = Count + 1 ;
            Mat_G(I, K) = abs(estimates(Count)) ;
        elseif Fix1orNot0orAbs2_Mat_G(I, K) == 0
            Count = Count + 1 ;
            Mat_G(I, K) = estimates(Count) ;
        elseif Fix1orNot0orAbs2_Mat_G(I, K) == 3
            if I <= K
                Count = Count + 1 ;
                Mat_G(I, K) = abs(estimates(Count)) ;
                Mat_G(K, I) = abs(estimates(Count)) ;
            end
        end
    end
end

if Fix1orNot0orAbs2_y0(1) ~= 1
    Count = Count + 1 ;
    % InitialPara = [Tau_distribution at n=1], [Tau_dis. at n=2], [Tau_dis. at n=3],...
    if Fix1orNot0orAbs2_y0(1) == 2
        y0 = abs(estimates(Count)) ;
    else
        y0 = estimates(Count) ;
    end
else
    y0 = Initial_y0 ;
end


%% model 2D map at dT = shortest

Dif_Mat_2DFDC_It = Mat_2DFDC_It ;
Var = length(Mat_2DFDC_It) ;
Dif_Mat_2DFDC_It(2:Var) = Mat_2DFDC_It(2:Var) - Mat_2DFDC_It(1:Var-1) ;
Dif_Mat_2DFDC_It(1) = Dif_Mat_2DFDC_It(2) ;
Dif_Mat_2DFDC_It = kron(Dif_Mat_2DFDC_It, Dif_Mat_2DFDC_It') ;

Mat_M_2DFLC = Mat_A * Mat_G * Mat_A' ;
Mat_M_Model = ExpCurve * Mat_M_2DFLC * ExpCurve' + y0*Dif_Mat_2DFDC_It ;

%% Kai^2 estimation
%MatCorIn = Mat0_2DFDC_cor ; % linear Mat is used
%MatIn = Mat0_2DFDC ;        % linear Mat is used

Var1 = mean(mean(Mat_2DFDC)) ;
%Var1 = max(max(Mat_2DFDC)) ;
% Var = Mat_2DFDC + (Mat_2DFDC==0)*max(max(Mat_2DFDC)) ;
% Var1 = min(min(Var)) / 1 ;
Var = ((Mat_2DFDC_cor - Mat_M_Model).^2) ./ (Mat_2DFDC+Var1) ;  % Var1 is just for avoiding error
% Var = Area .* Var ;
Imax = size(Mat_2DFDC) ;
Imax = Imax(1,1) ;
Kai2 = sum(sum(Var)) / (Imax*Imax) ;

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