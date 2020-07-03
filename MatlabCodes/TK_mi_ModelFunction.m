% TK_mi_ModelFunction

function [mi] = TK_mi_ModelFunction(...
    Mat_A, Tau, tMax, tMin, NumOfState, mi_TypeSelect)

mi = Mat_A ;
if mi_TypeSelect == 0   % mi = average constant for each state    
    mi(:,:) = mean(mean(Mat_A)) ;
    Var = mean(Mat_A) ;
    VarA = (1 - exp(-1*(tMax-tMin) ./ Tau)).^2 ;
    K = 0 ;
    while K < NumOfState
        K = K + 1 ;
        mi(:,K) = mi(:,K) / Var(K) ;
        mi(:,K) = mi(:,K) .* VarA ;
    end
    
    
elseif mi_TypeSelect == 1   % mi = initial Tau distribution for each state
    K = 0 ;
    while K < NumOfState
        K = K + 1 ;
        Var = mean(Mat_A(:,K)) ;
        mi(:,K) = Mat_A(:,K) ./ Var ;
    end

    
elseif mi_TypeSelect == 2   % mi = initial Tau distribution for each state
    mi(:,:) = mean(mean(Mat_A)) ;
    Var = mean(Mat_A) ;
    VarA = (1 - exp(-1*(tMax-tMin) ./ Tau)).^2 ;
    K = 0 ;
    while K < NumOfState
        K = K + 1 ;
        mi(:,K) = mi(:,K) / Var(K) ;
        mi(:,K) = mi(:,K) .* VarA ;
    end
    mi1 = mi ;
    
    K = 0 ;
    while K < NumOfState
        K = K + 1 ;
        Var = mean(Mat_A(:,K)) ;
        mi(:,K) = Mat_A(:,K) ./ Var ;
    end
    mi2 = mi ;
    
    mi = (mi1 + mi2.*3) ./4 ;
    
    
elseif mi_TypeSelect == 3   % mi = initial Tau distribution for each state
    Ivec = [1:length(Tau)] ;
    K = 0 ;
    while K < NumOfState
        K = K + 1 ;
        Var1 = Ivec * Mat_A(:,K) ;
        Var2 = sum(Mat_A(:,K)) ;
        AveI = Var1 / Var2 ;
        Var1 = (Ivec - AveI) .^ 2 ;
        Var1 = Var1 * Mat_A(:,K) ;
        StdI = sqrt(Var1 / Var2) ;
        
        Var1 = 1/StdI/sqrt(2*pi) * exp(-1*((Ivec-AveI).^2) ./ (2*StdI*StdI));
        mi(:,K) = Var1 .* Var2 ;
    end
end








