% TK_RateEq_MakeExpMatrix

function [EigenValues, ExpMatrix, EigenVector, InvEigenVector] =...
    TK_RateEq_MakeExpMatrix(RateMatrix)

%RateMatrix = [1,1,1;2,2,2;3,3,3] ;

% Make ModifiedRateMatrix
Var = size(RateMatrix) ;
Imax = Var(1,1) ;
Kmax = Var(1,2) ;
Ho = zeros(Imax, Kmax) ;

I = 0 ;
while I < Imax
    I = I + 1 ;
    K = 0 ;
    while K < Kmax
        K = K + 1 ;
        Ho(I, I) = Ho(I, I) + RateMatrix(K, I) ;
    end
    Ho(I, I) = Ho(I, I) + RateMatrix(I, I) ;
end
ModifiedRateMatrix = RateMatrix - Ho ;


% Diagonalization of ModifiedRateMagtrix
[EigenVector, EigenValues] = eig(ModifiedRateMatrix) ;
InvEigenVector = inv(EigenVector) ;
Check1 = isreal(EigenVector) ;
Check2 = isreal(EigenValues) ;
if (Check1*Check2 == 0)
    'Error: EigenMatrix contains complex value'
end


% Make exp-Matrix
ExpMatrix = zeros(Imax, Kmax) ;
I = 0 ;
while I < Imax
    I = I + 1 ;
    ExpMatrix(I, I) = exp(EigenValues(I, I)) ;
end













