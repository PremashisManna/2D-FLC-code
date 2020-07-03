% TK_MyMain_Exp_FFT_FWHM 

T = 0.004 ;%* 10^-9;             % sec Sampling period 
Fs = 1/T;            % Sampling frequency

Tau = [0.01:0.01:5.01]' ;
%IRF(:)= 0 ; IRF(300) = 1 ;

Imax = 50000 ;
xdata = ( [0:Imax-1]*T )' ;

[ExpCurve] = TK_ExpMultiDeco_For2DFLC(Tau, xdata,...
    IRF, RisePoint_FL, RisePoint_IRF, RangePoint_IRF_min, RangePoint_IRF_max) ;

Var = size(ExpCurve) ;
Imax = Var(2) ;
Kmax = Var(1) ;

n = 2^nextpow2(Kmax);
f = ( Fs*(0:(n/2))/n )' ;%* 10^-9; % GHz
P = zeros(n/2+1, Imax) ;
HWHM_GHz = zeros(Imax, 1) ;

I = 0 ;
while I < Imax
    I = I + 1 ;
    Y = fft(ExpCurve(:,I), n) ;
    Var = abs(Y) .^2 /n;
    P(:,I) = Var(1:n/2+1) ;
    
    Var1 = max(P(:,I))/2 ;
    Var = P(:,I) > Var1 ;
    Var2 = sum(Var) ;
    Var = (P(Var2,I) - Var1) / (P(Var2,I) - P(Var2+1,I)) ;
    
    HWHM_GHz(I) = f(Var2) + Var*(f(Var2+1)-f(Var2)) ;
    
    %I/Imax * 100
end

HWHM_ns = 1./(2*pi .* HWHM_GHz) ;

%figure;plot(f,P(:,100)) 
