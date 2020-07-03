function [Out1, Out2] = TK_Histgram1D(Data1D, BinDiv)%, From0orMin)
% Binngin Data1D per BinDiv
% From0orMin=0 => Binned from 0 to Max
% From0orMin=1 => Binned from Min to Max

HistPoints = ceil(max(Data1D)/BinDiv)+1 - floor(min(Data1D)/BinDiv);
MaxValue = ceil(max(Data1D)/BinDiv)*BinDiv ;
MinValue = floor(min(Data1D)/BinDiv)*BinDiv ;
I_HistData = floor(Data1D./BinDiv) - floor(min(Data1D)/BinDiv) + 1;

Hist = zeros(HistPoints, 1) ;
HistX = zeros(HistPoints, 1) ;
% HistX(:,1) = linspace(0,MaxValue,HistPoints) ;
HistX(:,1) = linspace(MinValue,MaxValue,HistPoints) ;

Imax = length(Data1D) ;
I = 0 ;
while I < Imax
   I = I+1 ;
   I_Hist = I_HistData(I) ;
   Hist(I_Hist, 1) = Hist(I_Hist, 1) + 1 ;
end

% if (From0orMin == 1)
%    I_Start = abs(floor(min(Data1D)/BinDiv)) + 1 ;
%    Hist = Hist(I_Start:HistPoints, 1) ;
%    HistX = HistX(I_Start:HistPoints, 1) ;
% end

Out1 = HistX ;
Out2 = Hist ;
%plot(HistX,Hist(:,1))
end
