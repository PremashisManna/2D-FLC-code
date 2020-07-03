% TK_MyMain_Analyze_02

%% Input parameters

% 3 states 2D fitting
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%

% Ana_estimates_Select = 0 ;
% TrialNumber = 1 ;
% 
% if Ana_estimates_Select == 0
%     
%     Var = size(Ana_estimates) ;
%     StateAssign = Ana_estimates(1,:) ;
%     
%     FitA = Ana_estimates(2,1) ;
%     Free0orFix1orAbs2_A = [1] ;
%     
%     FitE = Ana_estimates(3,:) ;
%     Free0orFix1orAbs2_E = [1, 1, 1, 1, 1, 1, 1] ;
%     
%     FitQ = Ana_estimates(4,:) ;
%     Free0orFix1orAbs2_Q = [2, 2, 2, 2, 2, 2, 2]./1 ;
%     
%     FitRateM = Ana_estimates(5:4+Var(2),:) ;
%     Free0orFix1orAbs2_RateM = [1, 2, 2, 1, 1, 1, 1 ;...
%                                2, 1, 2, 1, 1, 1, 1 ;...
%                                2, 2, 1, 1, 1, 1, 1 ;...
%                                1, 1, 1, 1, 2, 1, 1 ;...
%                                1, 1, 1, 2, 1, 1, 1 ;...
%                                1, 1, 1, 1, 1, 1, 2 ;...
%                                1, 1, 1, 1, 1, 2, 1 ] ;
%                                Free0orFix1orAbs2_RateM = [1, 1, 1, 1, 1, 1, 1 ;...
%                                1, 1, 1, 1, 1, 1, 1 ;...
%                                1, 1, 1, 1, 1, 1, 1 ;...
%                                1, 1, 1, 1, 1, 1, 1 ;...
%                                1, 1, 1, 1, 1, 1, 1 ;...
%                                1, 1, 1, 1, 1, 1, 1 ;...
%                                1, 1, 1, 1, 1, 1, 1 ] ;
% 
%     Fity0(2,1,:) = Ana_estimates_Free0orFix1orAbs2_y0(3:5,1) ;
%     Fity0(3,1,:) = Ana_estimates_Free0orFix1orAbs2_y0(3:5,2) ;
%     Fity0(3,2,:) = Ana_estimates_Free0orFix1orAbs2_y0(3:5,3) ;
%     Fity0(1,2,:) = Fity0(2,1,:) ;
%     Fity0(1,3,:) = Fity0(3,1,:) ;
%     Fity0(2,3,:) = Fity0(3,2,:) ;
%     Free0orFix1orAbs2_y0(2,1,:) = [2, 2, 2] ;
%     Free0orFix1orAbs2_y0(3,1,:) = [2, 2, 2] ;
%     Free0orFix1orAbs2_y0(3,2,:) = [2, 2, 2] ;
%     Free0orFix1orAbs2_y0(1,2,:) = Free0orFix1orAbs2_y0(2,1,:) ;
%     Free0orFix1orAbs2_y0(1,3,:) = Free0orFix1orAbs2_y0(3,1,:) ;
%     Free0orFix1orAbs2_y0(2,3,:) = Free0orFix1orAbs2_y0(3,2,:) ;
% else
%     %%%
%     StateAssign         = [1, 2, 3, 1, 3, 1, 3] ;
%     
%     FitA                = [0.0953] ;%[0.278] ;
%     Free0orFix1orAbs2_A = [1] ;
%     
%     % epsilon
%     FitE                = [1, 1, 1, 1, 1, 1, 1] ;
%     Free0orFix1orAbs2_E = [1, 1, 1, 1, 1, 1, 1] ;
%     
%     % quantum yield
%     FitQ                = [1, 1, 1, 1, 1, 1, 1] ;
%     Free0orFix1orAbs2_Q = [2, 2, 2, 2, 2, 2, 2] ;
%     
%     FitRateM                = [0, 0.5, 0.5, 0, 0, 0, 0 ;...
%                                0.5, 0, 0, 0, 0, 0, 0 ;...
%                                0.5, 0, 0, 0, 0, 0, 0 ;...
%                                0, 0, 0, 0, 30, 0, 0 ;...
%                                0, 0, 0, 30, 0, 0, 0 ;...
%                                0, 0, 0, 0, 0, 0, 150 ;...
%                                0, 0, 0, 0, 0, 150, 0 ] ;
%     Free0orFix1orAbs2_RateM = [1, 1, 1, 1, 1, 1, 1 ;...
%                                1, 1, 1, 1, 1, 1, 1 ;...
%                                1, 1, 1, 1, 1, 1, 1 ;...
%                                1, 1, 1, 1, 1, 1, 1 ;...
%                                1, 1, 1, 1, 1, 1, 1 ;...
%                                1, 1, 1, 1, 1, 1, 1 ;...
%                                1, 1, 1, 1, 1, 1, 1 ] ;
% 
%     Var = max(StateAssign) ;
%     Fity0 = zeros(Var,Var,3) ;
%     Free0orFix1orAbs2_y0 = zeros(Var,Var,3) ;
%     Fity0(2,1,:)                = [0.04342, 0, 0] ;%[2.6*10^5, 0.9961, 1.299*10^5] * 1 ;
%     Fity0(3,1,:)                = [0.13, 0.2369, 0.1732] ;%[0.0982, 0.1909, 0.338] * 1 ;
%     Fity0(3,2,:)                = [0, 0, 0] ;%[6.416*10^-9, 1.283*10^-8, 7.558*10^-10] *1 ;
%     Fity0(1,2,:) = Fity0(2,1,:) ;
%     Fity0(1,3,:) = Fity0(3,1,:) ;
%     Fity0(2,3,:) = Fity0(3,2,:) ;
%     Free0orFix1orAbs2_y0(2,1,:) = [1, 1, 1] ;
%     Free0orFix1orAbs2_y0(3,1,:) = [1, 1, 1] ;
%     Free0orFix1orAbs2_y0(3,2,:) = [1, 1, 1] ;
%     Free0orFix1orAbs2_y0(1,2,:) = Free0orFix1orAbs2_y0(2,1,:) ;
%     Free0orFix1orAbs2_y0(1,3,:) = Free0orFix1orAbs2_y0(3,1,:) ;
%     Free0orFix1orAbs2_y0(2,3,:) = Free0orFix1orAbs2_y0(3,2,:) ;
%     
% end


% Ana_estimates_Select = 0 ;
% TrialNumber = 3 ;
% 
% if Ana_estimates_Select == 0
%     
%     Var = size(Ana_estimates) ;
%     StateAssign = Ana_estimates(1,:) ;
%     
%     FitA = Ana_estimates(2,1) ;
%     Free0orFix1orAbs2_A = [1] ;
%     
%     FitE = Ana_estimates(3,:) ;
%     Free0orFix1orAbs2_E = [1, 1, 1, 1, 1] ;
%     
%     FitQ = Ana_estimates(4,:) ;
%     Free0orFix1orAbs2_Q = [2, 2, 2, 2, 2]./1 ;
%     
%     FitRateM = Ana_estimates(5:4+Var(2),:) ;
%     Free0orFix1orAbs2_RateM = [1, 2, 2, 2, 2 ;...
%                                2, 1, 2, 2, 2 ;...
%                                2, 2, 1, 1, 1 ;...
%                                2, 2, 1, 1, 1 ;...
%                                2, 2, 1, 1, 1 ] ;
% 
%     Fity0(2,1,:) = Ana_estimates_Free0orFix1orAbs2_y0(3:5,1) ;
%     Fity0(3,1,:) = Ana_estimates_Free0orFix1orAbs2_y0(3:5,2) ;
%     Fity0(3,2,:) = Ana_estimates_Free0orFix1orAbs2_y0(3:5,3) ;
%     Fity0(1,2,:) = Fity0(2,1,:) ;
%     Fity0(1,3,:) = Fity0(3,1,:) ;
%     Fity0(2,3,:) = Fity0(3,2,:) ;
%     Free0orFix1orAbs2_y0(2,1,:) = [2, 2, 2] ;
%     Free0orFix1orAbs2_y0(3,1,:) = [2, 2, 2] ;
%     Free0orFix1orAbs2_y0(3,2,:) = [2, 2, 2] ;
%     Free0orFix1orAbs2_y0(1,2,:) = Free0orFix1orAbs2_y0(2,1,:) ;
%     Free0orFix1orAbs2_y0(1,3,:) = Free0orFix1orAbs2_y0(3,1,:) ;
%     Free0orFix1orAbs2_y0(2,3,:) = Free0orFix1orAbs2_y0(3,2,:) ;
% else
%     %%%
%     StateAssign         = [1, 2, 3, 3, 3] ;
%     
%     FitA                = [0.0953] ;%[0.278] ;
%     Free0orFix1orAbs2_A = [1] ;
%     
%     % epsilon
%     FitE                = [1, 1, 1, 1, 1] ;
%     Free0orFix1orAbs2_E = [1, 1, 1, 1, 1] ;
%     
%     % quantum yield
%     FitQ                = [1, 1, 1, 1, 1] ;
%     Free0orFix1orAbs2_Q = [2, 2, 2, 2, 2] ;
%     
%     FitRateM                = [0, 0.5, 0.3, 30, 300 ;...
%                                0.4, 0, 0, 0, 0 ;...
%                                0.4, 0, 0, 0, 0 ;...
%                                30, 0, 0, 0, 0 ;...
%                                300, 0, 0, 0, 0 ] ;
%     Free0orFix1orAbs2_RateM = [1, 2, 2, 1, 1 ;...
%                                2, 1, 1, 1, 1 ;...
%                                2, 1, 1, 1, 1 ;...
%                                1, 1, 1, 1, 1 ;...
%                                1, 1, 1, 1, 1 ] ;
% 
%     Var = max(StateAssign) ;
%     Fity0 = zeros(Var,Var,3) ;
%     Free0orFix1orAbs2_y0 = zeros(Var,Var,3) ;
%     Fity0(2,1,:)                = [0.04342, 0, 0] ;%[2.6*10^5, 0.9961, 1.299*10^5] * 1 ;
%     Fity0(3,1,:)                = [0.13, 0.2369, 0.1732] ;%[0.0982, 0.1909, 0.338] * 1 ;
%     Fity0(3,2,:)                = [0, 0, 0] ;%[6.416*10^-9, 1.283*10^-8, 7.558*10^-10] *1 ;
%     Fity0(1,2,:) = Fity0(2,1,:) ;
%     Fity0(1,3,:) = Fity0(3,1,:) ;
%     Fity0(2,3,:) = Fity0(3,2,:) ;
%     Free0orFix1orAbs2_y0(2,1,:) = [1, 1, 1] ;
%     Free0orFix1orAbs2_y0(3,1,:) = [1, 1, 1] ;
%     Free0orFix1orAbs2_y0(3,2,:) = [1, 1, 1] ;
%     Free0orFix1orAbs2_y0(1,2,:) = Free0orFix1orAbs2_y0(2,1,:) ;
%     Free0orFix1orAbs2_y0(1,3,:) = Free0orFix1orAbs2_y0(3,1,:) ;
%     Free0orFix1orAbs2_y0(2,3,:) = Free0orFix1orAbs2_y0(3,2,:) ;
%     
% end


% Ana_estimates_Select = 0 ;
% TrialNumber = 1 ;
% 
% if Ana_estimates_Select == 0
%     
%     Var = size(Ana_estimates) ;
%     StateAssign = Ana_estimates(1,:) ;
%     
%     FitA = Ana_estimates(2,1) ;
%     Free0orFix1orAbs2_A = [2] ;
%     
%     FitE = Ana_estimates(3,:) ;
%     Free0orFix1orAbs2_E = [1, 1, 1, 1, 1] ;
%     
%     FitQ = Ana_estimates(4,:) ;
%     Free0orFix1orAbs2_Q = [2, 2, 2, 2, 2] ;
%     
%     FitRateM = Ana_estimates(5:4+Var(2),:) ;
%     Free0orFix1orAbs2_RateM = [1, 2, 2, 1, 1 ;...
%                                2, 1, 2, 1, 1 ;...
%                                2, 2, 1, 1, 1 ;...
%                                1, 1, 1, 1, 2 ;...
%                                1, 1, 1, 2, 1 ] ;
% 
%     Fity0(2,1,:) = Ana_estimates_Free0orFix1orAbs2_y0(3:5,1) ;
%     Fity0(3,1,:) = Ana_estimates_Free0orFix1orAbs2_y0(3:5,2) ;
%     Fity0(3,2,:) = Ana_estimates_Free0orFix1orAbs2_y0(3:5,3) ;
%     Fity0(1,2,:) = Fity0(2,1,:) ;
%     Fity0(1,3,:) = Fity0(3,1,:) ;
%     Fity0(2,3,:) = Fity0(3,2,:) ;
%     Free0orFix1orAbs2_y0(2,1,:) = [2, 2, 2] ;
%     Free0orFix1orAbs2_y0(3,1,:) = [2, 2, 2] ;
%     Free0orFix1orAbs2_y0(3,2,:) = [2, 2, 2] ;
%     Free0orFix1orAbs2_y0(1,2,:) = Free0orFix1orAbs2_y0(2,1,:) ;
%     Free0orFix1orAbs2_y0(1,3,:) = Free0orFix1orAbs2_y0(3,1,:) ;
%     Free0orFix1orAbs2_y0(2,3,:) = Free0orFix1orAbs2_y0(3,2,:) ;
% else
%     %%%
%     StateAssign         = [1, 2, 3, 1, 3] ;
%     
%     FitA                = [0.0953] ;%[0.278] ;
%     Free0orFix1orAbs2_A = [1] ;
%     
%     % epsilon
%     FitE                = [1, 1, 1, 1, 1] ;
%     Free0orFix1orAbs2_E = [1, 1, 1, 1, 1] ;
%     
%     % quantum yield
%     FitQ                = [1, 1, 1, 1, 1] ;
%     Free0orFix1orAbs2_Q = [2, 2, 2, 2, 2] ;
%     
%     FitRateM                = [0, 1, 0.5, 0, 0 ;...
%                                1, 0, 0, 0, 0 ;...
%                                0.5, 0, 0, 0, 0 ;...
%                                0, 0, 0, 0, 150 ;...
%                                0, 0, 0, 150, 0 ] ;
%     Free0orFix1orAbs2_RateM = [1, 1, 1, 1, 1 ;...
%                                1, 1, 1, 1, 1 ;...
%                                1, 1, 1, 1, 1 ;...
%                                1, 1, 1, 1, 1 ;...
%                                1, 1, 1, 1, 1 ] ;
% 
%     Var = max(StateAssign) ;
%     Fity0 = zeros(Var,Var,3) ;
%     Free0orFix1orAbs2_y0 = zeros(Var,Var,3) ;
%     Fity0(2,1,:)                = [0.04342, 0, 0] ;%[2.6*10^5, 0.9961, 1.299*10^5] * 1 ;
%     Fity0(3,1,:)                = [0.13, 0.2369, 0.1732] ;%[0.0982, 0.1909, 0.338] * 1 ;
%     Fity0(3,2,:)                = [0, 0, 0] ;%[6.416*10^-9, 1.283*10^-8, 7.558*10^-10] *1 ;
%     Fity0(1,2,:) = Fity0(2,1,:) ;
%     Fity0(1,3,:) = Fity0(3,1,:) ;
%     Fity0(2,3,:) = Fity0(3,2,:) ;
%     Free0orFix1orAbs2_y0(2,1,:) = [1, 1, 1] ;
%     Free0orFix1orAbs2_y0(3,1,:) = [1, 1, 1] ;
%     Free0orFix1orAbs2_y0(3,2,:) = [1, 1, 1] ;
%     Free0orFix1orAbs2_y0(1,2,:) = Free0orFix1orAbs2_y0(2,1,:) ;
%     Free0orFix1orAbs2_y0(1,3,:) = Free0orFix1orAbs2_y0(3,1,:) ;
%     Free0orFix1orAbs2_y0(2,3,:) = Free0orFix1orAbs2_y0(3,2,:) ;
%     
% end

% Ana_estimates_Select = 0 ;
% TrialNumber = 5 ;
% 
% if Ana_estimates_Select == 0
%     
%     Var = size(Ana_estimates) ;
%     StateAssign = Ana_estimates(1,:) ;
%     
%     FitA = Ana_estimates(2,1) ;
%     Free0orFix1orAbs2_A = [2] ;
%     
%     FitE = Ana_estimates(3,:) ;
%     Free0orFix1orAbs2_E = [1, 1, 1, 1] ;
%     
%     FitQ = Ana_estimates(4,:) ;
%     Free0orFix1orAbs2_Q = [2, 2, 2, 2]./1 ;
%     
%     FitRateM = Ana_estimates(5:4+Var(2),:) ;
%     Free0orFix1orAbs2_RateM = [1, 2, 2, 1 ;...
%                                2, 1, 2, 1 ;...
%                                2, 2, 1, 2 ;...
%                                1, 1, 2, 1 ] ;
% 
%     Fity0(2,1,:) = Ana_estimates_Free0orFix1orAbs2_y0(3:5,1) ;
%     Fity0(3,1,:) = Ana_estimates_Free0orFix1orAbs2_y0(3:5,2) ;
%     Fity0(3,2,:) = Ana_estimates_Free0orFix1orAbs2_y0(3:5,3) ;
%     Fity0(1,2,:) = Fity0(2,1,:) ;
%     Fity0(1,3,:) = Fity0(3,1,:) ;
%     Fity0(2,3,:) = Fity0(3,2,:) ;
%     Free0orFix1orAbs2_y0(2,1,:) = [2, 2, 2] ;
%     Free0orFix1orAbs2_y0(3,1,:) = [2, 2, 2] ;
%     Free0orFix1orAbs2_y0(3,2,:) = [2, 2, 2] ;
%     Free0orFix1orAbs2_y0(1,2,:) = Free0orFix1orAbs2_y0(2,1,:) ;
%     Free0orFix1orAbs2_y0(1,3,:) = Free0orFix1orAbs2_y0(3,1,:) ;
%     Free0orFix1orAbs2_y0(2,3,:) = Free0orFix1orAbs2_y0(3,2,:) ;
% else
%     %%%
%     StateAssign         = [1, 2, 3, 3] ;
%     
%     FitA                = [0.0953] ;%[0.278] ;
%     Free0orFix1orAbs2_A = [1] ;
%     
%     % epsilon
%     FitE                = [1, 1, 1, 1] ;
%     Free0orFix1orAbs2_E = [1, 1, 1, 1] ;
%     
%     % quantum yield
%     FitQ                = [1, 1, 1, 1] ;
%     Free0orFix1orAbs2_Q = [2, 2, 2, 2] ;
%     
%     FitRateM                = [0, 0, 20, 300 ;...
%                                0.5, 0, 0, 0 ;...
%                                20, 0, 0, 0 ;...
%                                300, 0, 0, 0] ;
%     Free0orFix1orAbs2_RateM = [1, 1, 1, 1 ;...
%                                1, 1, 1, 1 ;...
%                                1, 1, 1, 1 ;...
%                                1, 1, 1, 1 ] ;
% 
%     Var = max(StateAssign) ;
%     Fity0 = zeros(Var,Var,3) ;
%     Free0orFix1orAbs2_y0 = zeros(Var,Var,3) ;
%     Fity0(2,1,:)                = [0.04342, 0, 0] ;%[2.6*10^5, 0.9961, 1.299*10^5] * 1 ;
%     Fity0(3,1,:)                = [0.13, 0.2369, 0.1732] ;%[0.0982, 0.1909, 0.338] * 1 ;
%     Fity0(3,2,:)                = [0, 0, 0] ;%[6.416*10^-9, 1.283*10^-8, 7.558*10^-10] *1 ;
%     Fity0(1,2,:) = Fity0(2,1,:) ;
%     Fity0(1,3,:) = Fity0(3,1,:) ;
%     Fity0(2,3,:) = Fity0(3,2,:) ;
%     Free0orFix1orAbs2_y0(2,1,:) = [1, 1, 1] ;
%     Free0orFix1orAbs2_y0(3,1,:) = [1, 1, 1] ;
%     Free0orFix1orAbs2_y0(3,2,:) = [1, 1, 1] ;
%     Free0orFix1orAbs2_y0(1,2,:) = Free0orFix1orAbs2_y0(2,1,:) ;
%     Free0orFix1orAbs2_y0(1,3,:) = Free0orFix1orAbs2_y0(3,1,:) ;
%     Free0orFix1orAbs2_y0(2,3,:) = Free0orFix1orAbs2_y0(3,2,:) ;
%     
% end

% Ana_estimates_Select = 0 ;
% TrialNumber = 5 ;
% 
% if Ana_estimates_Select == 0
%     
%     Var = size(Ana_estimates) ;
%     StateAssign = Ana_estimates(1,:) ;
%     
%     FitA = Ana_estimates(2,1) ;
%     Free0orFix1orAbs2_A = [2] ;
%     
%     FitE = Ana_estimates(3,:) ;
%     Free0orFix1orAbs2_E = [1, 1, 1] ;
%     
%     FitQ = Ana_estimates(4,:) ;
%     Free0orFix1orAbs2_Q = [2, 2, 2]./1 ;
%     
%     FitRateM = Ana_estimates(5:4+Var(2),:) ;
%     Free0orFix1orAbs2_RateM = [1, 2, 2 ;...
%                                2, 1, 2 ;...
%                                2, 2, 1 ] ;
% 
%     Fity0(2,1,:) = Ana_estimates_Free0orFix1orAbs2_y0(3:5,1) ;
%     Fity0(3,1,:) = Ana_estimates_Free0orFix1orAbs2_y0(3:5,2) ;
%     Fity0(3,2,:) = Ana_estimates_Free0orFix1orAbs2_y0(3:5,3) ;
%     Fity0(1,2,:) = Fity0(2,1,:) ;
%     Fity0(1,3,:) = Fity0(3,1,:) ;
%     Fity0(2,3,:) = Fity0(3,2,:) ;
%     Free0orFix1orAbs2_y0(2,1,:) = [2, 2, 2] ;
%     Free0orFix1orAbs2_y0(3,1,:) = [2, 2, 2] ;
%     Free0orFix1orAbs2_y0(3,2,:) = [2, 2, 2] ;
%     Free0orFix1orAbs2_y0(1,2,:) = Free0orFix1orAbs2_y0(2,1,:) ;
%     Free0orFix1orAbs2_y0(1,3,:) = Free0orFix1orAbs2_y0(3,1,:) ;
%     Free0orFix1orAbs2_y0(2,3,:) = Free0orFix1orAbs2_y0(3,2,:) ;
% else
%     %%%
%     StateAssign         = [1, 2, 3] ;
%     
%     FitA                = [0.0953] ;
%     Free0orFix1orAbs2_A = [1] ;
%     
%     % epsilon
%     FitE                = [1, 1, 1] ;
%     Free0orFix1orAbs2_E = [1, 1, 1] ;
%     
%     % quantum yield
%     FitQ                = [1, 1, 1] ;
%     Free0orFix1orAbs2_Q = [2, 2, 2] ;
%     
%     FitRateM                = [0, 3, 30 ;...
%                                3, 0, 0 ;...
%                                30, 0, 0 ] ;
%     Free0orFix1orAbs2_RateM = [1, 2, 2 ;...
%                                2, 1, 1 ;...
%                                2, 1, 1 ] ;
% 
%     Var = max(StateAssign) ;
%     Fity0 = zeros(Var,Var,3) ;
%     Free0orFix1orAbs2_y0 = zeros(Var,Var,3) ;
%     Fity0(2,1,:)                = [0.04342, 0, 0] ;%[2.6*10^5, 0.9961, 1.299*10^5] * 1 ;
%     Fity0(3,1,:)                = [0.13, 0.2369, 0.1732] ;%[0.0982, 0.1909, 0.338] * 1 ;
%     Fity0(3,2,:)                = [0, 0, 0] ;%[6.416*10^-9, 1.283*10^-8, 7.558*10^-10] *1 ;
%     Fity0(1,2,:) = Fity0(2,1,:) ;
%     Fity0(1,3,:) = Fity0(3,1,:) ;
%     Fity0(2,3,:) = Fity0(3,2,:) ;
%     Free0orFix1orAbs2_y0(2,1,:) = [1, 1, 1] ;
%     Free0orFix1orAbs2_y0(3,1,:) = [1, 1, 1] ;
%     Free0orFix1orAbs2_y0(3,2,:) = [1, 1, 1] ;
%     Free0orFix1orAbs2_y0(1,2,:) = Free0orFix1orAbs2_y0(2,1,:) ;
%     Free0orFix1orAbs2_y0(1,3,:) = Free0orFix1orAbs2_y0(3,1,:) ;
%     Free0orFix1orAbs2_y0(2,3,:) = Free0orFix1orAbs2_y0(3,2,:) ;
%     
% end


%%
% 2 states 2D fitting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Ana_estimates_Select = 0 ;
% TrialNumber = 1 ;
% 
% if Ana_estimates_Select == 0
%     
%     Var = size(Ana_estimates) ;
%     StateAssign = Ana_estimates(1,:) ;
%     
%     FitA = Ana_estimates(2,1) ;
%     Free0orFix1orAbs2_A = [1] ;
%     
%     FitE = Ana_estimates(3,:) ;
%     Free0orFix1orAbs2_E = [1, 1, 1, 1, 1, 1] ;
%     
%     FitQ = Ana_estimates(4,:) ;
%     Free0orFix1orAbs2_Q = [2, 2, 2, 2, 2, 2]./1 ;
%     
%     FitRateM = Ana_estimates(5:4+Var(2),:) ;
%     Free0orFix1orAbs2_RateM = [1, 2, 1, 1, 1, 1 ;...
%                                2, 1, 1, 1, 1, 1 ;...
%                                1, 1, 1, 2, 1, 1 ;...
%                                1, 1, 2, 1, 1, 1 ;...
%                                1, 1, 1, 1, 1, 2 ;...
%                                1, 1, 1, 1, 2, 1 ] ;
% 
%     Fity0(2,1,:) = Ana_estimates_Free0orFix1orAbs2_y0(3:5) ;
%     Fity0(1,2,:) = Fity0(2,1,:) ;
%     Free0orFix1orAbs2_y0(2,1,:) = [2, 2, 2]./1 ; 
%     Free0orFix1orAbs2_y0(1,2,:) = Free0orFix1orAbs2_y0(2,1,:) ;
% else
%     %%%
%     StateAssign         = [1, 2, 1, 2, 1, 2] ;
%     
%     FitA                = [1] ;
%     Free0orFix1orAbs2_A = [1] ;
%     
%     % epsilon
%     FitE                = [1, 1, 1, 1, 1, 1] ;
%     Free0orFix1orAbs2_E = [1, 1, 1, 1, 1, 1] ;
%     
%     % quantum yield
%     FitQ                = [1, 1, 1, 1, 1, 1] ;
%     Free0orFix1orAbs2_Q = [2, 2, 2, 2, 2, 2] ;
%     
%     FitRateM                = [0, 0.3, 0, 0, 0, 0 ;...
%                                0.3, 0, 0, 0, 0, 0 ;...
%                                0, 0, 0, 30, 0, 0 ;...
%                                0, 0, 30, 0, 0, 0 ;...
%                                0, 0, 0, 0, 0, 300 ;...
%                                0, 0, 0, 0, 300, 0] ;
%     Free0orFix1orAbs2_RateM = [1, 2, 1, 1, 1, 1 ;...
%                                2, 1, 1, 1, 1, 1 ;...
%                                1, 1, 1, 1, 1, 1 ;...
%                                1, 1, 1, 1, 1, 1 ;...
%                                1, 1, 1, 1, 1, 1 ;...
%                                1, 1, 1, 1, 1, 1 ] ;
% 
%     Var = max(StateAssign) ;
%     Fity0 = zeros(Var,Var,3) ;
%     Free0orFix1orAbs2_y0 = zeros(Var,Var,3) ;
%     Fity0(2,1,:)                = [0.162239332, 0.288602364, 0.224694895] ;
%     Fity0(1,2,:) = Fity0(2,1,:) ;
%     Free0orFix1orAbs2_y0(2,1,:) = [1, 1, 1] ; 
%     Free0orFix1orAbs2_y0(1,2,:) = Free0orFix1orAbs2_y0(2,1,:) ;   
% end

% Ana_estimates_Select = 0 ;
% TrialNumber = 5 ;
% 
% if Ana_estimates_Select == 0
%     
%     Var = size(Ana_estimates) ;
%     StateAssign = Ana_estimates(1,:) ;
%     
%     FitA = Ana_estimates(2,1) ;
%     Free0orFix1orAbs2_A = [1] ;
%     
%     FitE = Ana_estimates(3,:) ;
%     Free0orFix1orAbs2_E = [1, 1, 1, 1] ;
%     
%     FitQ = Ana_estimates(4,:) ;
%     Free0orFix1orAbs2_Q = [2, 2, 2, 2]./1 ;
%     
%     FitRateM = Ana_estimates(5:4+Var(2),:) ;
%     Free0orFix1orAbs2_RateM = [1, 1, 2, 2 ;...
%                                2, 1, 1, 1 ;...
%                                2, 1, 1, 1 ;...
%                                2, 1, 1, 1] ;
% 
%     Fity0(2,1,:) = Ana_estimates_Free0orFix1orAbs2_y0(3:5) ;
%     Fity0(1,2,:) = Fity0(2,1,:) ;
%     Free0orFix1orAbs2_y0(2,1,:) = [2, 2, 2] ; 
%     Free0orFix1orAbs2_y0(1,2,:) = Free0orFix1orAbs2_y0(2,1,:) ;
% else
% 
%     StateAssign         = [1, 2, 2, 2] ;
%     
%     FitA                = [1] ;
%     Free0orFix1orAbs2_A = [1] ;
%     
%     % epsilon
%     FitE                = [1, 1, 1, 1] ;
%     Free0orFix1orAbs2_E = [1, 1, 1, 1] ;
%     
%     % quantum yield
%     FitQ                = [1, 1, 1, 1] ;
%     Free0orFix1orAbs2_Q = [2, 2, 2, 2] ;
%     
%     FitRateM                = [0, 0.3, 30, 300 ;...
%                                0.3, 0, 0, 0 ;...
%                                30, 0, 0, 0 ;...
%                                300, 0, 0, 0] ;
%     Free0orFix1orAbs2_RateM = [1, 2, 1, 1 ;...
%                                2, 1, 1, 1 ;...
%                                1, 1, 1, 1 ;...
%                                1, 1, 1, 1] ;
% 
%     Var = max(StateAssign) ;
%     Fity0 = zeros(Var,Var,3) ;
%     Free0orFix1orAbs2_y0 = zeros(Var,Var,3) ;
%     Fity0(2,1,:)                = [0.162239332, 0.288602364, 0.224694895] ;
%     Fity0(1,2,:) = Fity0(2,1,:) ;
%     Free0orFix1orAbs2_y0(2,1,:) = [1, 1, 1] ; 
%     Free0orFix1orAbs2_y0(1,2,:) = Free0orFix1orAbs2_y0(2,1,:) ;
%     
% end

%%

% Ana_estimates_Select = 0 ;
% TrialNumber = 1 ;
% 
% if Ana_estimates_Select == 0
%     
%     Var = size(Ana_estimates) ;
%     StateAssign = Ana_estimates(1,:) ;
%     
%     FitA = Ana_estimates(2,1) ;
%     Free0orFix1orAbs2_A = [1] ;
%     
%     FitE = Ana_estimates(3,:) ;
%     Free0orFix1orAbs2_E = [1, 1, 1, 1] ;
%     
%     FitQ = Ana_estimates(4,:) ;
%     Free0orFix1orAbs2_Q = [2, 2, 2, 2]./1 ;
%     
%     FitRateM = Ana_estimates(5:4+Var(2),:) ;
%     Free0orFix1orAbs2_RateM = [1, 2, 1, 1 ;...
%                                2, 1, 1, 1 ;...
%                                1, 1, 1, 2 ;...
%                                1, 1, 2, 1 ] ;
% 
%     Fity0(2,1,:) = Ana_estimates_Free0orFix1orAbs2_y0(3:5) ;
%     Fity0(1,2,:) = Fity0(2,1,:) ;
%     Free0orFix1orAbs2_y0(2,1,:) = [2, 2, 2] ; 
%     Free0orFix1orAbs2_y0(1,2,:) = Free0orFix1orAbs2_y0(2,1,:) ;
% else
%     %%%
%     StateAssign         = [1, 2, 1, 2] ;
%     
%     FitA                = [1] ;
%     Free0orFix1orAbs2_A = [1] ;
%     
%     % epsilon
%     FitE                = [1, 1, 1, 1] ;
%     Free0orFix1orAbs2_E = [1, 1, 1, 1] ;
%     
%     % quantum yield
%     FitQ                = [1, 1, 1, 1] ;
%     Free0orFix1orAbs2_Q = [2, 2, 2, 2] ;
%     
%     FitRateM                = [0, 0.3, 0, 0 ;...
%                                0.3, 0, 0, 0 ;...
%                                0, 0, 0, 300 ;...
%                                0, 0, 300, 0] ;
%     Free0orFix1orAbs2_RateM = [1, 1, 1, 1 ;...
%                                1, 1, 1, 1 ;...
%                                1, 1, 1, 1 ;...
%                                1, 1, 1, 1 ] ;
% 
%     Var = max(StateAssign) ;
%     Fity0 = zeros(Var,Var,3) ;
%     Free0orFix1orAbs2_y0 = zeros(Var,Var,3) ;
%     Fity0(2,1,:)                = [0.1771, 0.3005, 0.1755] ;
%     Fity0(1,2,:) = Fity0(2,1,:) ;
%     Free0orFix1orAbs2_y0(2,1,:) = [1, 1, 1] ; 
%     Free0orFix1orAbs2_y0(1,2,:) = Free0orFix1orAbs2_y0(2,1,:) ;
%     
% end

Ana_estimates_Select = 0 ;
TrialNumber = 1 ;

if Ana_estimates_Select == 0
    
    Var = size(Ana_estimates) ;
    StateAssign = Ana_estimates(1,:) ;
    
    FitA = Ana_estimates(2,1) ;
    Free0orFix1orAbs2_A = [1] ;
    
    FitE = Ana_estimates(3,:) ;
    Free0orFix1orAbs2_E = [1, 1] ;
    
    FitQ = Ana_estimates(4,:) ;
    Free0orFix1orAbs2_Q = [2, 2] ;
    
    FitRateM = Ana_estimates(5:4+Var(2),:) ;
    Free0orFix1orAbs2_RateM = [1, 2 ;...
                               2, 1 ] ;

    Fity0(2,1,:) = Ana_estimates_Free0orFix1orAbs2_y0(3:5) ;
    Fity0(1,2,:) = Fity0(2,1,:) ;
    Free0orFix1orAbs2_y0(2,1,:) = [1, 1, 1]*2 ; 
    Free0orFix1orAbs2_y0(1,2,:) = Free0orFix1orAbs2_y0(2,1,:) ;
else

    StateAssign         = [1, 2] ;
    
    FitA                = [1] ;
    Free0orFix1orAbs2_A = [1] ;
    
    % epsilon
    FitE                = [1, 1] ;
    Free0orFix1orAbs2_E = [1, 1] ;
    
    % quantum yield
    FitQ                = [1, 1] ;
    Free0orFix1orAbs2_Q = [2, 2] ;
    
    FitRateM                = [0, 5 ;...
                               20, 0 ] ;
    Free0orFix1orAbs2_RateM = [1, 2 ;...
                               2, 1 ] ;

    Var = max(StateAssign) ;
    Fity0 = zeros(Var,Var,3) ;
    Free0orFix1orAbs2_y0 = zeros(Var,Var,3) ;
    Fity0(2,1,:)                = [0, 0, 0] ;
    Fity0(1,2,:) = Fity0(2,1,:) ;
    Free0orFix1orAbs2_y0(2,1,:) = [1, 1, 1] ; 
    Free0orFix1orAbs2_y0(1,2,:) = Free0orFix1orAbs2_y0(2,1,:) ;
    
end
% Ana_estimates_Select = 1 ;
% TrialNumber = 1 ;
% 
% if Ana_estimates_Select == 0
%     
%     Var = size(Ana_estimates) ;
%     StateAssign = Ana_estimates(1,:) ;
%     
%     FitA = Ana_estimates(2,1) ;
%     Free0orFix1orAbs2_A = [1] ;
%     
%     FitE = Ana_estimates(3,:) ;
%     Free0orFix1orAbs2_E = [1, 1] ;
%     
%     FitQ = Ana_estimates(4,:) ;
%     Free0orFix1orAbs2_Q = [2, 2] ;
%     
%     FitRateM = Ana_estimates(5:4+Var(2),:) ;
%     Free0orFix1orAbs2_RateM = [1, 2 ;...
%                                2, 1 ] ;
%                            
%     Fity0 = Ana_estimates_y0 ;
%     Free0orFix1orAbs2_y0 = [2, 2 ;...
%                             2, 2 ] ;
% else
% 
%     StateAssign         = [1, 2] ;
%     
%     FitA                = [1] ;
%     Free0orFix1orAbs2_A = [1] ;
%     
%     % epsilon
%     FitE                = [1, 1] ;
%     Free0orFix1orAbs2_E = [1, 1] ;
%     
%     % quantum yield
%     FitQ                = [1, 1] ;
%     Free0orFix1orAbs2_Q = [2, 2] ;
%     
%     FitRateM                = [0, 50 ;...
%                                5, 0 ] ;
%     Free0orFix1orAbs2_RateM = [1, 2 ;...
%                                2, 1 ] ;
% 
%     Fity0 = [0, 0 ;...
%              0, 0 ] ;
%     Free0orFix1orAbs2_y0 = [1, 1 ;...
%                             1, 1 ] ;
% end
FittedPoints = 100000 ;

%%

Mat_PhotonSum = zeros(NumOfState) ;
N = 0 ;
while N < NumOfState
    N = N + 1 ;
    Var = Gresult_Mat_A_Global(:,N) .* Tau ;
    Var = sum(Var) ;
    Mat_PhotonSum(N,N) = Var ;
    Var = Gresult_Mat_A_Global(:,N) ;
    Var = sum(Var) ;
    Mat_AmpSum(N,N) = Var ;
end

Ana_Correlation = Gresult_Mat_G_dTfun ;

T = 0 ;
while T < length(dT)
	T = T + 1 ; 
    
    Ana_Correlation(:,:,T) = (Mat_PhotonSum) * Gresult_Mat_G_dTfun(:,:,T) * (Mat_PhotonSum) ;
    %Ana_Correlation(:,:,T) = (Mat_AmpSum) * Gresult_Mat_G_dTfun(:,:,T) * (Mat_AmpSum) ;
    %Ana_Correlation(:,:,T) = 1 * Gresult_Mat_G_dTfun(:,:,T) * 1 ;
end


Tmax = length(dT) ;
% if FuncType ~= 4
%     CorrelationRatio = zeros(NumOfState, NumOfState, Tmax) ;
% else
    Ana_CorrelationRatio = zeros(NumOfState, NumOfState, Tmax, 3) ;
% end
%CorrelationRatio_N1N2 = zeros(Tmax, 1) ;
T = 0 ;
while T < Tmax
    T = T + 1 ;
        
    %% estimate CorrelationRatio at each dT
    N1 = 0 ;
    while N1 < NumOfState
        N1 = N1 + 1 ;      
        N2 = 0 ;
        while N2 < NumOfState
            N2 = N2 + 1 ;
%             if FuncType == 1
%                 %CorrelationRatio(N1, N2, T) = Correlation(N2, N2, T) / Correlation(N1, N1, T) ;
%                 CorrelationRatio(N1, N2, T) = Correlation(N2, N1, T) / Correlation(N2, N2, T) ;
%             elseif FuncType == 2
%                 CorrelationRatio(N1, N2, T) = (Correlation(N1, N2, T)+Correlation(N2, N1, T)) / (Correlation(N1, N1, T)+Correlation(N2, N2, T)) ;
%             elseif FuncType == 3
%                 CorrelationRatio(N1, N2, T) = (Correlation(N1, N2, T)*Correlation(N2, N1, T)) / (Correlation(N1, N1, T)*Correlation(N2, N2, T)) ;
%             elseif FuncType == 4
                %CorrelationRatio(N1, N2, T, 1) = Correlation(N2, N2, T) / Correlation(N1, N1, T) ;
                Ana_CorrelationRatio(N1, N2, T, 1) = Ana_Correlation(N2, N1, T) / Ana_Correlation(N2, N2, T) ;
                Ana_CorrelationRatio(N1, N2, T, 2) = (Ana_Correlation(N1, N2, T)+Ana_Correlation(N2, N1, T)) / (Ana_Correlation(N1, N1, T)+Ana_Correlation(N2, N2, T)) ;
                Ana_CorrelationRatio(N1, N2, T, 3) = (Ana_Correlation(N1, N2, T)*Ana_Correlation(N2, N1, T)) / (Ana_Correlation(N1, N1, T)*Ana_Correlation(N2, N2, T)) ;
%             end
        end
    end
    
    %CorrelationRatio_N1N2(T) = CorrelationRatio(Ana_TargetN1N2(1), Ana_TargetN1N2(2), T) ;
    %%
end

%Ana_TauRange = TauRange ;
%Ana_Correlation = Correlation ;
%Ana_CorrelationRatio = CorrelationRatio ;
%Ana_CorrelationRatio_N1N2 = CorrelationRatio_N1N2 ;



clear Var Var1 Tau_I_Min Tau_I_Max N1 N2 T Tmax K I TargetMat Extract_Range_I
clear Extract_Range_K Sum_ExtractMat TauRange Correlation CorrelationRatio CorrelationRatio_N1N2 Ana_CorrelationRatio_N1N2



