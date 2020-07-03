% TK_MyMain_Analyze_02_NotRatio

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
%     Free0orFix1orAbs2_RateM = [3, 2, 2, 1, 1, 1, 1 ;...
%                                2, 3, 2, 1, 1, 1, 1 ;...
%                                2, 2, 3, 1, 1, 1, 1 ;...
%                                1, 1, 1, 3, 2, 1, 1 ;...
%                                1, 1, 1, 2, 3, 1, 1 ;...
%                                1, 1, 1, 1, 1, 3, 2 ;...
%                                1, 1, 1, 1, 1, 2, 3 ] ;
% %     Free0orFix1orAbs2_RateM = [1, 1, 1, 1, 1, 1, 1 ;...
% %                                1, 1, 1, 1, 1, 1, 1 ;...
% %                                1, 1, 1, 1, 1, 1, 1 ;...
% %                                1, 1, 1, 1, 1, 1, 1 ;...
% %                                1, 1, 1, 1, 1, 1, 1 ;...
% %                                1, 1, 1, 1, 1, 1, 1 ;...
% %                                1, 1, 1, 1, 1, 1, 1 ] ;
% 
%     Fity0 = Ana_estimates_y0 ;
%     Free0orFix1orAbs2_y0 = [2, 2, 2 ;...
%                             2, 2, 2 ;...
%                             2, 2, 2 ] ;
% else
%     %%%
%     StateAssign         = [1, 2, 3, 1, 3, 1, 3] ;
%     
%     FitA                = [1] ;%[0.278] ;
%     Free0orFix1orAbs2_A = [1] ;
%     
%     % epsilon
%     FitE                = [1, 1, 1, 1, 1, 1, 1] ;
%     Free0orFix1orAbs2_E = [1, 1, 1, 1, 1, 1, 1] ;
%     
%     % quantum yield
%     FitQ                = [1, 2, 3, 1, 3, 1, 3] ;
%     Free0orFix1orAbs2_Q = [2, 2, 2, 2, 2, 2, 2] ;
%     
%     FitRateM                = [0.5, 0.4, 50, 0, 0, 0, 0 ;...
%                                20, 0.5, 0.05, 0, 0, 0, 0 ;...
%                                200, 0.01, 0.5, 0, 0, 0, 0 ;...
%                                0, 0, 0, 0.5, 10, 0, 0 ;...
%                                0, 0, 0, 0.1, 0.5, 0, 0 ;...
%                                0, 0, 0, 0, 0, 0.5, 2 ;...
%                                0, 0, 0, 0, 0, 0.1, 0.5 ] ;
%                            FitRateM = FitRateM ./5 ;
%     Free0orFix1orAbs2_RateM = [1, 1, 1, 1, 1, 1, 1 ;...
%                                1, 1, 1, 1, 1, 1, 1 ;...
%                                1, 1, 1, 1, 1, 1, 1 ;...
%                                1, 1, 1, 1, 1, 1, 1 ;...
%                                1, 1, 1, 1, 1, 1, 1 ;...
%                                1, 1, 1, 1, 1, 1, 1 ;...
%                                1, 1, 1, 1, 1, 1, 1 ] ;
% 
%     Fity0 = [10, 0, 18.24 ;...
%              0, 0, 0 ;...
%              18.24, 0, 100 ] ;
%     Free0orFix1orAbs2_y0 = [2, 1, 2 ;...
%                             1, 2, 1 ;...
%                             2, 1, 2 ] ;
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
%     Free0orFix1orAbs2_A = [1] ;
%     
%     FitE = Ana_estimates(3,:) ;
%     Free0orFix1orAbs2_E = [1, 1, 1, 1, 1] ;
%     
%     FitQ = Ana_estimates(4,:) ;
%     Free0orFix1orAbs2_Q = [2, 2, 2, 2, 2]./1 ;
%     
%     FitRateM = Ana_estimates(5:4+Var(2),:) ;
%     Free0orFix1orAbs2_RateM = [3, 2, 2, 2, 2 ;...
%                                2, 3, 2, 2, 2 ;...
%                                2, 2, 3, 1, 1 ;...
%                                2, 2, 1, 3, 1 ;...
%                                2, 2, 1, 1, 3 ] ;
% 
%     Fity0 = Ana_estimates_y0 ;
%     Free0orFix1orAbs2_y0 = [2, 2, 2 ;...
%                             2, 2, 2 ;...
%                             2, 2, 2 ] ;
%                         
% else
%     %%%
%     StateAssign         = [1, 2, 3, 3, 3] ;
%     
%     FitA                = [1] ;%[0.278] ;
%     Free0orFix1orAbs2_A = [1] ;
%     
%     % epsilon
%     FitE                = [1, 1, 1, 1, 1] ;
%     Free0orFix1orAbs2_E = [1, 1, 1, 1, 1] ;
%     
%     % quantum yield
%     FitQ                = [1, 1, 3, 3, 3] ;
%     Free0orFix1orAbs2_Q = [2, 2, 2, 2, 2] ;
%     
%     FitRateM                = [0.05, 0.2, 0.001, 5, 100 ;...
%                                10, 0.05, 0, 0, 0 ;...
%                                0.001, 0, 0.05, 0, 0 ;...
%                                5, 0, 0, 0.05, 0 ;...
%                                5, 0, 0, 0, 0.05 ] ;
%     Free0orFix1orAbs2_RateM = [3, 1, 1, 1, 1 ;...
%                                1, 3, 1, 1, 1 ;...
%                                1, 1, 3, 1, 1 ;...
%                                1, 1, 1, 3, 1 ;...
%                                1, 1, 1, 1, 3 ] ;
% 
%     Fity0 = [10, 0, 18.24 ;...
%              0, 0, 0 ;...
%              18.24, 0, 100 ] ;
%     Free0orFix1orAbs2_y0 = [2, 1, 2 ;...
%                             1, 2, 1 ;...
%                             2, 1, 2 ] ;
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
%     Free0orFix1orAbs2_A = [1] ;
%     
%     FitE = Ana_estimates(3,:) ;
%     Free0orFix1orAbs2_E = [1, 1, 1, 1, 1] ;
%     
%     FitQ = Ana_estimates(4,:) ;
%     Free0orFix1orAbs2_Q = [2, 2, 2, 2, 2]./1 ;
%     
%     FitRateM = Ana_estimates(5:4+Var(2),:) ;
%     Free0orFix1orAbs2_RateM = [3, 2, 2, 1, 1 ;...
%                                2, 3, 2, 1, 1 ;...
%                                2, 2, 3, 1, 1 ;...
%                                1, 1, 1, 3, 2 ;...
%                                1, 1, 1, 2, 3 ] ;
% 
%     Fity0 = Ana_estimates_y0 ;
%     Free0orFix1orAbs2_y0 = [2, 2, 2 ;...
%                             2, 2, 2 ;...
%                             2, 2, 2 ] ;
% else
%     %%%
%     StateAssign         = [1, 2, 3, 1, 3] ;
%     
%     FitA                = [1] ;%[0.278] ;
%     Free0orFix1orAbs2_A = [1] ;
%     
%     % epsilon
%     FitE                = [1, 1, 1, 1, 1] ;
%     Free0orFix1orAbs2_E = [1, 1, 1, 1, 1] ;
%     
%     % quantum yield
%     FitQ                = [1, 2, 3, 1, 3] ;
%     Free0orFix1orAbs2_Q = [1, 1, 1, 1, 1] ;
%     
%     FitRateM                = [0.2, 0.2, 50, 0, 0 ;...
%                                0.2, 0.2, 0, 0, 0 ;...
%                                50, 0, 0.2, 0, 0 ;...
%                                0, 0, 0, 0.2, 100 ;...
%                                0, 0, 0, 100, 0.2 ] ;
%     Free0orFix1orAbs2_RateM = [1, 1, 1, 1, 1 ;...
%                                1, 1, 1, 1, 1 ;...
%                                1, 1, 1, 1, 1 ;...
%                                1, 1, 1, 1, 1 ;...
%                                1, 1, 1, 1, 1 ] ;
% 
%     Fity0 = [10, 0, 18.24 ;...
%              0, 0, 0 ;...
%              18.24, 0, 100 ] ;
%     Free0orFix1orAbs2_y0 = [2, 1, 2 ;...
%                             1, 2, 1 ;...
%                             2, 1, 2 ] ;
%     
% end

% Ana_estimates_Select = 0 ;
% TrialNumber = 4 ;
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
%     Free0orFix1orAbs2_RateM = [3, 2, 2, 2 ;...
%                                2, 3, 2, 2 ;...
%                                2, 2, 3, 1 ;...
%                                2, 2, 1, 3 ] ;
% 
%     Fity0 = Ana_estimates_y0 ;
%     Free0orFix1orAbs2_y0 = [2, 2, 2 ;...
%                             2, 2, 2 ;...
%                             2, 2, 2 ] ;
% else
%     %%%
%     StateAssign         = [1, 2, 3, 3] ;
%     
%     FitA                = [1] ;%[0.278] ;
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
%     FitRateM                = [0.05, 0, 20, 100 ;...
%                                0.5, 0.05, 0, 0 ;...
%                                20, 0, 0.05, 0 ;...
%                                100, 0, 0, 0.05] ;
%     Free0orFix1orAbs2_RateM = [1, 1, 1, 1 ;...
%                                1, 1, 1, 1 ;...
%                                1, 1, 1, 1 ;...
%                                1, 1, 1, 1 ] ;
% 
%     Fity0 = [10, 0, 18.24 ;...
%              0, 0, 0 ;...
%              18.24, 0, 100 ] ;
%     Free0orFix1orAbs2_y0 = [2, 1, 2 ;...
%                             1, 2, 1 ;...
%                             2, 1, 2 ] ;
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
%     Free0orFix1orAbs2_A = [1] ;
%     
%     FitE = Ana_estimates(3,:) ;
%     Free0orFix1orAbs2_E = [1, 1, 1] ;
%     
%     FitQ = Ana_estimates(4,:) ;
%     Free0orFix1orAbs2_Q = [2, 2, 2]./1 ;
%     
%     FitRateM = Ana_estimates(5:4+Var(2),:) ;
%     Free0orFix1orAbs2_RateM = [3, 2, 2 ;...
%                                2, 3, 2 ;...
%                                2, 2, 3 ] ;
% 
%     Fity0 = Ana_estimates_y0 ;
%     Free0orFix1orAbs2_y0 = [2, 2, 2 ;...
%                             2, 2, 2 ;...
%                             2, 2, 2 ] ;
% else
%     %%%
%     StateAssign         = [1, 2, 3] ;
%     
%     FitA                = [1] ;
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
%     FitRateM                = [0.05, 3, 30 ;...
%                                3, 0.05, 0 ;...
%                                30, 0, 0.05 ] ;
%     Free0orFix1orAbs2_RateM = [1, 1, 1 ;...
%                                1, 1, 1 ;...
%                                1, 1, 1 ] ;
% 
%     Fity0 = [1.4, 0.9406, 6.5 ;...
%              0.9406, 1, 0 ;...
%              6.5, 0, 54 ] ;
%     Free0orFix1orAbs2_y0 = [2, 1, 1 ;...
%                             1, 2, 1 ;...
%                             1, 1, 2 ] ;
%     
% end

%%

% Ana_estimates_Select = 0 ;
% TrialNumber = 1 ;
% FitNormalizedCor_y1n0 = 1 ;
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
%     Free0orFix1orAbs2_E = [1, 1, 1, 1, 1, 1, 1, 1] ;
%     
%     FitQ = Ana_estimates(4,:) ;
%     Free0orFix1orAbs2_Q = [2, 2, 2, 2, 2, 2, 2, 1].*1 ;
%     
%     FitRateM = Ana_estimates(5:4+Var(2),:) ;
%     Free0orFix1orAbs2_RateM = [3, 2, 1, 1, 1, 1, 1, 1 ;...
%                                2, 3, 1, 1, 1, 1, 1, 1 ;...
%                                1, 1, 3, 2, 1, 1, 1, 1 ;...
%                                1, 1, 2, 3, 1, 1, 1, 1 ;...
%                                1, 1, 1, 1, 3, 1, 1, 1 ;...
%                                1, 1, 1, 1, 1, 3, 1, 1 ;...
%                                1, 1, 1, 1, 1, 1, 3, 2 ;...
%                                1, 1, 1, 1, 1, 1, 2, 3 ] ;
% %Free0orFix1orAbs2_RateM = Free0orFix1orAbs2_RateM ./Free0orFix1orAbs2_RateM ;
%     Fity0 = Ana_estimates_y0 ;
%     Free0orFix1orAbs2_y0 = [1, 1 ;...
%                             1, 1 ] ;
% else
%     %%%
%     StateAssign         = [1, 2, 1, 2, 1, 2, 1, 2] ;
%     
%     FitA                = [1] ;
%     Free0orFix1orAbs2_A = [1] ;
%     
%     % epsilon
%     FitE                = [1, 1, 1, 1, 1, 1, 1, 1] ;
%     Free0orFix1orAbs2_E = [1, 1, 1, 1, 1, 1, 1, 1] ;
%     
%     % quantum yield
%     %FitQ                = [2, 1, 2, 9, 9, 2, 9] ;
%     FitQ                = [1, 7.9, 2.38, 15.6, 11.2, 26, 2, 9] ;
%     Free0orFix1orAbs2_Q = [2, 2, 2, 2, 2, 2, 2, 2]./1 ;
%     
%     FitRateM                = [0.01, 90, 0, 0, 0, 0, 0, 0 ;...
%                                190, 0.01, 0, 0, 0, 0, 0, 0 ;...
%                                0, 0, 0.01, 5.7, 0, 0, 0, 0 ;...
%                                0, 0, 47, 0.01, 0, 0, 0, 0 ;...
%                                0, 0, 0, 0, 0.01, 10000, 0, 0 ;...
%                                0, 0, 0, 0, 10000, 0.01, 0, 0 ;...
%                                0, 0, 0, 0, 0, 0, 0.01, 100 ;...
%                                0, 0, 0, 0, 0, 0, 1500, 0.01 ] ;
%     Free0orFix1orAbs2_RateM = [1, 1, 1, 1, 1, 1, 1, 1 ;...
%                                1, 1, 1, 1, 1, 1, 1, 1 ;...
%                                1, 1, 1, 1, 1, 1, 1, 1 ;...
%                                1, 1, 1, 1, 1, 1, 1, 1 ;...
%                                1, 1, 1, 1, 1, 1, 1, 1 ;...
%                                1, 1, 1, 1, 1, 1, 1, 1 ;...
%                                1, 1, 1, 1, 1, 1, 1, 1 ;...
%                                1, 1, 1, 1, 1, 1, 1, 1 ] ;
% 
%     Fity0 = [0, 0 ;...
%              0, 0 ] ;
%     Free0orFix1orAbs2_y0 = [1, 1 ;...
%                             1, 1 ] ;  
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
%     Free0orFix1orAbs2_A = [1] ;
%     
%     FitE = Ana_estimates(3,:) ;
%     Free0orFix1orAbs2_E = [1, 1, 1, 1, 1] ;
%     
%     FitQ = Ana_estimates(4,:) ;
%     Free0orFix1orAbs2_Q = [1, 1, 2, 2, 2].*1 ;
%     
%     FitRateM = Ana_estimates(5:4+Var(2),:) ;
%     Free0orFix1orAbs2_RateM = [1, 1, 1, 1, 1 ;...
%                                1, 1, 1, 1, 1 ;...
%                                1, 1, 1, 2, 1 ;...
%                                1, 1, 2, 1, 1 ;...
%                                1, 1, 1, 1, 2 ] ;
% %Free0orFix1orAbs2_RateM = Free0orFix1orAbs2_RateM ./Free0orFix1orAbs2_RateM ;
%     Fity0 = Ana_estimates_y0 ;
%     Free0orFix1orAbs2_y0 = [1, 2 ;...
%                             2, 1 ] ;
% else
%     %%%
%     StateAssign         = [1, 2, 1, 2, 1] ;
%     
%     FitA                = [1] ;
%     Free0orFix1orAbs2_A = [1] ;
%     
%     % epsilon
%     FitE                = [1, 1, 1, 1, 1] ;
%     Free0orFix1orAbs2_E = [1, 1, 1, 1, 1] ;
%     
%     % quantum yield
%     FitQ                = [12 11, 4, 4, 6] ;
%     Free0orFix1orAbs2_Q = [1, 1, 1, 1, 1]./1 ;
%     
%     FitRateM                = [0.2, 50, 0, 0, 0 ;...
%                                5000, 0.2, 0, 0, 0 ;...
%                                0, 0, 0.2, 10, 0 ;...
%                                0, 0, 50, 0.2, 0 ;...
%                                0, 0, 0, 0, 0.001 ] ;
%     Free0orFix1orAbs2_RateM = [1, 1, 1, 1, 1 ;...
%                                1, 1, 1, 1, 1;...
%                                1, 1, 1, 1, 1 ;...
%                                1, 1, 1, 1, 1 ;...
%                                1, 1, 1, 1, 1 ] ;
% 
%     Fity0 = [0, 12 ;...
%              12, 0 ] ;
%     Free0orFix1orAbs2_y0 = [3, 2 ;...
%                             2, 3 ] ;  
% end

%%
% 2 states 2D fitting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Ana_estimates_Select = 0 ;
% TrialNumber = 10 ;
% FitNormalizedCor_y1n0 = 0 ;
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
%     Free0orFix1orAbs2_RateM = [3, 2, 1, 1, 1, 1 ;...
%                                2, 3, 1, 1, 1, 1 ;...
%                                1, 1, 3, 2, 1, 1 ;...
%                                1, 1, 2, 3, 1, 1 ;...
%                                1, 1, 1, 1, 3, 1 ;...
%                                1, 1, 1, 1, 1, 3 ] ;
% %Free0orFix1orAbs2_RateM = Free0orFix1orAbs2_RateM ./Free0orFix1orAbs2_RateM ;
%     Fity0 = Ana_estimates_y0 ;
%     Free0orFix1orAbs2_y0 = [1, 1 ;...
%                             1, 1 ] ;
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
%     FitQ                = [10, 10, 10, 10, 10, 10] ;
%     Free0orFix1orAbs2_Q = [2, 2, 2, 2, 2, 2]./1 ;
%     
%     FitRateM                = [0, 30, 0, 0, 0, 0 ;...
%                                10, 0, 0, 0, 0, 0 ;...
%                                0, 0, 0, 100, 0, 0 ;...
%                                0, 0, 300, 0, 0, 0 ;...
%                                0, 0, 0, 0, 0, 10000 ;...
%                                0, 0, 0, 0, 10000, 0] ;
%     Free0orFix1orAbs2_RateM = [1, 1, 1, 1, 1, 1 ;...
%                                1, 1, 1, 1, 1, 1 ;...
%                                1, 1, 1, 1, 1, 1 ;...
%                                1, 1, 1, 1, 1, 1 ;...
%                                1, 1, 1, 1, 1, 1 ;...
%                                1, 1, 1, 1, 1, 1 ] ;
% 
%     Fity0 = [0, 0 ;...
%              0, 0 ] ;
%     Free0orFix1orAbs2_y0 = [1, 1 ;...
%                             1, 1 ] ;  
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
%     Free0orFix1orAbs2_A = [1] ;
%     
%     FitE = Ana_estimates(3,:) ;
%     Free0orFix1orAbs2_E = [1, 1, 1, 1] ;
%     
%     FitQ = Ana_estimates(4,:) ;
%     Free0orFix1orAbs2_Q = [2, 2, 2, 2]./1 ;
%     
%     FitRateM = Ana_estimates(5:4+Var(2),:) ;
%     Free0orFix1orAbs2_RateM = [3, 2, 2, 2 ;...
%                                2, 3, 1, 1 ;...
%                                2, 1, 3, 1 ;...
%                                2, 1, 1, 3] ;
% Free0orFix1orAbs2_RateM = Free0orFix1orAbs2_RateM ./Free0orFix1orAbs2_RateM ;
% 
%     Fity0 = Ana_estimates_y0 ;
%     Free0orFix1orAbs2_y0 = [3, 2 ;...
%                             2, 3 ] ;
% 
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
%     FitQ                = [0.4, 1, 1, 1] ;
%     Free0orFix1orAbs2_Q = [1, 2, 2, 2] ;
%     
%     FitRateM                = [0.05, 0.5, 30, 300 ;...
%                                0.5, 0.05, 0, 0 ;...
%                                30, 0, 0.05, 0 ;...
%                                300, 0, 0, 0.05] ;
%     Free0orFix1orAbs2_RateM = [1, 2, 1, 1 ;...
%                                2, 1, 1, 1 ;...
%                                1, 1, 1, 1 ;...
%                                1, 1, 1, 1] ;
% 
%     Fity0 = [2.1, 2.272 ;...
%              2.272, 9 ] ;
%     Free0orFix1orAbs2_y0 = [2, 1 ;...
%                             1, 2 ] ;
%     
% end

%%

% Ana_estimates_Select = 0 ;
% TrialNumber = 5 ;
% FitNormalizedCor_y1n0 = 0 ;
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
%                                1, 1, 1, 1 ;...
%                                1, 1, 1, 1 ] ;
% %Free0orFix1orAbs2_RateM = Free0orFix1orAbs2_RateM ./Free0orFix1orAbs2_RateM ;
% 
%     Fity0 = Ana_estimates_y0 ;
%     Free0orFix1orAbs2_y0 = [1, 1 ;...
%                             1, 1 ] ;
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
%     FitRateM                = [0, 10, 0, 0 ;...
%                                10, 0, 0, 0 ;...
%                                0, 0, 0, 10^4 ;...
%                                0, 0, 10^4, 0] ;
%     Free0orFix1orAbs2_RateM = [1, 1, 1, 1 ;...
%                                1, 1, 1, 1 ;...
%                                1, 1, 1, 1 ;...
%                                1, 1, 1, 1 ] ;
% 
%     Fity0 = [0, 0 ;...
%              0, 0 ] ;
%     Free0orFix1orAbs2_y0 = [1, 1 ;...
%                             1, 1 ] ;
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
%     Free0orFix1orAbs2_E = [1, 1, 1] ;
%     
%     FitQ = Ana_estimates(4,:) ;
%     Free0orFix1orAbs2_Q = [2, 2, 2]./1 ;
%     
%     FitRateM = Ana_estimates(5:4+Var(2),:) ;
%     Free0orFix1orAbs2_RateM = [3, 2, 2 ;...
%                                2, 3, 1 ;...
%                                2, 1, 3 ] ;
% %Free0orFix1orAbs2_RateM = Free0orFix1orAbs2_RateM ./Free0orFix1orAbs2_RateM ;
% 
%     Fity0 = Ana_estimates_y0 ;
%     Free0orFix1orAbs2_y0 = [3, 2 ;...
%                             2, 3 ] ;
% else
% 
%     StateAssign         = [1, 2, 2] ;
%     
%     FitA                = [1] ;
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
%     FitRateM                = [0.05, 10, 100 ;...
%                                10, 0.05, 0 ;...
%                                100, 0, 0.05 ] ;
%     Free0orFix1orAbs2_RateM = [1, 2, 1 ;...
%                                2, 1, 1 ;...
%                                1, 1, 1 ] ;
% 
%     Fity0 = [2.1, 2.272 ;...
%              2.272, 9 ] ;
%     Free0orFix1orAbs2_y0 = [2, 1 ;...
%                             1, 2 ] ;
% end


Ana_estimates_Select = 0 ;
TrialNumber = 5 ;
FitNormalizedCor_y1n0 = 0 ;

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
    Free0orFix1orAbs2_RateM = [3, 2 ;...
                               2, 3 ] ;
                           
    Fity0 = Ana_estimates_y0 ;
    Free0orFix1orAbs2_y0 = [1, 1 ;...
                            1, 1 ] ;
else

    StateAssign         = [1, 2] ;
    
    FitA                = [1] ;
    Free0orFix1orAbs2_A = [1] ;
    
    % epsilon
    FitE                = [1, 1] ;
    Free0orFix1orAbs2_E = [1, 1] ;
    
    % quantum yield
    FitQ                = [1, 3] ;
    Free0orFix1orAbs2_Q = [2, 2] ;
    
    FitRateM                = [0, 30 ;...
                               10, 0 ] ;
    Free0orFix1orAbs2_RateM = [1, 1 ;...
                               1, 1 ] ;

    Fity0 = [0, 0 ;...
             0, 0 ] ;
    Free0orFix1orAbs2_y0 = [1, 1 ;...
                            1, 1 ] ;
end


%%
% load('FixA_Fit_All_Independent_3states_OnlyResult.mat')
% Ana_estimates_Select = 0 ;
% TrialNumber = 10 ;
% FitNormalizedCor_y1n0 = 1 ;
% Var_Ana_estimates = Ana_estimates_V_7 ;
% 
% if Ana_estimates_Select == 0  % 3state 
%     
% %     Var = size(Ana_estimates) ;
% %     StateAssign = Ana_estimates(1,:) ;
% %     
% %     FitA = Ana_estimates(2,1) ;
% %     Free0orFix1orAbs2_A = [2] ;
% %     
% %     FitE = Ana_estimates(3,:) ;
% %     Free0orFix1orAbs2_E = [1, 1, 1, 1, 1, 1, 1, 1] ;
% %     
% %     FitQ = Ana_estimates(4,:) ;
% %     Free0orFix1orAbs2_Q = [1, 1, 1, 1, 1, 1, 2, 2].*1 ;
% %     
% %     FitRateM = Ana_estimates(5:4+Var(2),:) ;
% %     Free0orFix1orAbs2_RateM = [1, 1, 1, 1, 1, 1, 1, 1 ;...
% %                                1, 1, 1, 1, 1, 1, 1, 1 ;...
% %                                1, 1, 1, 1, 1, 1, 1, 1 ;...
% %                                1, 1, 1, 1, 1, 1, 1, 1 ;...
% %                                1, 1, 1, 1, 1, 1, 1, 1 ;...
% %                                1, 1, 1, 1, 1, 1, 1, 1 ;...
% %                                1, 1, 1, 1, 1, 1, 2, 1 ;...
% %                                1, 1, 1, 1, 1, 1, 1, 2 ] ;
% %                            
% %     Fity0 = Ana_estimates_y0 ;
% %     Free0orFix1orAbs2_y0 = [1, 1 ;...
% %                             1, 1 ] ;
%                         
%     Var = size(Ana_estimates) ;
%     StateAssign = Ana_estimates(1,:) ;
%     
%     FitA = Ana_estimates(2,1) ;
%     Free0orFix1orAbs2_A = [2] ;
%     
%     FitE = Ana_estimates(3,:) ;
%     Free0orFix1orAbs2_E = [1, 1, 1, 1, 1, 1] ;
%     
%     FitQ = Ana_estimates(4,:) ;
%     Free0orFix1orAbs2_Q = [1, 1, 1, 1, 2, 2].*1 ;
%     
%     FitRateM = Ana_estimates(5:4+Var(2),:) ;
%     Free0orFix1orAbs2_RateM = [1, 1, 1, 1, 1, 1 ;...
%                                1, 1, 1, 1, 1, 1 ;...
%                                1, 1, 1, 1, 1, 1 ;...
%                                1, 1, 1, 1, 1, 1 ;...
%                                1, 1, 1, 1, 1, 1 ;...
%                                1, 1, 1, 1, 1, 1 ] ;
%                            
%     Fity0 = Ana_estimates_y0 ;
%     Free0orFix1orAbs2_y0 = [3, 2 ;...
%                             2, 3 ] ;
%                         
% elseif Ana_estimates_Select == 1
%     
% %     StateAssign         = [1, 2, 1, 2, 1, 2, 1, 2] ;
% %     
% %     FitA                = [1] ;
% %     Free0orFix1orAbs2_A = [2] ;
% %     
% %     % epsilon
% %     FitE                = [1, 1, 1, 1, 1, 1, 1, 1] ;
% %     Free0orFix1orAbs2_E = [1, 1, 1, 1, 1, 1, 1, 1] ;
% %     
% %     % quantum yield
% %     FitQ             = [Var_Ana_estimates(4 , :),Var_Ana_estimates(4 , 5:6)] ;
% %     Free0orFix1orAbs2_Q = [1, 1, 1, 1, 1, 1, 2, 2].*1 ;
% %     
% %     Var = length(StateAssign) ;
% %     FitRateM = zeros(Var) ;
% %     FitRateM(1:Var-2,1:Var-2)   = Var_Ana_estimates(5:10 , :) ;
% %     FitRateM(Var-1:Var,Var-1:Var) = Var_Ana_estimates(9:10 , 5:6) ;
% %     Free0orFix1orAbs2_RateM = [1, 1, 1, 1, 1, 1, 1, 1 ;...
% %                                1, 1, 1, 1, 1, 1, 1, 1 ;...
% %                                1, 1, 1, 1, 1, 1, 1, 1 ;...
% %                                1, 1, 1, 1, 1, 1, 1, 1 ;...
% %                                1, 1, 1, 1, 1, 1, 1, 1 ;...
% %                                1, 1, 1, 1, 1, 1, 1, 1 ;...
% %                                1, 1, 1, 1, 1, 1, 1, 1 ;...
% %                                1, 1, 1, 1, 1, 1, 1, 1 ] ;
% % 
% %     Fity0 = [0, 0 ;...
% %              0, 0 ] ;
% %     Free0orFix1orAbs2_y0 = [1, 1 ;...
% %                             1, 1 ] ;
%                         
%     StateAssign         = [1, 2, 1, 2, 1, 2] ;
%     
%     FitA                = [1] ;
%     Free0orFix1orAbs2_A = [2] ;
%     
%     % epsilon
%     FitE                = [1, 1, 1, 1, 1, 1] ;
%     Free0orFix1orAbs2_E = [1, 1, 1, 1, 1, 1] ;
%     
%     % quantum yield
%     FitQ                = Var_Ana_estimates(4 , :) ;
%     Free0orFix1orAbs2_Q = [1, 1, 1, 1, 1, 1].*1 ;
%     
%     FitRateM                = Var_Ana_estimates(5:10 , :) ;
%     Free0orFix1orAbs2_RateM = [1, 1, 1, 1, 1, 1 ;...
%                                1, 1, 1, 1, 1, 1 ;...
%                                1, 1, 1, 1, 1, 1 ;...
%                                1, 1, 1, 1, 1, 1 ;...
%                                1, 1, 1, 1, 1, 1 ;...
%                                1, 1, 1, 1, 1, 1 ] ;
% 
%     Fity0 = [0, 0 ;...
%              0, 0 ] ;
%     Free0orFix1orAbs2_y0 = [1, 1 ;...
%                             1, 1 ] ;
%                         
% elseif Ana_estimates_Select == 100   % 3state Fix
%     
%     Var = size(Ana_estimates) ;
%     StateAssign = Ana_estimates(1,:) ;
%     
%     FitA = Ana_estimates(2,1) ;
%     Free0orFix1orAbs2_A = [2] ;
%     
%     FitE = Ana_estimates(3,:) ;
%     Free0orFix1orAbs2_E = [1, 1, 1, 1, 1, 1] ;
%     
%     FitQ = Ana_estimates(4,:) ;
%     Free0orFix1orAbs2_Q = [1, 1, 1, 1, 1, 1] ;
%     
%     FitRateM = Ana_estimates(5:4+Var(2),:) ;
%     Free0orFix1orAbs2_RateM = [1, 1, 1, 1, 1, 1 ;...
%                                1, 1, 1, 1, 1, 1 ;...
%                                1, 1, 1, 1, 1, 1 ;...
%                                1, 1, 1, 1, 1, 1 ;...
%                                1, 1, 1, 1, 1, 1 ;...
%                                1, 1, 1, 1, 1, 1 ] ;
%                            
%     Fity0 = Ana_estimates_y0 ;
%     Free0orFix1orAbs2_y0 = [1, 1 ;...
%                             1, 1 ] ;
%                         
% elseif Ana_estimates_Select == 101
% 
%     StateAssign         = [1, 2, 1, 2, 1, 2] ;
%     
%     FitA                = [1] ;
%     Free0orFix1orAbs2_A = [2] ;
%     
%     % epsilon
%     FitE                = [1, 1, 1, 1, 1, 1] ;
%     Free0orFix1orAbs2_E = [1, 1, 1, 1, 1, 1] ;
%     
%     % quantum yield
%     FitQ                = Var_Ana_estimates(4 , :) ;
%     Free0orFix1orAbs2_Q = [1, 1, 1, 1, 1, 1] ;
%     
%     FitRateM                = Var_Ana_estimates(5:10 , :) ;
%     Free0orFix1orAbs2_RateM = [1, 1, 1, 1, 1, 1 ;...
%                                1, 1, 1, 1, 1, 1 ;...
%                                1, 1, 1, 1, 1, 1 ;...
%                                1, 1, 1, 1, 1, 1 ;...
%                                1, 1, 1, 1, 1, 1 ;...
%                                1, 1, 1, 1, 1, 1 ] ;
% 
%     Fity0 = [0, 0 ;...
%              0, 0 ] ;
%     Free0orFix1orAbs2_y0 = [1, 1 ;...
%                             1, 1 ] ;
%                         
% elseif Ana_estimates_Select == 200  % 2state Fix
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
%     Free0orFix1orAbs2_Q = [1, 1, 1, 1] ;
%     
%     FitRateM = Ana_estimates(5:4+Var(2),:) ;
%     Free0orFix1orAbs2_RateM = [1, 1, 1, 1 ;...
%                                1, 1, 1, 1 ;...
%                                1, 1, 1, 1 ;...
%                                1, 1, 1, 1 ] ;
%                            
%     Fity0 = Ana_estimates_y0 ;
%     Free0orFix1orAbs2_y0 = [1, 1 ;...
%                             1, 1 ] ;
%                         
% elseif Ana_estimates_Select == 201
%     
%     StateAssign         = [1, 2, 1, 2] ;
%     
%     FitA                = [1] ;
%     Free0orFix1orAbs2_A = [2] ;
%     
%     % epsilon
%     FitE                = [1, 1, 1, 1] ;
%     Free0orFix1orAbs2_E = [1, 1, 1, 1] ;
%     
%     % quantum yield
%     FitQ                = Var_Ana_estimates(4 , 1:4) ;
%     Free0orFix1orAbs2_Q = [1, 1, 1, 1] ;
%     
%     FitRateM                = Var_Ana_estimates(5:8 , 1:4) ;
%     Free0orFix1orAbs2_RateM = [1, 1, 1, 1 ;...
%                                1, 1, 1, 1 ;...
%                                1, 1, 1, 1 ;...
%                                1, 1, 1, 1 ] ;
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
    
    Ana_Correlation(:,:,T) = (Mat_PhotonSum)' * Gresult_Mat_G_dTfun(:,:,T) * (Mat_PhotonSum) ;
    %Ana_Correlation(:,:,T) = (Mat_AmpSum)' * Gresult_Mat_G_dTfun(:,:,T) * (Mat_AmpSum) ;
    %Ana_Correlation(:,:,T) = 1 * Gresult_Mat_G_dTfun(:,:,T) * 1 ;
end



clear Var Var1 Tau_I_Min Tau_I_Max N1 N2 T Tmax K I TargetMat Extract_Range_I
clear Extract_Range_K Sum_ExtractMat TauRange Correlation CorrelationRatio CorrelationRatio_N1N2 Ana_CorrelationRatio_N1N2



