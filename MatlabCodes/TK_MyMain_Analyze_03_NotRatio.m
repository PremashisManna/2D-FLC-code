% TK_MyMain_Analyze_03_NotRatio

%%% Input parameters
%%
Ana_estimates_Select = 0 ;
TrialNumber = 1;
FitNormalizedCor_y1n0 = 1 ;
FitWeightFactor = 1 ;

if Ana_estimates_Select == 1
    ComponentNumber     = 2 ;
    
    Var1 = length(StateAssign) ;
    StateAssign = Ana_estimates(2,1:Var1) ;

    FitA = Ana_estimates(3,1) ;
    Free0orFix1orAbs2_A = [1] ;
    
    FitE(1,:) = Ana_estimates(4,1:Var1) ;
 %  FitE(2,:) = Ana_estimates(4,Var1+1:2*Var1) ;
   % FitE(3,:) = Ana_estimates(4,2*Var1+1:3*Var1) ;
    Free0orFix1orAbs2_E = [1, 1, 1 ;...
                           1, 1, 1] ;...
                          % 1, 1, 1] ;
                       
Free0orFix1orAbs2_E = Free0orFix1orAbs2_E ./ Free0orFix1orAbs2_E ;

    FitQ(1,:) = Ana_estimates(5,1:Var1) ;
  %  FitQ(2,:) = Ana_estimates(5,Var1+1:2*Var1) ;
    %FitQ(3,:) = Ana_estimates(5,2*Var1+1:3*Var1) ;
    Free0orFix1orAbs2_Q = [2, 2, 2 ;...
                           2, 2, 2 ];...
                          % 2, 2, 2] ;
                       
%Free0orFix1orAbs2_Q = Free0orFix1orAbs2_Q ./ Free0orFix1orAbs2_Q ;

    FitRateM(:,:,1) = Ana_estimates(6:5+Var1,1:Var1) ;
    FitRateM(:,:,2) = Ana_estimates(6:5+Var1,Var1+1:2*Var1) ;
    %FitRateM(:,:,3) = Ana_estimates(6:5+Var1,2*Var1+1:3*Var1) ;
    Free0orFix1orAbs2_RateM(:,:,1) =   [2, 2, 2 ;...
                                        2, 2, 2;...
                                        2, 2, 2] ;
    Free0orFix1orAbs2_RateM(:,:,2) =   [2, 2, 2 ;...
                                        2, 2, 2;...
                                        2, 2, 2] ;
    %Free0orFix1orAbs2_RateM(:,:,3) =   [5, 2 ;...
                                     %   2, 5 ] ;    
                                    
%Free0orFix1orAbs2_RateM = Free0orFix1orAbs2_RateM ./Free0orFix1orAbs2_RateM ;

    Fity0 = Ana_estimates_y0 ;
    Free0orFix1orAbs2_y0 = [1, 1, 1;...
                            1, 1, 1;...
                            1, 1, 1 ] ;
else
    ComponentNumber     = 2 ;
    
    StateAssign         = [1, 2] ;
    Var = length(StateAssign) ;
    
    FitA                = [1] ;
    Free0orFix1orAbs2_A = [1] ;
    
    % epsilon
    FitE                = [1, 1 ;...
                           1, 1] ;...
                         %  1, 1, 1];
    Free0orFix1orAbs2_E = [1, 1;...
                           1, 1] ;...
                         % 1, 1, 1];
    
    % quantum yield
    
     FitQ                = [1, 1 ;...
                            1, 1];...
                           % 1.5, 7, 0];
    Free0orFix1orAbs2_Q = [2, 2;...
                           2, 2];...
                         %  2, 2, 2];
                       
    FitRateM = zeros(Var, Var, 3) ;
    FitRateM(:,:,1)         = [0, 10;...        % 1/s,  Rate Matrix Knm (n; row, m; columm)
                               10, 0] ;
                                
    FitRateM(:,:,2)         = [0, 5;...        % 1/s,  Rate Matrix Knm (n; row, m; columm)
                               5, 0.2];
                               
   % FitRateM(:,:,3)         = [0.02, 600 ;...
                             %  100, 0.02 ] ;
                           
    Free0orFix1orAbs2_RateM = zeros(Var, Var, 3) ;
    Free0orFix1orAbs2_RateM(:,:,1) = [2, 2;...
                           2, 2];
                           
    Free0orFix1orAbs2_RateM(:,:,2) =  [2, 2;...
                           2, 2];
                           
   % Free0orFix1orAbs2_RateM(:,:,3) =   [1, 1 ;...
                                      %  1, 1 ] ;
    Fity0 = [0, 0, 0;...
             0, 0, 0;...
             0, 0, 0];
    Free0orFix1orAbs2_y0 = [1, 1, 1 ;...
                           1, 1, 1 ;...
                           1, 1, 1];
end

%%
% Ana_estimates_Select = 0 ;
% TrialNumber = 1 ;
% FitNormalizedCor_y1n0 = 1 ;
% FitWeightFactor = 10 ;
% 
% if Ana_estimates_Select == 0
%     ComponentNumber     = 2 ;
% 
%     Var1 = length(StateAssign) ;
%     %Var = size(Ana_estimates) ;
%     StateAssign = Ana_estimates(2,1:Var1) ;
%     
%     FitA = Ana_estimates(3,1) ;
%     Free0orFix1orAbs2_A = [1] ;
%     
%     FitE(1,:) = Ana_estimates(4,1:Var1) ;
%     FitE(2,:) = Ana_estimates(4,Var1+1:2*Var1) ;
%     Free0orFix1orAbs2_E = [1, 1 ;...
%                            1, 1] ;
%     
%     FitQ(1,:) = Ana_estimates(5,1:Var1) ;
%     FitQ(2,:) = Ana_estimates(5,Var1+1:2*Var1) ;
%     Free0orFix1orAbs2_Q = [2, 2 ;...
%                            2, 2] ;
% %Free0orFix1orAbs2_Q = Free0orFix1orAbs2_Q ./ Free0orFix1orAbs2_Q ;
%                        
%     FitRateM(:,:,1) = Ana_estimates(6:5+Var1,1:Var1) ;
%     FitRateM(:,:,2) = Ana_estimates(6:5+Var1,Var1+1:2*Var1) ;
%     Free0orFix1orAbs2_RateM(:,:,1) =   [4, 2 ;...
%                                         2, 4 ] ;
%     Free0orFix1orAbs2_RateM(:,:,2) =   [3, 2 ;...
%                                         2, 3 ] ;
% %Free0orFix1orAbs2_RateM = Free0orFix1orAbs2_RateM ./Free0orFix1orAbs2_RateM ;
% 
%     Fity0 = Ana_estimates_y0 ;
%     Free0orFix1orAbs2_y0 = [1, 1 ;...
%                             1, 1 ] ;
% else
%     ComponentNumber     = 2 ;
% 
%     StateAssign         = [1, 2] ;
%     Var = length(StateAssign) ;
%     
%     FitA                = [1] ;
%     Free0orFix1orAbs2_A = [1] ;
%     
%     % epsilon
%     FitE                = [1, 1 ;...
%                            1, 1] ;
%     Free0orFix1orAbs2_E = [1, 1 ;...
%                            1, 1] ;
%     
%     % quantum yield
%     FitQ                = [2, 4 ;...
%                            1, 10] ;
%     Free0orFix1orAbs2_Q = [2, 2 ;...
%                            2, 2] ;
%     
%     FitRateM = zeros(Var, Var, 2) ;
%     FitRateM(:,:,1)         = [0.02, 0.01 ;...
%                                0.01, 0.02 ] ;
%     FitRateM(:,:,2)         = [0.01, 180 ;...
%                                1, 0.01 ] ;
%     Free0orFix1orAbs2_RateM = zeros(Var, Var, 2) ;
%     Free0orFix1orAbs2_RateM(:,:,1) =   [1, 1 ;...
%                                         1, 1 ] ;
%     Free0orFix1orAbs2_RateM(:,:,2) =   [1, 1 ;...
%                                         1, 1 ] ;
%     Fity0 = [0, 0 ;...
%              0, 0 ] ;
%     Free0orFix1orAbs2_y0 = [1, 1 ;...
%                             1, 1 ] ;
% 
% end
% %%
% 
% Ana_estimates_Select = 0 ;
% TrialNumber = 100;
% FitNormalizedCor_y1n0 = 1 ; %1
% FitWeightFactor = 1 ;
% 
% if Ana_estimates_Select == 1
%     ComponentNumber     = 3 ;
%     
%     Var1 = length(StateAssign) ;
%     StateAssign = Ana_estimates(2,1:Var1) ;
% 
%     FitA = Ana_estimates(3,1) ;
%     Free0orFix1orAbs2_A = [1] ;
%     
%     FitE(1,:) = Ana_estimates(4,1:Var1) ;
%    FitE(2,:) = Ana_estimates(4,Var1+1:2*Var1) ;
%     FitE(3,:) = Ana_estimates(4,2*Var1+1:3*Var1) ;
%     Free0orFix1orAbs2_E = [1, 1;%];% ;...
%                          1, 1 ;...
%                            1, 1] ;
%                        
% Free0orFix1orAbs2_E = Free0orFix1orAbs2_E ./ Free0orFix1orAbs2_E ;
% 
%     FitQ(1,:) = Ana_estimates(5,1:Var1) ;
%     FitQ(2,:) = Ana_estimates(5,Var1+1:2*Var1) ;
%     FitQ(3,:) = Ana_estimates(5,2*Var1+1:3*Var1) ;
%     Free0orFix1orAbs2_Q = [2, 2;% ;...
%                            2, 2 ;...
%                           2, 2] ;
%                        
% Free0orFix1orAbs2_Q = Free0orFix1orAbs2_Q ./ Free0orFix1orAbs2_Q ;
% 
%     FitRateM(:,:,1) = Ana_estimates(6:5+Var1,1:Var1) ;
%     FitRateM(:,:,2) = Ana_estimates(6:5+Var1,Var1+1:2*Var1) ;
%     FitRateM(:,:,3) = Ana_estimates(6:5+Var1,2*Var1+1:3*Var1) ;
%     Free0orFix1orAbs2_RateM(:,:,1) =   [3, 2 ;...
%                                         2, 3 ] ;
%     Free0orFix1orAbs2_RateM(:,:,2) =   [4, 2 ;...
%                                         2, 4 ] ;
%     Free0orFix1orAbs2_RateM(:,:,3) =   [5, 2 ;...
%                                         2, 5 ] ;    
%                                     
% Free0orFix1orAbs2_RateM = Free0orFix1orAbs2_RateM ./Free0orFix1orAbs2_RateM ;
% 
%     Fity0 = Ana_estimates_y0 ;
%     Free0orFix1orAbs2_y0 = [6, 2 ;...
%                             2, 6 ] ;
% else
%     ComponentNumber     = 3 ;
%     
%     StateAssign         = [1, 2] ;
%     Var = length(StateAssign) ;
%     
%     FitA                = [1] ;
%     Free0orFix1orAbs2_A = [1] ;
%     
%     % epsilon
%     FitE                = [1, 1 ;%;...
%                            1, 1 ;...
%                            1, 1] ;
%     Free0orFix1orAbs2_E = [1, 1 ;%;...
%                           1, 1 ;...
%                           1, 1] ;
%     
%     % quantum yield
%     FitQ                = [1, 1 ;...
%                            1, 1 ;...
%                            1, 1] ;
%     Free0orFix1orAbs2_Q = [2, 2 ;...
%                            2, 2 ;...
%                            2, 2] ;
%     
%     FitRateM = zeros(Var, Var, 1) ;
%     FitRateM(:,:,1)         = [0, 10 ;...
%                                30, 0 ] ;
%     FitRateM(:,:,2)         = [0, 100 ;...
%                                300, 0 ] ;
%     FitRateM(:,:,3)         = [0.02, 600 ;...
%                               100, 0.02 ] ;
%                            
%     Free0orFix1orAbs2_RateM = zeros(Var, Var, 1) ;
%     Free0orFix1orAbs2_RateM(:,:,1) =   [1, 1 ;...
%                                         1, 1 ] ;
%     Free0orFix1orAbs2_RateM(:,:,2) =   [1, 1 ;...
%                                         1, 1 ] ;
%     Free0orFix1orAbs2_RateM(:,:,3) =   [1, 1 ;...
%                                         1, 1 ] ;
%     Fity0 = [0, 0 ;...
%              0, 0 ] ;
%     Free0orFix1orAbs2_y0 = [3, 2 ;...
%                             2, 3 ] ;
% end

%%
% Ana_estimates_Select = 1 ;
% TrialNumber = 5;
% FitNormalizedCor_y1n0 = 1 ; %1
% FitWeightFactor = 1 ;
% 
% if Ana_estimates_Select == 1
%     ComponentNumber     = 1 ;
%     
%     Var1 = length(StateAssign) ;
%     StateAssign = Ana_estimates(2,1:Var1) ;
% 
%     FitA = Ana_estimates(3,1) ;
%     Free0orFix1orAbs2_A = [1] ;
%     
%     FitE(1,:) = Ana_estimates(4,1:Var1) ;
%   % FitE(2,:) = Ana_estimates(4,Var1+1:2*Var1) ;
%    % FitE(3,:) = Ana_estimates(4,2*Var1+1:3*Var1) ;
%     Free0orFix1orAbs2_E = [1, 1, 1 ;%;...
%                            1, 1, 1 ;...
%                            1, 1, 1] ;
%                        
% Free0orFix1orAbs2_E = Free0orFix1orAbs2_E ./ Free0orFix1orAbs2_E ;
% 
%     FitQ(1,:) = Ana_estimates(5,1:Var1) ;
%   %  FitQ(2,:) = Ana_estimates(5,Var1+1:2*Var1) ;
%    % FitQ(3,:) = Ana_estimates(5,2*Var1+1:3*Var1) ;
%     Free0orFix1orAbs2_Q = [2, 2, 2;...
%                            2, 2, 2 ;...
%                            2, 2, 2] ;
%                        
% Free0orFix1orAbs2_Q = Free0orFix1orAbs2_Q ./ Free0orFix1orAbs2_Q ;
% 
%     FitRateM(:,:,1) = Ana_estimates(6:5+Var1,1:Var1) ;
%    % FitRateM(:,:,2) = Ana_estimates(6:5+Var1,Var1+1:2*Var1) ;
%    % FitRateM(:,:,3) = Ana_estimates(6:5+Var1,2*Var1+1:3*Var1) ;
%     Free0orFix1orAbs2_RateM(:,:,1) =   [2, 2, 2 ;...
%                                         2, 2, 2;...
%                                         2, 2, 2];
%    % Free0orFix1orAbs2_RateM(:,:,2) =   [4, 2 ;...
%                                     %    2, 4 ] ;
%     %Free0orFix1orAbs2_RateM(:,:,3) =   [5, 2 ;...
%                                     %    2, 5 ] ;    
%                                     
% Free0orFix1orAbs2_RateM = Free0orFix1orAbs2_RateM ./Free0orFix1orAbs2_RateM ;
% 
%     Fity0 = Ana_estimates_y0 ;
%     Free0orFix1orAbs2_y0 = [4, 2, 2 ;...
%                             2, 4, 2;...
%                             2, 2, 4];
% else
%     ComponentNumber     = 1 ;
%     
%     StateAssign         = [1, 2, 3] ;
%     Var = length(StateAssign) ;
%     
%     FitA                = [1] ;
%     Free0orFix1orAbs2_A = [1] ;
%     
%     % epsilon
%     FitE                = [1, 1, 1 ;%;...
%                            1, 1, 1 ;...
%                            1, 1, 1] ;
%     Free0orFix1orAbs2_E = [1, 1, 1 ;%;...
%                            1, 1, 1 ;...
%                            1, 1, 1] ;
%     
%     % quantum yield
%     FitQ                = [1, 1, 1 ;%;...
%                            1, 1, 1 ;...
%                            1, 1, 1] ;
%     Free0orFix1orAbs2_Q = [2, 2, 2 ;...
%                            2, 2, 2 ;...
%                            2, 2, 2] ;
%     
%     FitRateM = zeros(Var, Var, 1) ;
%     FitRateM(:,:,1)         = [0, 10, 200 ;...
%                                20, 0, 20 ;...
%                                800, 20, 0];
% %     FitRateM(:,:,2)         = [0, 700 ;...
% %                                70, 0 ] ;
%     %FitRateM(:,:,3)         = [0.02, 600 ;...
%                              % 100, 0.02 ] ;
%                            
%     Free0orFix1orAbs2_RateM = zeros(Var, Var, 1) ;
%     Free0orFix1orAbs2_RateM(:,:,1) =   [3, 2, 2 ;...
%                                         2, 3, 2;...
%                                         2, 2, 3] ;
% %     Free0orFix1orAbs2_RateM(:,:,2) =   [1, 1 ;...
% %                                         1, 1 ] ;
%     %Free0orFix1orAbs2_RateM(:,:,3) =   [1, 1 ;...
%                                         %1, 1 ] ;
%     Fity0 = [0, 0, 0;...
%              0, 0, 0;...
%              0, 0, 0 ] ;
%     Free0orFix1orAbs2_y0 = [2, 2, 2 ;...
%                             2, 2, 2 ;...
%                             2, 2, 2];
% end
% 
% 

%%

% Ana_estimates_Select = 0 ;
% TrialNumber = 10;
% FitNormalizedCor_y1n0 = 1 ; %1
% FitWeightFactor = 1 ;
% 
% if Ana_estimates_Select == 1
%     ComponentNumber     = 1 ;
%     
%     Var1 = length(StateAssign) ;
%     StateAssign = Ana_estimates(2,1:Var1) ;
% 
%     FitA = Ana_estimates(3,1) ;
%     Free0orFix1orAbs2_A = [1] ;
%     
%     FitE(1,:) = Ana_estimates(4,1:Var1) ;
%   % FitE(2,:) = Ana_estimates(4,Var1+1:2*Var1) ;
%    % FitE(3,:) = Ana_estimates(4,2*Var1+1:3*Var1) ;
%     Free0orFix1orAbs2_E = [1, 1, 1, 1, 1, 1];%];% ;...
%                        %  1, 1 ;...
%                         %   1, 1] ;
%                        
% Free0orFix1orAbs2_E = Free0orFix1orAbs2_E ./ Free0orFix1orAbs2_E ;
% 
%     FitQ(1,:) = Ana_estimates(5,1:Var1) ;
% %    FitQ(2,:) = Ana_estimates(5,Var1+1:2*Var1) ;
%   %  FitQ(3,:) = Ana_estimates(5,2*Var1+1:3*Var1) ;
%     Free0orFix1orAbs2_Q = [2, 2, 2, 2, 2, 2];% ;...
%                          %  2, 2 ;...
%                           % 2, 2] ;
%                        
% Free0orFix1orAbs2_Q = Free0orFix1orAbs2_Q ./ Free0orFix1orAbs2_Q ;
% 
%     FitRateM(:,:,1) = Ana_estimates(6:5+Var1,1:Var1) ;
%   %  FitRateM(:,:,2) = Ana_estimates(6:5+Var1,Var1+1:2*Var1) ;
%     %FitRateM(:,:,3) = Ana_estimates(6:5+Var1,2*Var1+1:3*Var1) ;
%     Free0orFix1orAbs2_RateM(:,:,1) =  [3, 2, 4, 2, 5, 2 ;...
%                                         2, 3, 2, 4, 2, 5 ;
%                                         6, 2, 7, 2, 8, 2 ;
%                                         2, 6, 2, 7, 2, 8 ;
%                                         9, 2, 10, 2, 11, 2 ;
%                                         2, 9, 2, 10, 2, 11] ;
%   %  Free0orFix1orAbs2_RateM(:,:,2) =   [4, 2 ;...
%                                   %      2, 4 ] ;
%   %  Free0orFix1orAbs2_RateM(:,:,3) =   [5, 2 ;...
%                                        % 2, 5 ] ;    
%                                     
% Free0orFix1orAbs2_RateM = Free0orFix1orAbs2_RateM ./Free0orFix1orAbs2_RateM ;
% 
%     Fity0 = Ana_estimates_y0 ;
%     Free0orFix1orAbs2_y0 = [1, 1, 1, 1, 1, 1 ;...
%                             1, 1, 1, 1, 1, 1 ;
%                             1, 1, 1, 1, 1, 1 ;
%                             1, 1 , 1 ,1, 1, 1;
%                             1, 1, 1, 1, 1, 1;
%                             1, 1, 1, 1, 1, 1] ;
% else
%     ComponentNumber     = 2 ;
%     
%     StateAssign         = [1, 2, 1, 3, 2, 3] ;
%     Var = length(StateAssign) ;
%     
%     FitA                = [1] ;
%     Free0orFix1orAbs2_A = [1] ;
%     
%     % epsilon
%     FitE                = [1, 1, 1, 1, 1, 1 ;%;...
%                           1, 1, 1, 1, 1, 1] ;...
%                           % 1, 1, 1] ;
%     Free0orFix1orAbs2_E = [1, 1, 1, 1, 1, 1 ;%;...
%                          1, 1, 1, 1, 1, 1] ;...
%                          % 1, 1, 1] ;
%     
%     % quantum yield
%     FitQ                = [1, 1, 1, 1, 1 ,1  ;...
%                           1, 1, 1, 1, 1, 1 ];...
%                           % 1, 1, 1] ;
%     Free0orFix1orAbs2_Q = [2, 2, 2, 2, 2, 2 ;...
%                           2, 2, 2, 2, 2, 2] ;...
%                            %2, 2, 2] ;
%     
%     FitRateM = zeros(Var, Var, 1) ;
% 
%     FitRateM(:,:,1) = [0, 20, 0, 0, 0, 0;
%                         10, 0, 0, 0, 0, 0;
%                        0, 0, 0, 50, 0, 0;
%                        0, 0, 100, 0, 0, 0;
%                         0, 0, 0, 0, 0, 500;
%                         0, 0, 0, 0, 100, 0];
%                     
%                     
%     FitRateM(:,:,2) = [0, 1, 0, 0, 0, 0;
%                         0.1, 0, 0, 0, 0, 0;
%                        0, 0, 0, 0.1, 0, 0;
%                        0, 0, 0.1, 0, 0, 0;
%                         0, 0, 0, 0, 0, 0.1;
%                         0, 0, 0, 0, 0.1, 0];
% %                         
% %     FitRateM(:,:,1)         = [0.4, 20, 0.01, 0.01, 0.01, 0.01 ;...
% %                                200, 0.4, 0.01, 0.01, 0.01, 0.01 ;
% %                                0.01, 0.01, 0.2, 20, 0.01, 0.01;
% %                                0.01, 0.01, 200, 0.2, 0.01, 0.01 ;
% %                                0.01, 0.01,  0.01, 0.01, 0.2, 100;
% %                                0.01, 0.01, 0.01, 0.01, 50, 0.2] ;
% %     FitRateM(:,:,2)         = [0.4, 20, 0.01, 0.01, 0.01, 0.01 ;...
% %                                1, 0.4, 0.01, 0.01, 0.01, 0.01 ;
% %                                0.01, 0.01, 0.2, 200, 0.01, 0.01;
% %                                0.01, 0.01, 20, 0.2, 0.01, 0.01 ;
% %                                0.01, 0.01,  0.01, 0.01, 0.2, 100;
% %                                0.01, 0.01, 0.01, 0.01, 50, 0.2] ;
%   %  FitRateM(:,:,3)         = [0.02, 600, 0 ;...
%                         %      30, 0, 0 ;
%                            %    0, 0, 0] ;
%                          %  
%     Free0orFix1orAbs2_RateM = zeros(Var, Var, 1) ;
%      
%     Free0orFix1orAbs2_RateM(:,:,1) =   [2, 2, 2, 2, 2, 2 ;...
%                                         2, 2, 2, 2, 2, 2 ;
%                                         2, 2, 2, 2, 2, 2 ;
%                                         2, 2, 2, 2, 2, 2 ;
%                                         2, 2, 2, 2, 2, 2 ;
%                                         2, 2, 2, 2, 2, 2] ;
%     Free0orFix1orAbs2_RateM(:,:,2) =  [2, 2, 2, 2, 2, 2 ;...
%                                         2, 2, 2, 2, 2, 2 ;
%                                         2, 2, 2, 2, 2, 2 ;
%                                         2, 2, 2, 2, 2, 2 ;
%                                         2, 2, 2, 2, 2, 2 ;
%                                         2, 2, 2, 2, 2, 2] ;
%    % Free0orFix1orAbs2_RateM(:,:,3) =   [2, 2, 2 ;...
%                                        % 2, 2, 2 ;
%                                         %2, 2, 2] ;
%     Fity0 = [0, 0, 0, 0, 0 ,0 ;...
%              0, 0, 0, 0, 0, 0;
%              0, 0, 0, 0, 0, 0;
%              0, 0, 0, 0, 0, 0;
%              0, 0, 0, 0, 0, 0;
%              0, 0, 0, 0, 0, 0] ;
%    
%      Free0orFix1orAbs2_y0 =[1, 1, 1, 1, 1, 1  ;...
%                             1, 1, 1, 1, 1, 1;
%                             1, 1, 1, 1, 1, 1;
%                             1, 1, 1, 1, 1, 1;
%                             1, 1, 1, 1, 1, 1
%                             1, 1, 1, 1, 1, 1];
% %     Free0orFix1orAbs2_y0 = [2, 2, 2, 2, 2, 2  ;...
% %                             2, 2, 2, 2, 2, 2 ; 
% %                             2, 2, 2, 2, 2, 2;
% %                             2, 2, 2, 2, 2, 2;
% %                             2, 2, 2, 2, 2, 2];
%                            
% end
%%

% Ana_estimates_Select = 0 ;
% %%%%%%%%%%%%%% TrialNumber = 1 ;
% %%%%%%%%%%%%%% FitNormalizedCor_y1n0 = 1 ;
% 
% if Ana_estimates_Select == 0
%     ComponentNumber     = 3 ;
%     
%     Var1 = length(StateAssign) ;
%     StateAssign = Ana_estimates(2,1:Var1) ;
% 
%     FitA = Ana_estimates(3,1) ;
%     Free0orFix1orAbs2_A = [1] ;
%     
%     FitE(1,:) = Ana_estimates(4,1:Var1) ;
%     FitE(2,:) = Ana_estimates(4,Var1+1:2*Var1) ;
%     FitE(3,:) = Ana_estimates(4,2*Var1+1:3*Var1) ;
%     
% %     Free0orFix1orAbs2_E = [2, 2 ; 3, 3 ; 4, 4] ;          
% %     Free0orFix1orAbs2_E = Free0orFix1orAbs2_E ./ Free0orFix1orAbs2_E ;
% 
%     
%     FitQ(1,:) = Ana_estimates(5,1:Var1) ;
%     FitQ(2,:) = Ana_estimates(5,Var1+1:2*Var1) ;
%     FitQ(3,:) = Ana_estimates(5,2*Var1+1:3*Var1) ;
%     
% %     Free0orFix1orAbs2_Q = [2, 2 ; 2, 2 ; 2, 2] ;                   
%     Free0orFix1orAbs2_Q = Free0orFix1orAbs2_Q ./ Free0orFix1orAbs2_Q ;

    
%     FitRateM(:,:,1) = Ana_estimates(6:5+Var1,1:Var1) ;
%     FitRateM(:,:,2) = Ana_estimates(6:5+Var1,Var1+1:2*Var1) ;
%     FitRateM(:,:,3) = Ana_estimates(6:5+Var1,2*Var1+1:3*Var1) ;
%     
%     Free0orFix1orAbs2_RateM(:,:,1) =   [5, 2 ; 2, 5 ] ;
%     Free0orFix1orAbs2_RateM(:,:,2) =   [3, 2 ; 2, 3 ] ;
%     Free0orFix1orAbs2_RateM(:,:,3) =   [4, 2 ; 2, 4 ] ;                                      
%     Free0orFix1orAbs2_RateM = Free0orFix1orAbs2_RateM ./Free0orFix1orAbs2_RateM ;

    
%     Fity0 = Ana_estimates_y0 ;
%     Free0orFix1orAbs2_y0 = [1, 1 ;...
%                             1, 1 ] ;
% end



% Ana_estimates_Select = 0 ;
% %%%%%%%%%%%%%% TrialNumber = 1 ;
% %%%%%%%%%%%%%% FitNormalizedCor_y1n0 = 1 ;
% 
% if Ana_estimates_Select == 0
%     ComponentNumber     = 2 ;
%     
%     Var1 = length(StateAssign) ;
%     StateAssign = Ana_estimates(2,1:Var1) ;
%     
%     FitA = Ana_estimates(3,1) ;
%     Free0orFix1orAbs2_A = [1] ;
%     
%     FitE(1,:) = Ana_estimates(4,1:Var1) ;
%     FitE(2,:) = Ana_estimates(4,Var1+1:2*Var1) ;
%     %Free0orFix1orAbs2_E = [1, 1 ;...
%     %                       1, 1] ;
%     
%     FitQ(1,:) = Ana_estimates(5,1:Var1) ;
%     FitQ(2,:) = Ana_estimates(5,Var1+1:2*Var1) ;
%     %Free0orFix1orAbs2_Q = [2, 2 ;...
%     %                       2, 2] ;
% %Free0orFix1orAbs2_Q = Free0orFix1orAbs2_Q ./ Free0orFix1orAbs2_Q ;
%                        
%     FitRateM(:,:,1) = Ana_estimates(6:5+Var1,1:Var1) ;
%     FitRateM(:,:,2) = Ana_estimates(6:5+Var1,Var1+1:2*Var1) ;
%     %Free0orFix1orAbs2_RateM(:,:,1) =   [3, 2 ;...
%     %                                    2, 3 ] ;
%     %Free0orFix1orAbs2_RateM(:,:,2) =   [4, 2 ;...
%     %                                    2, 4 ] ;
% %Free0orFix1orAbs2_RateM = Free0orFix1orAbs2_RateM ./Free0orFix1orAbs2_RateM ;
% 
%     Fity0 = Ana_estimates_y0 ;
%     Free0orFix1orAbs2_y0 = [1, 1 ;...
%                             1, 1 ] ;
% end


FittedPoints = 100000 ;

%%

Mat_PhotonSum = zeros(NumOfState) ;
N = 0 ;
while N < NumOfState
    N = N + 1 ;
    
    Var = Gresult_Mat_A_Global(10:end,N) .* Tau(10:end) ;
    Var = sum(Var) ;
    Mat_PhotonSum(N,N) = Var ;
    
    Var = Gresult_Mat_A_Global(10:end,N) ;
    Var = sum(Var) ;
    Mat_AmpSum(N,N) = Var ;
end

Ana_Correlation = Gresult_Mat_G_dTfun ;

T = 0 ;
while T < length(dT)
	T = T + 1 ; 
    
    %Ana_Correlation(:,:,T) = (Mat_PhotonSum)' * Gresult_Mat_G_dTfun(:,:,T) * (Mat_PhotonSum) ;
    Ana_Correlation(:,:,T) = (Mat_AmpSum)' * Gresult_Mat_G_dTfun(:,:,T) * (Mat_AmpSum) ;
    %Ana_Correlation(:,:,T) = 1 * Gresult_Mat_G_dTfun(:,:,T) * 1 ;
end



clear Var Var1 Tau_I_Min Tau_I_Max N1 N2 T Tmax K I TargetMat Extract_Range_I
clear Extract_Range_K Sum_ExtractMat TauRange Correlation CorrelationRatio CorrelationRatio_N1N2 Ana_CorrelationRatio_N1N2



