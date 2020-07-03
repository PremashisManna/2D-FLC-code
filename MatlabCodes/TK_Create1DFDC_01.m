% test

function [Mat_1DFDC_lin, Mat_2DFDC_lint, Mat_1DFDC_log, Mat_2DFDC_logt] = TK_Create1DFDC_01...
    (tt1, kin1, Tstart, Tend, tMin, tMax, tStep, lint_BinFactor, logt_Imax)
%%
% % macrotime should be in "tt1"
% % microtime should be in "kin1"

% dT = 0.1 ;   % sec
% ddT = 0.05 ;   % sec
% 
% Tstart = 0 ;    % sec
% Tend = 10^8 ;   % sec
% 
% tMin = 1;%1.9    ;   % ns
% tMax = 12 ;   % ns
% tStep = 0.004 ; % ns
% 
% lint_BinFactor = 2 ;   % data points in linear scale = t_Imax / lint_BinFactor
% logt_Imax = 100 ;   % data points in log scale
% 
Seekbar_Step = 0.90 ; % seek bar
%%

Imax = length(tt1) ;
logt_Imax = logt_Imax + 1 ;

t_Imax = ceil((tMax-tMin)/tStep) + lint_BinFactor ;
lint_Imax = ceil(t_Imax / lint_BinFactor) ;
t_Imax = lint_Imax * lint_BinFactor ;
Mat_2DFDC_lin = zeros(lint_Imax, lint_Imax) ;
Mat_2DFDC_lint = (lint_BinFactor*tStep) * [0:lint_Imax-1]' ;
Mat_2DFDC_log = zeros(logt_Imax, logt_Imax) ;
Mat_2DFDC_logt = t_Imax .^ ([0:logt_Imax-1]'/(logt_Imax-1)) * tStep - tStep ;

% Var = Tend - (dT + ddT/2) ;
Var = Tend ;
Seekbar_Imax = sum((Tstart<=tt1)&(tt1<=Var)) ;
Seekbar = Seekbar_Imax * Seekbar_Step ;
Seekbar_Count = 1 ;

NonCount_I = 0 ;
NonCount_K = 0 ;
Count_I = 0 ;
Count_K = 0 ;

Check_K = 0 ;
I = 0 ;

tic
while I < Imax % t1 search
    I = I + 1 ;
    if (tt1(I) >= Tstart) && (tt1(I) <= Tend)
        tauI = kin1(I) - round(tMin/tStep) ;
        
        if (tauI <= 0) || (tauI >= t_Imax)
            NonCount_I = NonCount_I + 1 ;
            continue
        end
        
        lint_I = ceil(tauI/lint_BinFactor) ;
        
%         if (tauI > 1)
%             logt_I = logt_Imax * log(tauI) / log(t_Imax) ;   % data storage point for t1
%         else
%             logt_I = logt_Imax * log(tauI+0.0001) / log(t_Imax) ;   % data storage point for t1
%         end
%         logt_I = ceil(logt_I) ;
        T = 1 ;
        while T < logt_Imax
            T = T + 1 ;
            if (tauI*tStep) <= Mat_2DFDC_logt(T)
                logt_I = T - 1 ;
                break
            end
            if T == logt_Imax
                logt_I = -1 ;
            end
        end
        
%         dTstart = tt1(I) + dT - ddT/2 ;
%         dTend = tt1(I) + dT + ddT/2 ;
        dTstart = tt1(I) ; 
        dTend = tt1(I) ; 
        
        if (dTend > tt1(Imax)) || (dTend > Tend)
            NonCount_I = NonCount_I + 1 ;
            break
        end
        Count_I = Count_I + 1 ;
        
%         if Check_K == 0
%             K = I ;%- 1 ;
%         else
%             if Start_K > I
%                 K = Start_K - 1 ;
%             else
%                 K = I ;
%             end
%         end
%         Check_K = 0 ;

        K = I - 1 ;
        while K < Imax  % t2 search
            K = K + 1 ;
            
%             if (tt1(K) >= dTstart) && (Check_K == 0)
%                 Check_K = 1 ;
%                 Start_K = K ;
%             end
            
            if (tt1(K) >= dTstart) && (tt1(K) <= dTend)              
                tauK = kin1(K) - round(tMin/tStep) ;
                
                if (tauK <= 0) || (tauK >= t_Imax)
                    NonCount_K = NonCount_K + 1 ;
                    continue
                end
                lint_K = ceil(tauK/lint_BinFactor) ;
                Mat_2DFDC_lin(lint_I, lint_K) = Mat_2DFDC_lin(lint_I, lint_K) + 1 ;
                
%                 if (tauK > 1)
%                     logt_K = logt_Imax * log(tauK) / log(t_Imax) ;   % data storage point for t2
%                 else
%                     logt_K = logt_Imax * log(tauK+0.0001) / log(t_Imax) ;   % data storage point for t2
%                 end
%                 logt_K = ceil(logt_K) ;
                T = 1 ;
                while T < logt_Imax
                    T = T + 1 ;
                    if (tauK*tStep) <= Mat_2DFDC_logt(T)
                        logt_K = T - 1 ;
                        break
                    end
                    if T == logt_Imax
                        logt_K = -1 ;
                    end
                end

                if (logt_I > 0) && (logt_K > 0)
                    Mat_2DFDC_log(logt_I, logt_K) = Mat_2DFDC_log(logt_I, logt_K) + 1 ;
                end
                
                Count_K = Count_K + 1 ;
                         
            elseif (tt1(K) > dTend)
              %  [K, I]
              %  K = Imax ;
                break
            end
        end    
    
    
        % progress check
        if Count_I >= Seekbar
            Var = round(Seekbar_Count * Seekbar_Step * 100) ;
            display( strcat(num2str(Var), ' % ___ ', num2str(toc), ' sec') )  % done /%
            Seekbar_Count = Seekbar_Count + 1 ;
            Seekbar = Seekbar_Count * Seekbar_Step * Seekbar_Imax ;
            display(' ')
        end
    end
end


Mat_1DFDC_lin = diag(Mat_2DFDC_lin) ;
Var = length(Mat_1DFDC_lin)-1 ;
Mat_1DFDC_lin = Mat_1DFDC_lin(1:Var) ;
Mat_2DFDC_lint = Mat_2DFDC_lint(1:Var) ;

Mat_2DFDC_logt = Mat_2DFDC_logt(1:logt_Imax-1) ;
Mat_1DFDC_log = diag(Mat_2DFDC_log) ;
Mat_1DFDC_log = Mat_1DFDC_log(1:logt_Imax-1) ;

end


