%TK_MyMain_Create2DFDC_cor_02

%%
% macrotime should be in "tt1"
% microtime should be in "kin1"

% Fluorescence Decay Correlation parameter
UseCor1orNot0 = 1 ; % Use or not cor-matrix (after subtraction of longest-matrix to remove uncorrelated contribution)
%%%%%%%%%%%%%% dT = [10^-5, 10^-4, 10^-3, 10^-2, 10^-1, 10^0, 10^1] ;    % sec 
%%%%%%%%%%%%%% ddT = 10^-4;   % sec
dT_long = max(dT) ;%10^1 ;  % sec

Tstart = 0 ;    % sec
Tend = 10^20 ;   % sec

tMin = 1;     % ns
tMax = 12.2 ;   % ns
tStep = 0.004 ; % ns

Equi0orNonequi1 = 0 ;

lint_BinFactor = 4 ;   % data points in linear scale = t_Imax / lint_BinFactor
logt_Imax = 100 ;   % data points in log scale

%%
if Tend > max(tt1)
   Actual_Tend = max(tt1) ;
else
   Actual_Tend = Tend ; 
end

if (UseCor1orNot0 == 1) && (dT_long > max(dT))
    dT_longest = dT_long ;
else
    dT_longest = max(dT) ;
end

%% Mat_2DFDC at dT = longest for uncorrelated-subtraction
if UseCor1orNot0 == 1
    Var_Tend = Actual_Tend - (dT_longest - dT_long) ;
    [Mat_2DFDC_unc_lin, Mat_2DFDC_lin_t, Mat_2DFDC_unc_log, Mat_2DFDC_log_t] = TK_Create2DFDC_04...
        (tt1, kin1, dT_long, ddT, Tstart, Var_Tend, tMin, tMax, tStep, lint_BinFactor, logt_Imax) ;
end

%% Mat_2DFDC at dT = shortest
dT_shortest = ddT/2 ;   % 1.001 is needed for removing contribution of AutoCorrelation
Var_Tend = Actual_Tend - (dT_longest - dT_shortest) ;
[temp_Mat_2DFDC_Short_lin, Mat_2DFDC_lin_t, temp_Mat_2DFDC_Short_log, Mat_2DFDC_log_t] = TK_Create2DFDC_04...
    (tt1, kin1, dT_shortest, ddT, Tstart, Var_Tend, tMin, tMax, tStep, lint_BinFactor, logt_Imax) ;

if UseCor1orNot0 == 1
    % Subtraction of un-correlated component
    Var_lin = temp_Mat_2DFDC_Short_lin - Mat_2DFDC_unc_lin ;
    Var_log = temp_Mat_2DFDC_Short_log - Mat_2DFDC_unc_log ;
    if Equi0orNonequi1 == 0
        % To eliminate anti-symmetric noise
        Mat_2DFDC_Short_cor_lin = (Var_lin + Var_lin.')/2 ;
        Mat_2DFDC_Short_cor_log = (Var_log + Var_log.')/2 ;
    elseif Equi0orNonequi1 == 1
        Mat_2DFDC_Short_cor_lin = Var_lin ;
        Mat_2DFDC_Short_cor_log = Var_log ;
    end
end

Var_lin = temp_Mat_2DFDC_Short_lin ;
Var_log = temp_Mat_2DFDC_Short_log ;
if Equi0orNonequi1 == 0
    % To eliminate anti-symmetric noise
    Mat_2DFDC_Short_lin = (Var_lin + Var_lin.')/2 ;
    Mat_2DFDC_Short_log = (Var_log + Var_log.')/2 ;
elseif Equi0orNonequi1 == 1
    Mat_2DFDC_Short_lin = Var_lin ;
    Mat_2DFDC_Short_log = Var_log ;
end
clear Var_lin Var_log

dT_shortest = 0 ;   % 1.001 is needed for removing contribution of AutoCorrelation
Var_Tend = Actual_Tend - (dT_longest - dT_shortest) ;
[Mat_1DFDC_lin, Mat_2DFDC_lin_t, Mat_1DFDC_log, Mat_2DFDC_log_t] = TK_Create1DFDC_01...
    (tt1, kin1, Tstart, Var_Tend, tMin, tMax, tStep, lint_BinFactor, logt_Imax) ;

%% Mat_2DFDC at dT

Imax = length(dT) ;
I = 0 ;
while I < Imax
    I = I + 1 ;
    
    Var = dT(I) ;
    Var_Tend = Actual_Tend - (dT_longest - Var) ;
    [temp_Mat_2DFDC_lin, Mat_2DFDC_lin_t, temp_Mat_2DFDC_log, Mat_2DFDC_log_t] = TK_Create2DFDC_04...
        (tt1, kin1, Var, ddT, Tstart, Var_Tend, tMin, tMax, tStep, lint_BinFactor, logt_Imax) ;
    
    if UseCor1orNot0 == 1
        % Subtraction of un-correlated component
        Var_lin = temp_Mat_2DFDC_lin - Mat_2DFDC_unc_lin ;
        Var_log = temp_Mat_2DFDC_log - Mat_2DFDC_unc_log ;
        if Equi0orNonequi1 == 0
            % To eliminate anti-symmetric noise
            temp_Mat_2DFDC_cor_lin = (Var_lin + Var_lin.')/2 ;
            temp_Mat_2DFDC_cor_log = (Var_log + Var_log.')/2 ;
        elseif Equi0orNonequi1 == 1
            temp_Mat_2DFDC_cor_lin = Var_lin  ;
            temp_Mat_2DFDC_cor_log = Var_log  ;
        end

        if I == 1
            Isize = size(temp_Mat_2DFDC_cor_lin) ;
            Mat_2DFDC_cor_lin = zeros(Isize(1), Isize(2), Imax) ;
            Mat_2DFDC_cor_lin(:, :, I) = temp_Mat_2DFDC_cor_lin(:, :) ;
            Isize = size(temp_Mat_2DFDC_cor_log) ;
            Mat_2DFDC_cor_log = zeros(Isize(1), Isize(2), Imax) ;
            Mat_2DFDC_cor_log(:, :, I) = temp_Mat_2DFDC_cor_log(:, :) ;
        else
            Mat_2DFDC_cor_lin(:, :, I) = temp_Mat_2DFDC_cor_lin(:, :) ;
            Mat_2DFDC_cor_log(:, :, I) = temp_Mat_2DFDC_cor_log(:, :) ;
        end
    end
    
    Var_lin = temp_Mat_2DFDC_lin ;
    Var_log = temp_Mat_2DFDC_log ;
    if Equi0orNonequi1 == 0
        % To eliminate anti-symmetric noise
        temp_Mat_2DFDC_lin = (Var_lin + Var_lin.')/2 ;
        temp_Mat_2DFDC_log = (Var_log + Var_log.')/2 ;
    elseif Equi0orNonequi1 == 1
        temp_Mat_2DFDC_lin = Var_lin ;
        temp_Mat_2DFDC_log = Var_log ;
    end

    if I == 1
        Isize = size(temp_Mat_2DFDC_lin) ;
        Mat_2DFDC_lin = zeros(Isize(1), Isize(2), Imax) ;
        Mat_2DFDC_lin(:, :, I) = temp_Mat_2DFDC_lin(:, :) ;
        Isize = size(temp_Mat_2DFDC_log) ;
        Mat_2DFDC_log = zeros(Isize(1), Isize(2), Imax) ;
        Mat_2DFDC_log(:, :, I) = temp_Mat_2DFDC_log(:, :) ;
    else
        Mat_2DFDC_lin(:, :, I) = temp_Mat_2DFDC_lin(:, :) ;
        Mat_2DFDC_log(:, :, I) = temp_Mat_2DFDC_log(:, :) ;
    end
    
end

clear Var_lin Var_log temp_Mat_2DFDC_lin temp_Mat_2DFDC_log temp_Mat_2DFDC_cor_lin temp_Mat_2DFDC_cor_log
clear temp_Mat_2DFDC_Short_lin temp_Mat_2DFDC_Short_log

clear event_num sync_num
clear I Imax Isize Need1orNot0_longest_shortest_dTs Var
clear dT_longest Actual_Tend Var_Tend


