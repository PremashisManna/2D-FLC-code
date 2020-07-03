% TK_MyMain_Create2DFDC_cor_SeparateData_BootStrap_v02

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

BootStrap_PhotonNum_factor = 1 ; % (default 1)   If set at N, pick up molecules to reach a photon number of N x TotalPhotonNum 
BootStrap_GroupNum_factor = 1 ; % (default 1)  If set at N(2), pick up molecules from (set of all molecules) x N

%DataFolderPath = '~/Dropbox (MIT)/2D_FLCS/SimulationData/PhotonNumber/' ;
DataFolderPath = 'C:\Users\Steph\Dropbox (MIT)\LiposomeAnalysis\2D_FLC_liposome\' ;
%DataFolderPath = 'C:\Users\ACER-PC\Dropbox (MIT)\Data_MIT\2D_FLCS_Data\SimulatedData\analysis_8\' ;

DataListFromExternalFile = 0;

if DataListFromExternalFile == 1
    FileName_DataList =  strcat(DataFolderPath, 'DataNameList_LHCSR_V_7p5_path.txt') ;
else
DataNameList = 'L7000_data_tr';
end

% Fluorescence Decay Correlation parameter
UseCor1orNot0 = 1 ; % Use or not cor-matrix (after subtraction of longest-matrix to remove uncorrelated contribution)
%%%%%%%%%%%%%% dT = [10^-5, 10^-4, 10^-3, 10^-2, 10^-1, 10^0, 10^1] ;    % sec 
%%%%%%%%%%%%%% ddT = 10^-4;   % sec
dT_long = max(dT) ;%10^1 ;  % sec

tMin = 0.5;     % ns
tMax = 12.2 ;   % ns
tStep = 0.004 ; % ns

tt1_Step = 1e-3; %10e-3; % s either 1 (for simulated data) 1e-3 for experimental data

Equi0orNonequi1 = 0 ;

lint_BinFactor = 4 ;   % data points in linear scale = t_Imax / lint_BinFactor
logt_Imax = 100 ;   % data points in log scale


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if DataListFromExternalFile == 1
    Var = strcat(FileName_DataList) ;
    Fid = fopen(Var,'rt') ;
    IIIImax = 0 ;
    while (feof(Fid) == 0)
        IIIImax = IIIImax + 1 ;
        Line = fgetl(Fid) ;
        DataNameLength(IIIImax) = length(Line) ;
        DataNameList(IIIImax, 1:DataNameLength(IIIImax)) = char(Line) ;
    end
    fclose(Fid) ;
else 
    Var = size(DataNameList) ;
    IIIImax = Var(1) ;
end

Var1 = 0 ;
IIII = 0 ;
while IIII < IIIImax
    IIII = IIII + 1 ;
    FileName = strcat(DataFolderPath, DataNameList(IIII, :), '.mat') ;
    load(FileName) ;

    Var = size(tt1) ;
    if Var(1) > Var(2)
        PhotonNum = Var(1) ;
    else
        PhotonNum = Var(2) ;
    end
    Var1 = Var1 + PhotonNum ;
end
TargetPhotonNum = Var1 * BootStrap_PhotonNum_factor ;

%%
MeasurementTime_s = zeros(IIIImax, 1) ;

if (BootStrap_PhotonNum_factor == 1) && (BootStrap_GroupNum_factor == 1)
    IIII_Random = [1:IIIImax] ;
else
    Var = datetime('now') ;
    [Var1,Var2,Var3] = hms(Var) ;
    Var = round(Var1)*round(Var2)*round(Var3) ;
    Var = randi(Var) ;
    Var = Var - floor(Var/1000)*1000 ;
    
    IIII = 0 ;
    while IIII < Var
        IIII = IIII + 1 ;
        Var1 = randperm(IIIImax*BootStrap_GroupNum_factor) ;
    end
    Var = (Var1./IIIImax) - floor(Var1./IIIImax) ;
    Var = Var + (Var<=0) ;
    IIII_Random = round(Var .* IIIImax) ;
end

Var_TotalPhotonNum = 0 ;
IIII = 0 ;

while IIII < (IIIImax*BootStrap_GroupNum_factor)
    IIII = IIII + 1 ;
    FileName = strcat(DataFolderPath, DataNameList(IIII_Random(IIII), :), '.mat') ;
    load(FileName) ;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%
    Var = size(tt1) ;
    if Var(1) > Var(2)
        tt1 = tt1' ;
        PhotonNum = Var(1) ;
    else
        PhotonNum = Var(2) ;
    end
    tt1 = tt1 .* tt1_Step ;
    Var_TotalPhotonNum = Var_TotalPhotonNum + PhotonNum ;
    
    Var = size(kin1);
    if Var(1) > Var(2)
        kin1 = kin1' ;
    end
    
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
    Var_Tend = Actual_Tend - (dT_longest - dT_shortest - ddT/2) ;
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
    
    %%%% 1D %%%%
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
        Var_Tend = Actual_Tend - (dT_longest - Var - ddT/2) ;
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
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%
    if IIII == 1
        Hozonn_Mat_1DFDC_lin = Mat_1DFDC_lin ;
        Hozonn_Mat_1DFDC_log = Mat_1DFDC_log ;
        
        Hozonn_Mat_2DFDC_Short_lin = Mat_2DFDC_Short_lin ;
        Hozonn_Mat_2DFDC_Short_log = Mat_2DFDC_Short_log ;
        
        Hozonn_Mat_2DFDC_Short_cor_lin = Mat_2DFDC_Short_cor_lin ;
        Hozonn_Mat_2DFDC_Short_cor_log = Mat_2DFDC_Short_cor_log ;
        
        Hozonn_Mat_2DFDC_lin = Mat_2DFDC_lin ;
        Hozonn_Mat_2DFDC_log = Mat_2DFDC_log ;
        
        Hozonn_Mat_2DFDC_cor_lin = Mat_2DFDC_cor_lin ;
        Hozonn_Mat_2DFDC_cor_log = Mat_2DFDC_cor_log ;
    else
        Hozonn_Mat_1DFDC_lin = Hozonn_Mat_1DFDC_lin + Mat_1DFDC_lin ;
        Hozonn_Mat_1DFDC_log = Hozonn_Mat_1DFDC_log + Mat_1DFDC_log ;
        
        Hozonn_Mat_2DFDC_Short_lin = Hozonn_Mat_2DFDC_Short_lin + Mat_2DFDC_Short_lin ;
        Hozonn_Mat_2DFDC_Short_log = Hozonn_Mat_2DFDC_Short_log + Mat_2DFDC_Short_log ;
        
        Hozonn_Mat_2DFDC_Short_cor_lin = Hozonn_Mat_2DFDC_Short_cor_lin + Mat_2DFDC_Short_cor_lin ;
        Hozonn_Mat_2DFDC_Short_cor_log = Hozonn_Mat_2DFDC_Short_cor_log + Mat_2DFDC_Short_cor_log ;
        
        Hozonn_Mat_2DFDC_lin = Hozonn_Mat_2DFDC_lin + Mat_2DFDC_lin ;
        Hozonn_Mat_2DFDC_log = Hozonn_Mat_2DFDC_log + Mat_2DFDC_log ;
        
        Hozonn_Mat_2DFDC_cor_lin = Hozonn_Mat_2DFDC_cor_lin + Mat_2DFDC_cor_lin ;
        Hozonn_Mat_2DFDC_cor_log = Hozonn_Mat_2DFDC_cor_log + Mat_2DFDC_cor_log ;    
    end
    
    %%
    
    Rmax = length(kin1) ;
    R = 0 ;
    CheckR = 0 ;
    while R < Rmax
        R = R + 1 ;
        CheckT = tt1(1, R) ;
        if (CheckT >= Tstart) && (CheckT <= Tend)
            %K = K + 1 ;
            %totalTime = totalTimeK + tt1(1, R) ;

            if CheckR == 0
                tt1min_R = R ;
                tt1max_R = R ;
                CheckR = 1 ;
            else
                tt1max_R = R ;
            end
        end
    end
    if CheckR == 0
        tt1min_R = 1 ;
        tt1max_R = 1 ;
    end
    %totalTimeK = totalTime + TimeInterval_BetweenData ;
    if R == 0
        MeasurementTime_s(IIII) = 0 ;
    else
        MeasurementTime_s(IIII) = tt1(1, tt1max_R) - tt1(1, tt1min_R) ;
    end
    
    display(FileName, '   done') 
    display('')
    
    if Var_TotalPhotonNum >= TargetPhotonNum
       TotalPhotonNum = Var_TotalPhotonNum ;
       MeasurementTime_s = MeasurementTime_s(1:IIII) ;
       break 
    end
end

Mat_1DFDC_lin = Hozonn_Mat_1DFDC_lin ;
Mat_1DFDC_log = Hozonn_Mat_1DFDC_log ;

Mat_2DFDC_Short_lin = Hozonn_Mat_2DFDC_Short_lin ;
Mat_2DFDC_Short_log = Hozonn_Mat_2DFDC_Short_log ;

Mat_2DFDC_Short_cor_lin = Hozonn_Mat_2DFDC_Short_cor_lin ;
Mat_2DFDC_Short_cor_log = Hozonn_Mat_2DFDC_Short_cor_log ;

Mat_2DFDC_lin = Hozonn_Mat_2DFDC_lin ;
Mat_2DFDC_log = Hozonn_Mat_2DFDC_log ;

Mat_2DFDC_cor_lin = Hozonn_Mat_2DFDC_cor_lin ;
Mat_2DFDC_cor_log = Hozonn_Mat_2DFDC_cor_log ;

clear Hozonn_Mat_1DFDC_lin Hozonn_Mat_1DFDC_log Hozonn_Mat_2DFDC_Short_lin Hozonn_Mat_2DFDC_Short_log...
        Hozonn_Mat_2DFDC_Short_cor_lin Hozonn_Mat_2DFDC_Short_cor_log Hozonn_Mat_2DFDC_lin Hozonn_Mat_2DFDC_log...
        Hozonn_Mat_2DFDC_cor_lin Hozonn_Mat_2DFDC_cor_log


clear event_num sync_num tt1 kin1 Var_TotalPhotonNum PhotonNum;

% Var = strcat('MeasurementTime_s', '.mat') ;
% save(Var, 'MeasurementTime_s') ;

clear Fid Line FileName DataFolderPath IIII IIIImax I Imax K Kmax R Rmax ans totalTime ;

clear Var_lin Var_log temp_Mat_2DFDC_lin temp_Mat_2DFDC_log temp_Mat_2DFDC_cor_lin temp_Mat_2DFDC_cor_log
clear temp_Mat_2DFDC_Short_lin temp_Mat_2DFDC_Short_log

clear event_num sync_num DataNameLength CheckR CheckT
clear I Imax Isize Need1orNot0_longest_shortest_dTs Var
clear dT_longest Actual_Tend Var_Tend