% TK_MyMain_OpenAllFiles_01

%%
DataFolderPath = 'C:\Users\tkondo\Desktop\Study\Data\2016\160325_LHCSR_Alldata\Data\' ;

FileName_DataList =  'DataNameList_LHCSR_V_7p5_path.txt' ;
%FileName_DataList =  'DataNameList_LHCSR_Z_7p5_path.txt' ;
%FileName_DataList =  'DataNameList_LHCSR_V_5_path.txt' ;
%FileName_DataList =  'DataNameList_LHCSR_Z_5_path.txt' ;
%FileName_DataList =  'DataNameList_LHCII_7p5_Remove_i005_all.txt' ;
%FileName_DataList =  'DataNameList_LHCSR_II_5_Remove_k003_path.txt' ;

DataFolderPath = 'C:\Users\tkondo\Desktop\Study\Data\180604_AllData_LHCSR3\LHCSR3\WT_Vio_7p5\' ;
FileName_DataList =  strcat(DataFolderPath, 'DataNameList_LHCSR3_w_V_7p5_path.txt') ;

% Memo
% T1:	11 - 31 s	(20 s)
% T2:	31 - 61 s	(30 s)
% T3:	61 - 111 s	(50 s)
% T4:	111 - 311 s	(200 s)
% T5:	311 - 611 s	(300 s)
% T6:	611 - /////	(// s)

% T_Range = [11, 31 ;... 
%           31, 61 ;...
%           61, 111 ;...
%           111, 311 ;...
%           311, 611 ;...
%           611, 10^20 ;...
%           0, 10^20 ;...
%           11, 311 ;...
%           311, 10^20] ;% sec
      
% T_Range = [11, 31 ;... 
%           31, 61 ;...
%           61, 111 ;...
%           111, 311 ;...
%           311, 611 ;...
%           611, 10^20 ;...
%           0, 10^20 ;...
%           11, 61 ;...
%           61, 311 ;...
%           311, 10^20] ;% sec

% T_Range = [11, 31 ;... 
%           31, 61 ;...
%           61, 111 ;...
%           111, 311 ;...
%           311, 611 ;...
%           611, 10^20 ] ;% sec
      
T_Range = [0, 10^20 ] ;% sec

%%
Var = strcat(FileName_DataList) ;
Fid = fopen(Var,'rt') ;
IIIImax = 0 ;
%Kmax = 0 ;
while (feof(Fid) == 0)
    IIIImax = IIIImax + 1 ;
    Line = fgetl(Fid) ;
    DataNameLength(IIIImax) = length(Line) ;
    DataNameList(IIIImax, 1:DataNameLength(IIIImax)) = char(Line) ;
end
fclose(Fid) ;


Var = size(T_Range) ;
Time_Imax = Var(1) ;


%% Estimate total photon number of each molecule at each time interval
PhotonNumberPerMol = zeros(IIIImax, Time_Imax) ;
MeasurementTimePerMol = zeros(IIIImax, Time_Imax) ;
MolNumber = zeros(Time_Imax, 1) ;


IIII = 0 ;
while IIII < IIIImax
    IIII = IIII + 1 ;
    FileName = strcat(DataFolderPath, DataNameList(IIII, :), '.mat') ;
    load(FileName) ;
    
    %%%%%%%%%% write function below %%%%%%%%%%%%%%%
    
    
    Imax = length(tt1) ;
    tt1Max = tt1(Imax) ;
    
    Time_I = 0 ;
    while Time_I < Time_Imax
        Time_I = Time_I + 1 ;
        Tstart = T_Range(Time_I, 1) ;
        Tend = T_Range(Time_I, 2) ;
        
        Check_ReachTend = -1 ;
        I = 0 ;
        while I < Imax
            I = I + 1 ;
            if (tt1(I) >= Tstart) && (tt1(I) < Tend)
                PhotonNumberPerMol(IIII, Time_I) = PhotonNumberPerMol(IIII, Time_I) + 1 ;
                Check_ReachTend = 0 ;
            elseif (tt1(I) >= Tend)
                Check_ReachTend = 1 ;
            end
        end
        
        if Check_ReachTend == -1
            Var = 0 ;
        elseif Check_ReachTend == 0
            Var = tt1Max - Tstart ;
        elseif Check_ReachTend == 1
            Var = Tend - Tstart ;
        end
        MeasurementTimePerMol(IIII, Time_I) = Var ;
        
        if (Check_ReachTend ~= -1)
            MolNumber(Time_I) = MolNumber(Time_I) + 1 ;
        end
    end
end


%% Estimate Intensity-Histogram at each time interval
PhotonNumberPerMol = zeros(IIIImax, Time_Imax) ;
MolNumber = zeros(Time_Imax, 1) ;

% dev_tt1 = 0.1;
% dev_Int = 10 ;   % count per dev_tt1 
% dev_tt1 = 1;
% dev_Int = 50 ;   % count per dev_tt1 
dev_tt1 = 0.01;
dev_Int = 1 ;   % count per dev_tt1 
% dev_tt1 = 0.05;
% dev_Int = 5 ;   % count per dev_tt1 

Nor_Imin = round(600*dev_tt1/dev_Int) ;
%Nor_Imax = round()

Histo_Int_Imax = round(50 * (1000 * dev_tt1 / dev_Int)) ;
Histo_Int_x = zeros(Histo_Int_Imax, 1) ;
Histo_Int_y = zeros(Histo_Int_Imax, Time_Imax) ;

%Histo_Int_y_Ch = zeros(Histo_Int_Imax, IIIImax) ;

Check_Int_x = 0 ;

IIII = 0 ;
while IIII < IIIImax
    IIII = IIII + 1 ;
    FileName = strcat(DataFolderPath, DataNameList(IIII, :), '.mat') ;
    load(FileName) ;
       
    Imax = length(tt1) ;
    tt1Max = tt1(Imax) ;
    
    Time_I = 0 ;
    while Time_I < Time_Imax
        Time_I = Time_I + 1 ;
        Tstart = T_Range(Time_I, 1) ;
        Tend = T_Range(Time_I, 2) ;
        
        Check_ReachTend = -1 ;
        I = 0 ;
        while I < Imax
            I = I + 1 ;
            if (tt1(I) >= Tstart) && (tt1(I) < Tend)
                PhotonNumberPerMol(IIII, Time_I) = PhotonNumberPerMol(IIII, Time_I) + 1 ;
                if Check_ReachTend == -1
                   VarIstart = I ; 
                end
                Check_ReachTend = 0 ;
            end
        end

        if Check_ReachTend ~= -1
            VarIend = VarIstart + PhotonNumberPerMol(IIII, Time_I) - 1 ;
            
            Var_tt1 = tt1(VarIstart:VarIend) ;
            [tempHisto_tt1_x, tempHisto_tt1_y] = TK_Histgram1D(Var_tt1, dev_tt1) ;
            tempHisto_tt1_y = vertcat(tempHisto_tt1_y, [0]) ;
            [tempHisto_Int_x, tempHisto_Int_y] = TK_Histgram1D(tempHisto_tt1_y, dev_Int) ;
            tempHisto_Int_y(1) = tempHisto_Int_y(1) -1 ;
            Var = length(tempHisto_Int_x) ;

            if max(tempHisto_Int_x) >= Check_Int_x
                Check_Int_x = max(tempHisto_Int_x) ;
                Check_Int_Imax = Var ;
                Histo_Int_x(1:Var) = tempHisto_Int_x ;
            end

            Histo_Int_y(1:Var, Time_I) = Histo_Int_y(1:Var, Time_I) + tempHisto_Int_y(1:Var) ;
            
            %Histo_Int_y_Ch(1:Var, IIII) = tempHisto_Int_y ;
        end
        
        if (Check_ReachTend ~= -1)
            MolNumber(Time_I) = MolNumber(Time_I) + 1 ;
        end
    end
    
end

Histo_Int_x = Histo_Int_x(1:Check_Int_Imax, :) ;
Histo_Int_y = Histo_Int_y(1:Check_Int_Imax, :) ;
Nor_Histo_Int_y = Histo_Int_y ;

Time_I = 0 ;
while Time_I < Time_Imax
    Time_I = Time_I + 1 ;
    %Var = max(Histo_Int_y(Nor_Imin:Check_Int_Imax, Time_I)) ;
    Var = sum(Histo_Int_y(Nor_Imin:Check_Int_Imax, Time_I)) ;
    Nor_Histo_Int_y(:, Time_I) = Histo_Int_y(:, Time_I) / Var ;
end
%%

clearvars -except PhotonNumberPerMol MeasurementTimePerMol MolNumber T_Range...
    Histo_Int_x Histo_Int_y Histo_Int_y_Ch Nor_Histo_Int_y








