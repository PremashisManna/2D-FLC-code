% TK_MyMain_ConcatenateMeasureTime

Fid = fopen('DataNameListForMeaTime.txt','rt') ;
Imax = 0 ;
Kmax = 0 ;
while (feof(Fid) == 0)
    Imax = Imax + 1 ;
    Line = fgetl(Fid) ;
    DataNameLength(Imax) = length(Line) ;
    DataNameList(Imax, 1:DataNameLength(Imax)) = char(Line) ;
    FileName = strcat(DataNameList(Imax, 1:DataNameLength(Imax)), '.mat') ;
    load(FileName) ;
    Kmax = Kmax + length(MeasurementTime_s) ;
end
fclose(Fid) ;

%Concatenated_event_num = zeros(1, Kmax) ;
%Concatenated_sync_num = zeros(1, Kmax) ;
%Concatenated_tt1 = zeros(1, Kmax) ;
%Concatenated_kin1 = zeros(1, Kmax) ;

Concatenated_MeasurementTime_s = zeros(Kmax, 1) ;

I = 0 ;
K = 0 ;
while I < Imax
    I = I + 1 ;
    FileName = strcat(DataNameList(I, 1:DataNameLength(I)), '.mat') ;
    load(FileName) ;
    
    %MeasurementTime_s(I) = max(tt1) ;
    %if I == 1
    %    totalTime = 0 ;
    %else
    %    totalTime = totalTime + MeasurementTime_s(I-1) ;
    %end
    
    Rmax = length(MeasurementTime_s) ;
    R = 0 ;
    while R < Rmax
        R = R + 1 ;
        K = K + 1 ;
        %Concatenated_event_num(1, K) = event_num(1, R) ;
        %Concatenated_sync_num(1, K) = sync_num(1, R) ;
        %Concatenated_tt1(1, K) = tt1(1, R) + totalTime ;
        %Concatenated_kin1(1, K) = kin1(1, R) ;
        Concatenated_MeasurementTime_s(K, 1) = MeasurementTime_s(R, 1) ;
    end
    strcat(FileName, '   done') 
end

clear MeasurementTime_s ;

%event_num = double(Concatenated_event_num) ;
%sync_num = double(Concatenated_sync_num) ;
%tt1 = double(Concatenated_tt1) ;
%kin1 = double(Concatenated_kin1) ;
MeasurementTime_s = double(Concatenated_MeasurementTime_s) ;
clear Concatenated_MeasurementTime_s ;

%SaveFileName = strcat('Concatenate','.mat') ;
%save(SaveFileName, 'event_num', 'sync_num', 'tt1', 'kin1') ;
SaveFileName = strcat('Concatenated_MeasurementTime_s','.mat') ;
save(SaveFileName, 'MeasurementTime_s') ;

clear Fid Line DataNameList FileName I Imax K Kmax R Rmax ans SaveFileName totalTime ;