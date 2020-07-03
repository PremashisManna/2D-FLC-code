% TK_MyMain_Analyze_1DMEM

%%

Sum_RisePointIRF_center = 300 ;
Sum_RisePointIRF_width = 13 ;



%%
Sum_RisePointIRF_start = Sum_RisePointIRF_center - floor(Sum_RisePointIRF_width/2) ;
Sum_RisePointIRF_end = Sum_RisePointIRF_center + floor(Sum_RisePointIRF_width/2) ;

Var = abs(RisePoint_Result_1DFDC_EstQ(:,2) - Sum_RisePointIRF_start) ;
[Var, Imin] = min(Var) ;
Var = abs(RisePoint_Result_1DFDC_EstQ(:,2) - Sum_RisePointIRF_end) ;
[Var, Imax] = min(Var) ;

Count = 0 ;
I = Imin-1 ;
while I < Imax
    I = I + 1 ;
    
    Var = sum(RisePoint_Result_1DFDC_Mat_A(:, I)) ;
    Var = RisePoint_Result_1DFDC_Mat_A(:, I) ./ Var ;
    
    Count = Count + 1 ;
    if Count == 1 
        VarSum = Var ;
    else
        VarSum = VarSum + Var ;
    end
end
VarSum = VarSum ./ Count ;

RisePoint_Result_1DFDC_Mat_A_Ave = VarSum ;

clear Var Imin Imax VarSum I Count Sum_RisePointIRF_start Sum_RisePointIRF_end ;



