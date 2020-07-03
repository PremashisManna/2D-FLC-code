% TK_MyMain_Simu_PhotonStream
clc; clear; 
%%     
%%%%%%
load('IRF.mat') % load irf
TargetPhotonNum = round(1*10^5) ;    % sampling photon number
TargetTotalT = 1000;
StateTransCheckEveryT_On1orOff0 = 1 ;
% % 
 NumOfState = 2 ;        % number of states
 Int = [10000/1, 10000/1] ;    % cps   Fluorescence intensity
 Life = [1, 2] ;       % ns    Fluorescence lifetime
% 
 TransRate = [0, 10 ;...        % 1/s,  Rate Matrix Knm (n; row, m; columm)
              10, 0 ];           % Knm = rate n -> m
% %%%%%%
IRF_forSimu = IRF ;         % Set IRF for simu
IRFxdata_forSimu = xdata ;  % Note: ns units

Tstep = 10^-6 ;         % sec     simulated T step
tstep = 0.004 ;         % ns    TCSPC time resolution

TafterPhotobleach = 10 ; % sec   interval time after photobleach
RateDiff = 10^-100 ;   % 1/s,  Diffusion rate

DoCheck1orNot0 = 0 ;    % Do check (1) or Not (0)
Dwell_Tstep = 0.01 ;   % sec for dwell histogram

% noise contribution is not included in the current ver.
%LifeMax = 12.5 ;        % ns
%IntNoise = 200     ;    % cps   Noise intensity

tic
%%
% InitialState = NumOfState ;
% Current_State = InitialState ;
%Current_Int = Int(InitialState) ;
%Current_Life = Life(InitialState) ;

Stream_State = zeros(TargetPhotonNum, 1) ;
Stream_Life = zeros(TargetPhotonNum, 1) ;
Stream_detectedT = zeros(TargetPhotonNum, 1) ;

Var = sum(IRF_forSimu) ;
IRF_forSimu = IRF_forSimu / Var ; % area normalize the irf 
Imax = length(IRF_forSimu); % length of the time window (default is 12.5 ns)
IRF_probability = zeros(Imax,1) ;
Var = 0 ;
I = 0 ;
while I < Imax
    I = I + 1 ;
    Var = Var + IRF_forSimu(I) ;
    IRF_probability(I) = Var ;
end

%% Population of each state in equilibrium state
Var = ones(NumOfState) - eye(NumOfState) ;
Var = Var .* TransRate ;    % avoid diffution/photobleach effects
[EigenValues, ExpMatrix, EigenVector, InvEigenVector] = TK_RateEq_MakeExpMatrix(Var) ;

Var = ones(NumOfState, 1) ;
InitialPop_V = Var ./ NumOfState ; % intial population
Var = TransRate + (TransRate==0)*max(max(TransRate)) ;
Var = min(min(Var)) ;
Var = ExpMatrix .^ (1/Var*(10^10)) ;
Var = EigenVector * Var * InvEigenVector ;
EquiribriumPop_V = Var * InitialPop_V ; % equilibrium population
% Var = ones(NumOfState, 1) ;
% InitialPop_V = Var ./ NumOfState ;
% 
% Var = mean(mean(TransRate)) + 10^-10 ;
% Var = ExpMatrix .^ (1/Var*(10^10)) ;
% Var = EigenVector * Var * InvEigenVector ;
% EquiribriumPop_V = Var * InitialPop_V ;

[EigenValues, ExpMatrix, EigenVector, InvEigenVector] = TK_RateEq_MakeExpMatrix(TransRate) ;

Var = ExpMatrix .^ Tstep ;
Var_Posi_Trans = EigenVector * Var * InvEigenVector ;
        
%% Initial state
Var1 = EquiribriumPop_V ./ sum(EquiribriumPop_V) ;
Var = rand ;

Thre = 0 ;
I = 0 ;
while I < NumOfState
    I = I + 1 ;
    Thre = Thre + Var1(I) ;
    if Var < Thre
        Current_State = I ;
        break
    end
end

%% Photon detection
TotalT = 0 ;
%DwellT = 0 ;
CheckT = 0 ;
%CheckTransT = 0 ;
IniPop = zeros(NumOfState, 1) ;
IniPop(Current_State) = 1 ;
CheckDiffT = 0 ;
ParticleNum = 1 ;
MissPhotonCount = 0 ;
PhotonCount = 0 ;
while PhotonCount < TargetPhotonNum && TotalT < TargetTotalT
    TotalT = TotalT + Tstep ;
    
    if StateTransCheckEveryT_On1orOff0 == 1
        % State transition check
        Posi_Trans = Var_Posi_Trans * IniPop ;
        
        Ns = rand ;
        Thre = 0 ;
        LostCheck = 0 ;
        I = 0 ;
        while I < NumOfState
            I = I + 1 ;
%            if I ~= Current_State
                Thre = Thre + Posi_Trans(I) ;
                if Ns <= Thre
                    Current_State = I ;
                    %CheckTransT = 0 ;
                    IniPop(:) = 0 ;
                    IniPop(Current_State) = 1 ;
                    break
                end
%            end
            if I == NumOfState
                LostCheck = 1 ;
            end
        end
        
        if LostCheck == 1 % photobleach or diffusion occurs
            TotalT = TotalT + TafterPhotobleach ;
            Var1 = EquiribriumPop_V ./ sum(EquiribriumPop_V) ;
            Var = rand ;
                        
            Thre = 0 ;
            I = 0 ;
            while I < NumOfState
                I = I + 1 ;
                Thre = Thre + Var1(I) ;
                if Var < Thre
                    Current_State = I ;
                    IniPop(:) = 0 ;
                    IniPop(Current_State) = 1 ;
                    ParticleNum = ParticleNum + 1 ;
                    break
                end
            end
        end
    end
    
    Posi = Int(Current_State) * Tstep ;
    Np = rand ;
    if Np <= Posi
        PhotonCount = PhotonCount + 1 ;
        
        if StateTransCheckEveryT_On1orOff0 == 0
            DwellT = TotalT - CheckT ;
            CheckT = TotalT ;
            
            % state transition
            CheckTrans = zeros(NumOfState, 1) ;
            I = 0 ;
            while I < NumOfState
                I = I + 1 ;
                if I ~= Current_State
                    RateK = TransRate(Current_State, I) ;
                    Posi_Trans = 1-exp(-1*RateK*DwellT) ;
                    Ns = rand ;
                    if Ns <= Posi_Trans
                        CheckTrans(I) = 1 ;
                    end
                end
            end
            VarSum = sum(CheckTrans) ;
            if VarSum ~= 0
                Var = TransRate(Current_State, :)' ;
                CheckTransRate = Var .* CheckTrans ;
                Var1 = CheckTransRate ./ sum(CheckTransRate) ;
                Ns = rand ;
                
                Thre = 0 ;
                I = 0 ;
                while I < NumOfState
                    I = I + 1 ;
                    Thre = Thre + Var1(I) ;
                    if Ns < Thre
                        Current_State = I ;
                        break
                    end
                end
            end
        end
        
        Stream_State(PhotonCount) = Current_State ;
        Stream_detectedT(PhotonCount) = TotalT ;
        
        % lifetime
        Nt = rand ;
        t = -1*Life(Current_State)*log(Nt) ;    % ns
        Nt = rand ;
        Var = IRF_probability(:) - Nt ;
        Var1 = Var <= 0 ; 
        VarIndex = sum(Var1) + 1 ;
        t_IRFoffset = IRFxdata_forSimu(VarIndex) ;  % ns
        Stream_Life(PhotonCount) = round((t+t_IRFoffset)/tstep) ;
    else
        MissPhotonCount = MissPhotonCount + 1 ;
    end
       
    % diffusion check
    CheckDiffT = CheckDiffT + Tstep ;
    Posi_Diff = 1-exp(-1*RateDiff*CheckDiffT) ;
    Ns = rand ;
    if Ns < Posi_Diff
        Var1 = EquiribriumPop_V ./ sum(EquiribriumPop_V) ;
        Var = rand ;
        
        Thre = 0 ;
        I = 0 ;
        while I < NumOfState
            I = I + 1 ;
            Thre = Thre + Var1(I) ;
            if Var < Thre
                Current_State = I ;
                CheckDiffT = 0 ;
                ParticleNum = ParticleNum + 1 ;
                break
            end
        end
    end
    
end

if PhotonCount >= TargetPhotonNum
    EndCheck = 'Reach_TargetPhotonNum' ;
elseif TotalT >= TargetTotalT
    EndCheck = 'Reach_TargetTotalT' ;
end

kin1 = Stream_Life(1:PhotonCount) ;
tt1 = Stream_detectedT(1:PhotonCount) ;

clear CheckT CheckTrans CheckTransRate Current_Int Current_Life Current_State DwellT I 
clear Np Ns Nt Posi Posi_Trans RateK t t_IRFoffset Thre Var Var1 VarSum IRF_probability
clear EigenValues ExpMatrix EigenVector InvEigenVector InitialPop_V Posi_Diff CheckDiffT VarIndex
clear Stream_detectedT Stream_Life Stream_State

%% Check
if DoCheck1orNot0 == 1
    % Lifetime check
    Kmax = max(kin1) ;
    kin1_EachState = zeros(Kmax, NumOfState) ;
    Imax = length(kin1) ;
    I = 0 ;
    while I < Imax
        I = I + 1 ;
        K = kin1(I) ;
        S = Stream_State(I) ;
        if kin1(I) ~= 0
            kin1_EachState(K, S) = kin1_EachState(K, S) + 1 ;
        end
    end
    
    % Intensity check
    Imax = length(kin1) ;
    Int_EachState = zeros(NumOfState, Imax) ;
    Count_EachState = zeros(NumOfState, 1) ;
    IniState = Stream_State(1) ;
    DwellT = 0 ;
    PhotonCount = 0 ;
    I = 0 ;
    while I < Imax-1
        I = I + 1 ;
        DwellT = DwellT + (tt1(I+1)-tt1(I)) ;
        PhotonCount = PhotonCount + 1 ;
        if IniState ~= Stream_State(I+1)
            Count_EachState(IniState) = Count_EachState(IniState) + 1 ;
            Int_EachState(IniState, Count_EachState(IniState)) = PhotonCount / DwellT ;
            
            DwellT = 0 ;
            PhotonCount = 0 ;
            IniState = Stream_State(I+1) ;
        end
    end
    IntAve_EachState = zeros(NumOfState, 1) ;
    I = 0 ;
    while I < NumOfState
        I = I + 1 ;
        Var = Int_EachState(I, 1:Count_EachState(I)) ;
        IntAve_EachState(I) = sum(Var)/Count_EachState(I) ;
    end
    
    % TransRate check
    Imax = length(kin1) ;
    DwellHist_ForEachTrans = zeros(NumOfState, NumOfState, Imax) ;
    IniState = Stream_State(1) ;
    DwellT = 0 ;
    DwellT_max = 0 ;
    I = 0 ;
    while I < Imax-1
        I = I + 1 ;
        DwellT = DwellT + (tt1(I+1)-tt1(I)) ;
        if IniState ~= Stream_State(I+1)
            n = IniState ;
            m = Stream_State(I+1) ;
            Var = round(DwellT/Dwell_Tstep) ;
            %[n, m, Var]
            if Var ~= 0
                DwellHist_ForEachTrans(n, m, Var) = DwellHist_ForEachTrans(n, m, Var) + 1 ;
            end
            if DwellT_max < DwellT
                DwellT_max = DwellT ;
            end
            DwellT = 0 ;
            IniState = Stream_State(I+1) ;
        end
    end
    
    Var = round(DwellT_max/Dwell_Tstep) ;
    % Dwell histogram
    DwellHist_ForEachTrans = DwellHist_ForEachTrans(:, :, 1:Var) ;
    
    A = zeros(Var, 1) ;
    I = 0 ;
    while I < Var
        I = I + 1 ;
        A(I) = DwellHist_ForEachTrans(2, 2, I) ;
    end
    
    
    clear Count_EachState Dwell_ForEachTrans DwellT DwellT_max I Imax IniState Int_EachState K Kmax m n PhotonCount S Var
end

toc






