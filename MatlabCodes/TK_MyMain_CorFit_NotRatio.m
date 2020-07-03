% TK_MyMain_CorFit_NotRatio
%%
TargetM = Gresult_Mat_M_2DFLC_dTfun ;
TargetDecay = Gresult_Mat_M_Model_lin_dTfun ;

DisplayAll0orMap1orCor2 = 2 ;
Dis_2D_min = 1 ;

DataLineWidth = 1.5 ;
FontSize = 18 ;
LineWidth = 1.5 ;

%% Normalization with total photon number in each state
Var = size(TargetM) ;
Tmax = Var(3) ;
Imax = Var(1) ;
TrialNumFor_RegulatorConst
Var_Max = max(max(max(TargetM(1:Imax,1:Imax,:))))  /3 ;

Mat_PhotonSum = zeros(NumOfState) ;
N = 0 ;
while N < NumOfState
    N = N + 1 ;
    Var = Gresult_Mat_A_Global(:,N) .* Tau ;
    Var = sum(Var) ;
    Mat_PhotonSum(N,N) = Var ;
end
Var = Gresult_Mat_A_Global * inv(Mat_PhotonSum) ;

%%

if DisplayAll0orMap1orCor2 == 0 || DisplayAll0orMap1orCor2 == 2
    
    Ana_CorAuto = zeros(length(dT),NumOfState) ;
    Var_Num = NumOfState*(NumOfState-1)/2 + NumOfState ;
    Ana_CorCross = zeros(length(dT),Var_Num) ;
    
    StrF_Legend_CorAuto = zeros(NumOfState, 1) ;
    StrF_Legend_CorCross = zeros(Var_Num, 3) ;
    TK_MyMain_Analyze_03_NotRatio
    
    Trial_I = 0 ;
    while Trial_I < TrialNumber
        Trial_I = Trial_I + 1 
        TK_MyMain_Analyze_03_NotRatio
        [Ana_estimates, Ana_Free0orFix1orAbs2, Ana_estimates_A, Ana_estimates_y0, ...
            Ana_estimates_E, Ana_estimates_Q, Ana_estimates_RateM, Ana_FittedX, Ana_FittedY, model] = ...
            TK_FitF_CorrelationDecay_RateMat_16_NotRatio(...
            dT, Ana_Correlation, ComponentNumber, StateAssign, FitA, Fity0, FitE, FitQ, FitRateM,...
            Free0orFix1orAbs2_A, Free0orFix1orAbs2_y0, Free0orFix1orAbs2_E, Free0orFix1orAbs2_Q, Free0orFix1orAbs2_RateM,...
            FittedPoints, FitNormalizedCor_y1n0, FitWeightFactor) ;
    end
    
    %%%%%% population %%%%%%%%
    ComNumMax = ComponentNumber ;
    Ana_EquPop_V_Substate = zeros(length(StateAssign), ComNumMax) ;
    Ana_EquPop_V = zeros(max(StateAssign), ComNumMax) ;
    ComNum = 0 ;
    while ComNum < ComNumMax
        ComNum = ComNum + 1 ;
        
        Var = Ana_estimates_RateM(:,:,ComNum) - diag(diag(Ana_estimates_RateM(:,:,ComNum))) ;   % remove photobleach component
        [EigenValues, ExpMatrix, EigenVector, InvEigenVector] = TK_RateEq_MakeExpMatrix(Var) ;
        
        Var = length(StateAssign) ;
        InitialPop_V = ones(Var, 1) ./ Var ;
        
        Var = ExpMatrix .^ (10^10) ;
        Var = EigenVector * Var * InvEigenVector ;
        Ana_EquPop_V_Substate(:,ComNum) = Var * InitialPop_V ;
        %Ana_EquPop_V = zeros(max(StateAssign),1) ;
        I = 0 ;
        while I < length(StateAssign)
            I = I + 1 ;
            Ana_EquPop_V(StateAssign(I),ComNum) = Ana_EquPop_V(StateAssign(I),ComNum) + Ana_EquPop_V_Substate(I,ComNum) ;
        end
        Var = sum(Ana_EquPop_V(:,ComNum)) ;
        Ana_EquPop_V(:,ComNum) = Ana_EquPop_V(:,ComNum) ./ Var ;
        Var = sum(Ana_EquPop_V_Substate(:,ComNum)) ;
        Ana_EquPop_V_Substate(:,ComNum) = Ana_EquPop_V_Substate(:,ComNum) ./ Var ;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%
    
    Count = 0 ;
    VarN1 = 0 ;
    while VarN1 < NumOfState
        VarN1 = VarN1 + 1 ;
        
        VarN2 = VarN1-1 ;
        while VarN2 < NumOfState
            VarN2 = VarN2 + 1 ;
            
            %%
            Var = size(Ana_Correlation) ;
            Var_Cor1 = zeros(Var(3),1) ;
            Var_Cor1(:) = Ana_Correlation(VarN1, VarN2, :) ;
            Var = size(Ana_FittedY) ;
            Var_Cor2 = zeros(Var(3),1) ;
            Var_Cor2(:) = Ana_FittedY(VarN1, VarN2, :) ;

            figure;
            semilogx(dT, Var_Cor1(:,1),'o',...
                Ana_FittedX, Var_Cor2(:,1),'r',...
                'LineWidth',DataLineWidth, 'MarkerSize', 10) ;
            StrF = strcat('CorRatio - ', num2str(VarN1), num2str(VarN2), ' (F1)') ;
            title(StrF)
            xlabel('dT (s)','FontSize',FontSize,'FontWeight','bold','Color','k')
            ylabel('Correlation ratio','FontSize',FontSize,'FontWeight','bold','Color','k')
            set(gca, 'tickdir', 'out') ; set(gca,'LineWidth',LineWidth); set(gca,'FontSize',FontSize); set(gca, 'ticklength', [0.01 0.01]) ;
 
            Count = Count + 1 ;
            Ana_CorCross(:,Count) = Ana_Correlation(VarN1, VarN2, :) ;
            StrF_Legend_CorCross(Count,:) = strcat(num2str(VarN1), '-', num2str(VarN2)) ;
        end
        
        Ana_CorAuto(:,VarN1) = Ana_Correlation(VarN1, VarN1, :) ;
        StrF_Legend_CorAuto(VarN1,:) = num2str(VarN1) ;
    end
end


if DisplayAll0orMap1orCor2 == 0 || DisplayAll0orMap1orCor2 == 1
    %%
    figure;
    semilogx(dT, Ana_CorAuto,'-o','LineWidth',DataLineWidth, 'MarkerSize', 10) ;
    title('Auto correation')
    %xlim([min(dT) max(dT)]) ;
    xlabel('dT (s)','FontSize',FontSize,'FontWeight','bold','Color','k')
    ylabel('Auto correlation','FontSize',FontSize,'FontWeight','bold','Color','k')
    set(gca, 'tickdir', 'out') ; set(gca,'LineWidth',LineWidth); set(gca,'FontSize',FontSize); set(gca, 'ticklength', [0.01 0.01]) ;
    legend(char(StrF_Legend_CorAuto)) ;
    
    StrF = strcat('Cor_Auto','.pdf') ; print('-r300','-dpdf',StrF);
    StrF = strcat('Cor_Auto') ; saveas(gcf, StrF, 'png') ;
    %close gcf ;
    %%
    figure;
    semilogx(dT, Ana_CorCross,'-o','LineWidth',DataLineWidth, 'MarkerSize', 10) ;
    title('Cross correation')
    %xlim([min(dT) max(dT)]) ;
    xlabel('dT (s)','FontSize',FontSize,'FontWeight','bold','Color','k')
    ylabel('Cross correlation','FontSize',FontSize,'FontWeight','bold','Color','k')
    set(gca, 'tickdir', 'out') ; set(gca,'LineWidth',LineWidth); set(gca,'FontSize',FontSize); set(gca, 'ticklength', [0.01 0.01]) ;
    legend(char(StrF_Legend_CorCross)) ;
    
    StrF = strcat('Cor_Cross','.pdf') ; print('-r300','-dpdf',StrF);
    StrF = strcat('Cor_Cross') ; saveas(gcf, StrF, 'png') ;
    %close gcf ;
    %%
    Var = size(Ana_CorRatio) ;
    I = 0 ;
    while I < Var(3)
        I = I + 1 ;
        figure;
        Var_A = Ana_CorRatio(:,:,I) ;
        semilogx(dT, Var_A,'-o','LineWidth',DataLineWidth, 'MarkerSize', 10) ;
        if I == 1
            title(strcat('Correation ratio',' - g21/g22'))
        elseif I == 2
            title(strcat('Correation ratio',' - (g12+g21) /(g11+g22)'))
        elseif I == 3
            title(strcat('Correation ratio',' - (g12xg21) /(g11xg22)'))
        end
        %xlim([min(dT) max(dT)]) ;
        xlabel('dT (s)','FontSize',FontSize,'FontWeight','bold','Color','k')
        ylabel('Correlation ratio','FontSize',FontSize,'FontWeight','bold','Color','k')
        set(gca, 'tickdir', 'out') ; set(gca,'LineWidth',LineWidth); set(gca,'FontSize',FontSize); set(gca, 'ticklength', [0.01 0.01]) ;
        legend(char(StrF_Legend_CorCross)) ;
        
        StrF = strcat('Cor_Ratio_',num2str(I),'.pdf') ; print('-r300','-dpdf',StrF);
        StrF = strcat('Cor_Ratio_',num2str(I)) ; saveas(gcf, StrF, 'png') ;
        %close gcf ;
    end

end

clear TargetM VarN1 VarN2 VarN StrF Count T Var Var_Cor StrF_Legend_CorAuto StrF_Legend_CorCross
