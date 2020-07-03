% TK_MyMain_Fig
%%
TargetM = Gresult_Mat_M_2DFLC_dTfun ;
TargetDecay = Gresult_Mat_M_Model_lin_dTfun ;

DisplayAll0orMap1orCor2 = 1 ;

DisplayAllFigure = 5 ;
%DisplayAllFigure = [1:length(dT)] ;

Dis_2D_min = 1 ;
Dis_2D_max = 80 ;

Dis_maxFactor = 70 ;

DataLineWidth = 1.5 ;
FontSize = 18 ;
LineWidth = 1.5 ;

if DisplayAll0orMap1orCor2 == 0 || DisplayAll0orMap1orCor2 == 2
    Var = size(size(Ana_FittedY)) ;
    if Var ~= 4
        Ana_FittedY_Hozzon = Ana_FittedY ;
    end
    Ana_FittedY = zeros(NumOfState,NumOfState,length(Ana_FittedX),4) ;
    Ana_FittedY(:,:,:,1) = Ana_FittedY_Hozzon ;
    Ana_FittedY(:,:,:,2) = Ana_FittedY_Hozzon ;
    Ana_FittedY(:,:,:,3) = Ana_FittedY_Hozzon ;
    Ana_FittedY(:,:,:,4) = Ana_FittedY_Hozzon ;
    
    % Ana_FittedY = zeros(NumOfState,NumOfState,length(Ana_FittedX),4) ;
    % Ana_FittedY(:,:,:,1) = Ana_C2_FittedY ;
    % Ana_FittedY(:,:,:,2) = Ana_I2_FittedY ;
    % Ana_FittedY(:,:,:,3) = Ana_C3_FittedY ;
    % Ana_FittedY(:,:,:,4) = Ana_I3_FittedY ;
    
    % Ana_FittedY(:,:,:,1) = Ana_C2_FittedY ;
    % Ana_FittedY(:,:,:,2) = Ana_I2_FittedY ;
    % Ana_FittedY(:,:,:,3) = Ana_1_FittedY ;
    % Ana_FittedY(:,:,:,4) = Ana_1_FittedY ;
end


%%
Var = size(TargetM) ;
Tmax = Var(3) ;
%Imax = Dis_2D_max ;%Var(1) ;

Var_Max = max(max(max(TargetM(Dis_2D_min:Dis_2D_max,Dis_2D_min:Dis_2D_max,:))))  /1.5 ;
% Var = size(TargetM(Dis_2D_min:Dis_2D_max,Dis_2D_min:Dis_2D_max,:)) ;
% Var = Var(1) * Var(2) * Var(3) ;
% Var_Max = sum(sum(sum(TargetM(Dis_2D_min:Dis_2D_max,Dis_2D_min:Dis_2D_max,:))))/ Var  *Dis_maxFactor ;

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
if DisplayAll0orMap1orCor2 == 0 || DisplayAll0orMap1orCor2 == 1
    figure;
    plot(Tau(Dis_2D_min:Dis_2D_max),Var(Dis_2D_min:Dis_2D_max,:),'LineWidth',DataLineWidth) ;
    xlabel('Tau (ns)','FontSize',FontSize,'FontWeight','bold','Color','k')
    ylabel('Amplitude','FontSize',FontSize,'FontWeight','bold','Color','k')
    set(gca, 'tickdir', 'out') ; set(gca,'LineWidth',LineWidth); set(gca,'FontSize',FontSize); set(gca, 'ticklength', [0.01 0.01]) ;
    StrF = strcat('A_Gresult_N.pdf') ; print('-r300','-dpdf',StrF);
    saveas(gcf, 'A_Gresult_N', 'png') ;
    close gcf ;
    
    Var = max(Gresult_Mat_A_Global) ;
    Var1 = size(Gresult_Mat_A_Global) ;
    Var2 = Gresult_Mat_A_Global ;
    I = 0 ;
    while I < Var1(2)
        I = I + 1 ;
        Var2(:,I) = Gresult_Mat_A_Global(:,I) / Var(I) ;
    end
    figure;
    plot(Tau(Dis_2D_min:Dis_2D_max),Var2(Dis_2D_min:Dis_2D_max,:),'LineWidth',DataLineWidth) ;
    xlabel('Tau (ns)','FontSize',FontSize,'FontWeight','bold','Color','k')
    ylabel('Amplitude','FontSize',FontSize,'FontWeight','bold','Color','k')
    set(gca, 'tickdir', 'out') ; set(gca,'LineWidth',LineWidth); set(gca,'FontSize',FontSize); set(gca, 'ticklength', [0.01 0.01]) ;
    StrF = strcat('A_Gresult_Namp.pdf') ; print('-r300','-dpdf',StrF);
    saveas(gcf, 'A_Gresult_Namp', 'png') ;
    close gcf ;
%     plot(Tau,Result_Short_Mat_A,'LineWidth',DataLineWidth) ;
%     xlabel('Tau (ns)','FontSize',FontSize,'FontWeight','bold','Color','k')
%     ylabel('Amplitude','FontSize',FontSize,'FontWeight','bold','Color','k')
%     set(gca, 'tickdir', 'out') ; set(gca,'LineWidth',LineWidth); set(gca,'FontSize',FontSize); set(gca, 'ticklength', [0.01 0.01]) ;
%     StrF = strcat('A_Result.pdf') ; print('-r300','-dpdf',StrF);
%     saveas(gcf, 'A_Result', 'png') ;
%     close gcf ;
    
    figure;
    plot(Tau(Dis_2D_min:Dis_2D_max),Gresult_Mat_A_Global(Dis_2D_min:Dis_2D_max,:),'LineWidth',DataLineWidth) ;
    xlabel('Tau (ns)','FontSize',FontSize,'FontWeight','bold','Color','k')
    ylabel('Amplitude','FontSize',FontSize,'FontWeight','bold','Color','k')
    set(gca, 'tickdir', 'out') ; set(gca,'LineWidth',LineWidth); set(gca,'FontSize',FontSize); set(gca, 'ticklength', [0.01 0.01]) ;
    StrF = strcat('A_Gresult.pdf') ; print('-r300','-dpdf',StrF);
    saveas(gcf, 'A_Gresult', 'png') ;
    close gcf ;
    
    T = 0 ;
    while T < Tmax
        T = T + 1 ;
        Var1 = Tau(Dis_2D_min:Dis_2D_max) ;
        Var2 = TargetM(Dis_2D_min:Dis_2D_max,Dis_2D_min:Dis_2D_max,T) ;
        %TK_DisIntLife2Dmap(Tau(Dis_2D_min:Dis_2D_max), Tau(Dis_2D_min:Dis_2D_max), TargetM(Dis_2D_min:Dis_2D_max,Dis_2D_min:Dis_2D_max,T),1,0) ;
        TK_DisIntLife2Dmap(Var1, Var1, Var2,1,0) ;
        
        %g2Dmax = max(max(Var2)) ;
        g2Dmax = Var_Max ;
        xlim([0 max(Var1)]) ; ylim([0 max(Var1)]) ; caxis([0 Var_Max]) ;
        xlabel('Tau (ns)','FontSize',FontSize,'FontWeight','bold','Color','k')
        ylabel('Tau (ns)','FontSize',FontSize,'FontWeight','bold','Color','k')
        pbaspect([1 1 1]) ;
        set(gca, 'tickdir', 'out') ; set(gca,'LineWidth',LineWidth); set(gca,'FontSize',FontSize); set(gca, 'ticklength', [0.01 0.01]) ;
        LabelHandle = get(colorbar, 'Label') ; set(LabelHandle,'Rotation',-90) ; set(LabelHandle,'String','Total dwell time /s') ; set(LabelHandle,'Position', [4 g2Dmax/2 0]) ;
        
        StrF = strcat('B',num2str(T),'.pdf') ; print('-r300','-dpdf',StrF);
        StrF = strcat('B',num2str(T)) ; saveas(gcf, StrF, 'png') ;
        close gcf ;
        
        if sum(T == DisplayAllFigure)
            Var_minmax = zeros(2,1) ;
            Var_minmax(1) =  max(max(max(Mat_2DFDC_lin(:,:,:)))) ; %max(max(Mat_2DFDC_lin(:,:,T))) ;   
            Var_minmax(2) =  max(max(max(Gresult_Mat_M_Model_lin_dTfun(:,:,:)))) ; %max(max(Gresult_Mat_M_Model_lin_dTfun(:,:,T))) ;    
            g2Dmax = max(Var_minmax) ;  % max value of 2D map
            Var_minmax(1) =  min(min(min(Mat_2DFDC_lin(:,:,:)))) ; %min(min(Mat_2DFDC_lin(:,:,T))) ;   
            Var_minmax(2) =  min(min(min(Gresult_Mat_M_Model_lin_dTfun(:,:,:)))) ; %min(min(Gresult_Mat_M_Model_lin_dTfun(:,:,T))) ;    
            g2Dmin = min(Var_minmax) ;  % max value of 2D map
            
            TK_DisIntLife2Dmap(Gresult_Mat_2DFDC_lin_t, Gresult_Mat_2DFDC_lin_t, Mat_2DFDC_lin(:,:,T),1,0) ;
            caxis([g2Dmin g2Dmax])
            %xlim([0 max(Gresult_Mat_2DFDC_lin_t)]) ; ylim([0 max(Gresult_Mat_2DFDC_lin_t)]) ;
            xlabel('t (ns)','FontSize',FontSize,'FontWeight','bold','Color','k')
            ylabel('t (ns)','FontSize',FontSize,'FontWeight','bold','Color','k')
            pbaspect([1 1 1]) ;
            set(gca, 'tickdir', 'out') ; set(gca,'LineWidth',LineWidth); set(gca,'FontSize',FontSize); set(gca, 'ticklength', [0.01 0.01]) ;
            LabelHandle = get(colorbar, 'Label') ; set(LabelHandle,'Rotation',-90) ; set(LabelHandle,'String','Total dwell time /s') ; set(LabelHandle,'Position', [4 g2Dmax/2 0]) ;
            StrF = strcat('B',num2str(T),'_2DFLD',num2str(T),'.pdf') ; print('-r300','-dpdf',StrF);
            StrF = strcat('B',num2str(T),'_2DFLD',num2str(T)) ; saveas(gcf, StrF, 'png') ;
            close gcf ;
            
            TK_DisIntLife2Dmap(Gresult_Mat_2DFDC_lin_t, Gresult_Mat_2DFDC_lin_t, Gresult_Mat_M_Model_lin_dTfun(:,:,T),1,0) ;
            caxis([g2Dmin g2Dmax])
            %xlim([0 max(Gresult_Mat_2DFDC_lin_t)]) ; ylim([0 max(Gresult_Mat_2DFDC_lin_t)]) ;
            xlabel('t (ns)','FontSize',FontSize,'FontWeight','bold','Color','k')
            ylabel('t (ns)','FontSize',FontSize,'FontWeight','bold','Color','k')
            pbaspect([1 1 1]) ;
            set(gca, 'tickdir', 'out') ; set(gca,'LineWidth',LineWidth); set(gca,'FontSize',FontSize); set(gca, 'ticklength', [0.01 0.01]) ;
            LabelHandle = get(colorbar, 'Label') ; set(LabelHandle,'Rotation',-90) ; set(LabelHandle,'String','Total dwell time /s') ; set(LabelHandle,'Position', [4 g2Dmax/2 0]) ;
            
            StrF = strcat('B',num2str(T),'_2DFLD_Fit',num2str(T),'.pdf') ; print('-r300','-dpdf',StrF);
            StrF = strcat('B',num2str(T),'_2DFLD_Fit',num2str(T)) ; saveas(gcf, StrF, 'png') ;
            close gcf ;
        end
        
        [Plot_Var1, Plot_Var2] = max(max(TargetDecay(:,:,T))) ;
        Plot_X = Gresult_Mat_2DFDC_lin_t ;
        if UseCor1orNot0 == 0
            Plot_Y1 = Mat_2DFDC_lin(Plot_Var2,:,T) ;
        elseif UseCor1orNot0 == 1
            Plot_Y1 = Mat_2DFDC_cor_lin(Plot_Var2,:,T) ;
        end
        Plot_Y2 = TargetDecay(Plot_Var2,:,T) ;
        figure; Plot_hp2 = plot(Plot_X,Plot_Y1,Plot_X,Plot_Y2) ;
        Plot_hp2(2).LineWidth = DataLineWidth*3 ;
        set(gca, 'tickdir', 'out') ; set(gca,'LineWidth',LineWidth); set(gca,'FontSize',FontSize); set(gca, 'ticklength', [0.01 0.01]) ;
        StrF = strcat('C',num2str(T),'.pdf') ; print('-r300','-dpdf',StrF);
        StrF = strcat('C',num2str(T)) ; saveas(gcf, StrF, 'png') ;
        close gcf ;
    end
end

if DisplayAll0orMap1orCor2 == 0 || DisplayAll0orMap1orCor2 == 2
    
    Ana_CorAuto = zeros(length(dT),NumOfState) ;
    Var_Num = NumOfState*(NumOfState-1)/2 ;
    Ana_CorCross = zeros(length(dT),Var_Num) ;
    Var_Num2 = NumOfState*(NumOfState+1)/2 ;
    Ana_CorAll = zeros(length(dT),Var_Num2) ;
    
    StrF_Legend_CorAuto = zeros(NumOfState, 1) ;
    StrF_Legend_CorCross = zeros(Var_Num, 3) ;
    StrF_Legend_CorAll = zeros(Var_Num2, 3) ;
    
    Count = 0 ;
    Count2 = 0 ;
    VarN1 = 0 ;
    while VarN1 < NumOfState
        VarN1 = VarN1 + 1 ;
        
        VarN2 = VarN1-1 ;
        while VarN2 < NumOfState
            VarN2 = VarN2 + 1 ;
            
                
            %%       
            Var = length(dT) ;
            Var_Cor1 = zeros(Var,1) ;
            Var_Cor1(:) = Ana_Correlation(VarN1, VarN2, :) ;
            
            Var = length(Ana_FittedX) ;
            Var_Cor2 = zeros(Var, 2) ;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            Var_Cor2(:, 1) = Ana_FittedY(VarN1,VarN2,:,1) ;
            Var_Cor2(:, 2) = Ana_FittedY(VarN1,VarN2,:,2) ;
            
            figure;
            semilogx(dT, Var_Cor1(:,1),'o',...
                Ana_FittedX, Var_Cor2(:,1),'r',...
                Ana_FittedX, Var_Cor2(:,2),'b',...
                'LineWidth',DataLineWidth, 'MarkerSize', 10) ;
            StrF = strcat('CorRatio - ', num2str(VarN1), num2str(VarN2), ' (F1)') ;
            title(StrF)
            xlabel('dT (s)','FontSize',FontSize,'FontWeight','bold','Color','k')
            ylabel('Correlation ratio','FontSize',FontSize,'FontWeight','bold','Color','k')
            set(gca, 'tickdir', 'out') ; set(gca,'LineWidth',LineWidth); set(gca,'FontSize',FontSize); set(gca, 'ticklength', [0.01 0.01]) ;
            legend('g21/g22','Com','Ind','Location','south') ;
            StrF = strcat('Ca1_',num2str(VarN1),num2str(VarN2),'.pdf') ; print('-r300','-dpdf',StrF);
            StrF = strcat('Ca1_',num2str(VarN1),num2str(VarN2)) ; saveas(gcf, StrF, 'png') ;
            close gcf ;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            Var_Cor2(:, 1) = Ana_FittedY(VarN1,VarN2,:,3) ;
            Var_Cor2(:, 2) = Ana_FittedY(VarN1,VarN2,:,4) ;
            
            figure;
            semilogx(dT, Var_Cor1(:,1),'o',...
                Ana_FittedX, Var_Cor2(:,1),'r',...
                Ana_FittedX, Var_Cor2(:,2),'b',...
                'LineWidth',DataLineWidth, 'MarkerSize', 10) ;
            StrF = strcat('CorRatio - ', num2str(VarN1), num2str(VarN2), ' (F1)') ;
            title(StrF)
            xlabel('dT (s)','FontSize',FontSize,'FontWeight','bold','Color','k')
            ylabel('Correlation ratio','FontSize',FontSize,'FontWeight','bold','Color','k')
            set(gca, 'tickdir', 'out') ; set(gca,'LineWidth',LineWidth); set(gca,'FontSize',FontSize); set(gca, 'ticklength', [0.01 0.01]) ;
            legend('g21/g22','Com','Ind','Location','south') ;
            StrF = strcat('Cb1_',num2str(VarN1),num2str(VarN2),'.pdf') ; print('-r300','-dpdf',StrF);
            StrF = strcat('Cb1_',num2str(VarN1),num2str(VarN2)) ; saveas(gcf, StrF, 'png') ;
            close gcf ;
            
            if VarN1 ~= VarN2
                Count = Count + 1 ;
                Ana_CorCross(:,Count) = Ana_Correlation(VarN1, VarN2, :) ;
                StrF_Legend_CorCross(Count,:) = strcat(num2str(VarN1), '-', num2str(VarN2)) ;
            end
            
            Count2 = Count2 + 1 ;
            Ana_CorAll(:,Count2) = Ana_Correlation(VarN1, VarN2, :) ;
            StrF_Legend_CorAll(Count2,:) = strcat(num2str(VarN1), '-', num2str(VarN2)) ;
        end
        
        Ana_CorAuto(:,VarN1) = Ana_Correlation(VarN1, VarN1, :) ;
        StrF_Legend_CorAuto(VarN1,:) = num2str(VarN1) ;
    end
end


if DisplayAll0orMap1orCor2 == 0 || DisplayAll0orMap1orCor2 == 2 
    %%
    figure;
    semilogx(dT, Ana_CorAuto,'-o','LineWidth',DataLineWidth, 'MarkerSize', 10) ;
    title('Auto correation')
    xlim([min(dT)*0.7 max(dT)*1.3]) ;
    Var = max(max(Ana_CorAuto)) ;
    ylim([-0.1*Var, Var*1.1]) ;
    xlabel('dT (s)','FontSize',FontSize,'FontWeight','bold','Color','k')
    ylabel('Auto correlation','FontSize',FontSize,'FontWeight','bold','Color','k')
    set(gca, 'tickdir', 'out') ; set(gca,'LineWidth',LineWidth); set(gca,'FontSize',FontSize); set(gca, 'ticklength', [0.01 0.01]) ;
    legend(char(StrF_Legend_CorAuto),'Location','east') ;
    pbaspect([8*1.3 10 1])
%     xticks([10^-4 10^-3 10^-2 10^-1 10^0 10^1])
%     xticklabels({'-4','-3','-2','-1','0','1'})
    
    StrF = strcat('Cor_Auto','.pdf') ; print('-r1200','-dpdf',StrF);
    StrF = strcat('Cor_Auto') ; saveas(gcf, StrF, 'png') ;
    %close gcf ;
    %%
    figure;
    semilogx(dT, Ana_CorCross,'-o','LineWidth',DataLineWidth, 'MarkerSize', 10) ;
    title('Cross correation')
    xlim([min(dT)*0.7 max(dT)*1.3]) ;
    Var = max(max(Ana_CorCross)) ;
    ylim([-0.1*Var, Var*1.1]) ;
    xlabel('dT (s)','FontSize',FontSize,'FontWeight','bold','Color','k')
    ylabel('Cross correlation','FontSize',FontSize,'FontWeight','bold','Color','k')
    set(gca, 'tickdir', 'out') ; set(gca,'LineWidth',LineWidth); set(gca,'FontSize',FontSize); set(gca, 'ticklength', [0.01 0.01]) ;
    legend(char(StrF_Legend_CorCross),'Location','east') ;
    pbaspect([8*1.3 10 1])
    
    StrF = strcat('Cor_Cross','.pdf') ; print('-r1200','-dpdf',StrF);
    StrF = strcat('Cor_Cross') ; saveas(gcf, StrF, 'png') ;
    %close gcf ;
% end
% if DisplayAll0orMap1orCor2 == 3
    %%
    figure;
    semilogx(dT, Ana_CorAll,'-o','LineWidth',DataLineWidth, 'MarkerSize', 10) ;
    title('All correation')
    xlim([min(dT)*0.7 max(dT)*1.3]) ;
    Var = max(max(Ana_CorAll)) ;
    ylim([-0.1*Var, Var*1.1]) ;
    xlabel('dT (s)','FontSize',FontSize,'FontWeight','bold','Color','k')
    ylabel('Correlation','FontSize',FontSize,'FontWeight','bold','Color','k')
    set(gca, 'tickdir', 'out') ; set(gca,'LineWidth',LineWidth); set(gca,'FontSize',FontSize); set(gca, 'ticklength', [0.01 0.01]) ;
    legend(char(StrF_Legend_CorAll),'Location','northeast') ; %southeast  east
    pbaspect([8*1.3 10 1])
    
    StrF = strcat('Cor_All','.pdf') ; print('-r1200','-dpdf',StrF);
    StrF = strcat('Cor_All') ; saveas(gcf, StrF, 'png') ;

end


clear TargetM VarN1 VarN2 VarN StrF Count T Var Var_Cor StrF_Legend_CorAuto StrF_Legend_CorCross
