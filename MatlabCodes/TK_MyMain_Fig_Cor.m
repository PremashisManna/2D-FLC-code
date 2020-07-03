% TK_MyMain_Fig_Cor

%%

DataLineWidth = 1.5 ;
FontSize = 18 ;
LineWidth = 1.5 ;

%Ana_FittedY = Ana_I3_FittedY ;

%Ana_FittedY = zeros(NumOfState,NumOfState,length(Ana_FittedX),4) ;
% Ana_FittedY(:,:,:,1) = Ana_C2_FittedY ;
% Ana_FittedY(:,:,:,2) = Ana_I2_FittedY ;
% Ana_FittedY(:,:,:,3) = Ana_C3_FittedY ;
% Ana_FittedY(:,:,:,4) = Ana_I3_FittedY ;
% Ana_FittedY(:,:,:,1) = Ana_C2_FittedY ;
% Ana_FittedY(:,:,:,2) = Ana_I2_FittedY ;
% Ana_FittedY(:,:,:,3) = Ana_1_FittedY ;
% Ana_FittedY(:,:,:,4) = Ana_1_FittedY ;

    Ana_CorAuto = zeros(length(dT),NumOfState) ;
    Ana_CorAuto_Fit = zeros(length(Ana_FittedX),NumOfState) ;
    Var_Num = NumOfState*(NumOfState-1)/2 ;
    Ana_CorCross = zeros(length(dT),Var_Num) ;
    Ana_CorCross_Fit = zeros(length(Ana_FittedX),Var_Num) ;
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
            
%             Var = length(Ana_FittedX) ;
%             Var_Cor2 = zeros(Var, 2) ;
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%
%             Var_Cor2(:, 1) = Ana_FittedY(VarN1,VarN2,:,1) ;
%             Var_Cor2(:, 2) = Ana_FittedY(VarN1,VarN2,:,2) ;
%             
%             figure;
%             semilogx(dT, Var_Cor1(:,1),'o',...
%                 Ana_FittedX, Var_Cor2(:,1),'r',...
%                 Ana_FittedX, Var_Cor2(:,2),'b',...
%                 'LineWidth',DataLineWidth, 'MarkerSize', 10) ;
%             StrF = strcat('CorRatio - ', num2str(VarN1), num2str(VarN2), ' (F1)') ;
%             title(StrF)
%             xlabel('dT (s)','FontSize',FontSize,'FontWeight','bold','Color','k')
%             ylabel('Correlation ratio','FontSize',FontSize,'FontWeight','bold','Color','k')
%             set(gca, 'tickdir', 'out') ; set(gca,'LineWidth',LineWidth); set(gca,'FontSize',FontSize); set(gca, 'ticklength', [0.01 0.01]) ;
%             legend('g21/g22','Com','Ind','Location','south') ;
%             StrF = strcat('Ca1_',num2str(VarN1),num2str(VarN2),'.pdf') ; print('-r300','-dpdf',StrF);
%             StrF = strcat('Ca1_',num2str(VarN1),num2str(VarN2)) ; saveas(gcf, StrF, 'png') ;
%             close gcf ;
%             
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%
%             Var_Cor2(:, 1) = Ana_FittedY(VarN1,VarN2,:,3) ;
%             Var_Cor2(:, 2) = Ana_FittedY(VarN1,VarN2,:,4) ;
%             
%             figure;
%             semilogx(dT, Var_Cor1(:,1),'o',...
%                 Ana_FittedX, Var_Cor2(:,1),'r',...
%                 Ana_FittedX, Var_Cor2(:,2),'b',...
%                 'LineWidth',DataLineWidth, 'MarkerSize', 10) ;
%             StrF = strcat('CorRatio - ', num2str(VarN1), num2str(VarN2), ' (F1)') ;
%             title(StrF)
%             xlabel('dT (s)','FontSize',FontSize,'FontWeight','bold','Color','k')
%             ylabel('Correlation ratio','FontSize',FontSize,'FontWeight','bold','Color','k')
%             set(gca, 'tickdir', 'out') ; set(gca,'LineWidth',LineWidth); set(gca,'FontSize',FontSize); set(gca, 'ticklength', [0.01 0.01]) ;
%             legend('g21/g22','Com','Ind','Location','south') ;
%             StrF = strcat('Cb1_',num2str(VarN1),num2str(VarN2),'.pdf') ; print('-r300','-dpdf',StrF);
%             StrF = strcat('Cb1_',num2str(VarN1),num2str(VarN2)) ; saveas(gcf, StrF, 'png') ;
%             close gcf ;
            
            if VarN1 ~= VarN2
                Count = Count + 1 ;
                Ana_CorCross(:,Count) = Ana_Correlation(VarN1, VarN2, :) ;
                StrF_Legend_CorCross(Count,:) = strcat(num2str(VarN1), '-', num2str(VarN2)) ;
                Ana_CorCross_Fit(:,Count) = Ana_FittedY(VarN1,VarN2,:) ;
            end
            
            Count2 = Count2 + 1 ;
            Ana_CorAll(:,Count2) = Ana_Correlation(VarN1, VarN2, :) ;
            StrF_Legend_CorAll(Count2,:) = strcat(num2str(VarN1), '-', num2str(VarN2)) ;
        end
        
        Ana_CorAuto(:,VarN1) = Ana_Correlation(VarN1, VarN1, :) ;
        StrF_Legend_CorAuto(VarN1,:) = num2str(VarN1) ;
        Ana_CorAuto_Fit(:,VarN1) = Ana_FittedY(VarN1, VarN1, :) ;
    end


%%
Fig1 = figure;
Fig1.OuterPosition = [500, 500, 800, 400] ;
Var1 = Ana_CorAuto(:,1) ;
VarAna = Ana_CorAuto_Fit(:,1) ;
Var = vertcat(Var1, VarAna) ;
semilogx(dT, Var1,'ok',...
    Ana_FittedX, VarAna, 'r',...
    'LineWidth',DataLineWidth, 'MarkerSize', 10) ;

% VarAna = Ana_CorAuto_Fit(:,1) ;
% VarAna1 = VarAna;  VarAna1(:) = Ana_CorAuto_Fit(:,1) ;
% VarAna2 = VarAna;  VarAna2(:) = Ana_CorAuto_Fit_free(:,1) ;
% Var = vertcat(Var1, VarAna1, VarAna2) ;
% semilogx(dT, Var1,'ok',...  
%     Ana_FittedX, VarAna1, 'r',...
%     Ana_FittedX, VarAna2, 'b',...
%     'LineWidth',DataLineWidth, 'MarkerSize', 10) ;

% VarAna = Ana_CorAuto_Fit(:,1) ;
% VarAna1 = VarAna;  VarAna1(:) = Fit2stateFix_Ana_FittedY(1,1,:) ;
% VarAna2 = VarAna;  VarAna2(:) = Fit3stateFix_Ana_FittedY(1,1,:) ;
% VarAna3 = VarAna;  VarAna3(:) = Fit3state_Ana_FittedY(1,1,:) ;
% VarAna4 = VarAna;  VarAna4(:) = Fit3stateFixBG_Ana_FittedY(1,1,:) ;
% %VarAna5 = VarAna;  VarAna5(:) = Fit3stateFreeBG_Ana_FittedY(1,1,:) ;
% Var = vertcat(Var1, VarAna2, VarAna3, VarAna4) ;
% semilogx(dT, Var1,'ok',...
%     Ana_FittedX, VarAna1, 'b',...    
%     Ana_FittedX, VarAna2, 'r',...
%     Ana_FittedX, VarAna3, 'm',...
%     Ana_FittedX, VarAna4, 'g',...%Ana_FittedX, VarAna5, 'g',...
%     'LineWidth',DataLineWidth, 'MarkerSize', 10) ;

% Fig1 = figure;
% Fig1.OuterPosition = [500, 500, 800, 400] ;
% Var1 = Ana_CorAutoTTT1(:,1) ;
% VarAna1 = Ana_CorAutoTTT1_Fit(:,1) ;
% Var = vertcat(VarAna1) ;   Var = vertcat(Var1, VarAna1, Ana_CorAutoTTT2(:,1),Ana_CorAutoTTT2_Fit(:,1)) ;
% VarVar0 = max(max(Var)) - min(min(Var)) ;
% Var1 = (Var1 - min(min(Var))) ./ VarVar0 ;
% VarAna1 = (VarAna1 - min(min(Var))) ./ VarVar0 ;
% 
% Var2 = Ana_CorAutoTTT2(:,1) ;
% VarAna2 = Ana_CorAutoTTT2_Fit(:,1) ;
% Var = vertcat(VarAna2) ;   Var = vertcat(Ana_CorAutoTTT1(:,1),Ana_CorAutoTTT1_Fit(:,1), Var2, VarAna2) ;
% VarVar0 = max(max(Var)) - min(min(Var)) ;
% Var2 = (Var2 - min(min(Var))) ./ VarVar0 ;
% VarAna2 = (VarAna2 - min(min(Var))) ./ VarVar0 ;
% 
% Var = vertcat(Var1, VarAna1, Var2, VarAna2) ;
% semilogx(dT, Var1,'oc', dT, Var2,'om',...
%     Ana_FittedX, VarAna1, 'c', Ana_FittedX, VarAna2, 'm',...
%     'LineWidth',DataLineWidth, 'MarkerSize', 10) ;
% 
% % Var3 = Ana_CorAutoTTT3(:,1) ;
% % VarAna3 = Ana_CorAutoTTT3_Fit(:,1) ;
% % Var = vertcat(VarAna3) ; 
% % VarVar0 = max(max(Var)) - min(min(Var)) ;
% % Var3 = (Var3 - min(min(Var))) ./ VarVar0 ;
% % VarAna3 = (VarAna3 - min(min(Var))) ./ VarVar0 ;
% % 
% % Var = vertcat(Var1, VarAna1, Var2, VarAna2, Var3, VarAna3) ;
% % semilogx(dT, Var1,'oc', dT, Var2,'om', dT, Var3,'oy',...
% %     Ana_FittedX, VarAna1, 'c', Ana_FittedX, VarAna2, 'm', Ana_FittedX, VarAna3, 'y',...
% %     'LineWidth',DataLineWidth, 'MarkerSize', 10) ;

% Fig1 = figure;
% Fig1.OuterPosition = [500, 500, 800, 400] ;
% Var = vertcat(Ana_CorAutoTTT1(:,1), Ana_CorAutoTTT1_Fit(:,1), Ana_CorAutoTTT2_Fit(:,1)) ;
% semilogx(dT, Ana_CorAutoTTT1(:,1),'ok',...
%     Ana_FittedX,  Ana_CorAutoTTT1_Fit(:,1), 'r', Ana_FittedX,  Ana_CorAutoTTT2_Fit(:,1), 'b',...
%     'LineWidth',DataLineWidth, 'MarkerSize', 10) ;


title('Auto correation 1')
xlim([min(dT)*0.7 max(dT)*1.3]) ;
Var0 = max(max(Var)) - min(min(Var)) ;
Var1 = min(min(Var)) - Var0 *0.15 ;
Var2 = max(max(Var)) + Var0 *0.15 ;
ylim([Var1, Var2]) ;
xlabel('dT (s)','FontSize',FontSize,'FontWeight','bold','Color','k')
ylabel('Correlation','FontSize',FontSize,'FontWeight','bold','Color','k')
set(gca, 'tickdir', 'out') ; set(gca,'LineWidth',LineWidth); set(gca,'FontSize',FontSize); set(gca, 'ticklength', [0.01 0.01]) ;
legend(char(StrF_Legend_CorAuto),'Location','eastoutside') ;

pbaspect([25 10 1])
ax = gca; 
ax.TickLength = [0.02 0.02] ;
ax.XTick = [10^-4, 10^-3, 10^-2, 10^-1, 10^0, 10^1] ;
ax.XMinorTick = 'off' ;
ax.YMinorTick = 'off' ;
ax.XTickMode = 'manual'; 
ax.YTickMode = 'manual';
ax.ZTickMode = 'manual';
ax.XLimMode = 'manual';
ax.YLimMode = 'manual';
ax.ZLimMode = 'manual';
StrF = strcat('Cor_Auto1_WithFit','.pdf') ; print('-r1200','-dpdf',StrF);
StrF = strcat('Cor_Auto1_WithFit') ; saveas(gcf, StrF, 'png') ;
%close gcf ;

%%
Fig2 = figure;
Fig2.OuterPosition = [500, 500, 800, 400] ;
Var1 = Ana_CorAuto(:,2) ;
VarAna = Ana_CorAuto_Fit(:,2) ;
Var = vertcat(Var1, VarAna) ;
semilogx(dT, Var1,'ok',...
    Ana_FittedX, VarAna, 'r',...
    'LineWidth',DataLineWidth, 'MarkerSize', 10) ;

% VarAna = Ana_CorAuto_Fit(:,2) ;
% VarAna1 = VarAna;  VarAna1(:) = Ana_CorAuto_Fit(:,2) ;
% VarAna2 = VarAna;  VarAna2(:) = Ana_CorAuto_Fit_free(:,2) ;
% Var = vertcat(Var1, VarAna1, VarAna2) ;
% semilogx(dT, Var1,'ok',...  
%     Ana_FittedX, VarAna1, 'r',...
%     Ana_FittedX, VarAna2, 'b',...
%     'LineWidth',DataLineWidth, 'MarkerSize', 10) ;

% VarAna = Ana_CorAuto_Fit(:,2) ;
% VarAna1 = VarAna;  VarAna1(:) = Fit2stateFix_Ana_FittedY(2,2,:) ;
% VarAna2 = VarAna;  VarAna2(:) = Fit3stateFix_Ana_FittedY(2,2,:) ;
% VarAna3 = VarAna;  VarAna3(:) = Fit3state_Ana_FittedY(2,2,:) ;
% VarAna4 = VarAna;  VarAna4(:) = Fit3stateFixBG_Ana_FittedY(2,2,:) ;
% %VarAna5 = VarAna;  VarAna5(:) = Fit3stateFreeBG_Ana_FittedY(2,2,:) ;
% Var = vertcat(Var1, VarAna2, VarAna3, VarAna4) ;
% semilogx(dT, Var1,'ok',...
%     Ana_FittedX, VarAna1, 'b',...    
%     Ana_FittedX, VarAna2, 'r',...
%     Ana_FittedX, VarAna3, 'm',...
%     Ana_FittedX, VarAna4, 'g',...%Ana_FittedX, VarAna5, 'g',...
%     'LineWidth',DataLineWidth, 'MarkerSize', 10) ;

% Fig2 = figure;
% Fig2.OuterPosition = [500, 500, 800, 400] ;
% Var1 = Ana_CorAutoTTT1(:,2) ;
% VarAna1 = Ana_CorAutoTTT1_Fit(:,2) ;
% Var = vertcat(VarAna1) ;   Var = vertcat(Var1, VarAna1, Ana_CorAutoTTT2(:,2),Ana_CorAutoTTT2_Fit(:,2)) ;
% VarVar0 = max(max(Var)) - min(min(Var)) ;
% Var1 = (Var1 - min(min(Var))) ./ VarVar0 ;
% VarAna1 = (VarAna1 - min(min(Var))) ./ VarVar0 ;
% 
% Var2 = Ana_CorAutoTTT2(:,2) ;
% VarAna2 = Ana_CorAutoTTT2_Fit(:,2) ;
% Var = vertcat(VarAna2) ;   Var = vertcat(Ana_CorAutoTTT1(:,2),Ana_CorAutoTTT1_Fit(:,2), Var2, VarAna2) ;
% VarVar0 = max(max(Var)) - min(min(Var)) ;
% Var2 = (Var2 - min(min(Var))) ./ VarVar0 ;
% VarAna2 = (VarAna2 - min(min(Var))) ./ VarVar0 ;
% 
% Var = vertcat(Var1, VarAna1, Var2, VarAna2) ;
% semilogx(dT, Var1,'oc', dT, Var2,'om',...
%     Ana_FittedX, VarAna1, 'c', Ana_FittedX, VarAna2, 'm',...
%     'LineWidth',DataLineWidth, 'MarkerSize', 10) ;
% 
% % Var3 = Ana_CorAutoTTT3(:,2) ;
% % VarAna3 = Ana_CorAutoTTT3_Fit(:,2) ;
% % Var = vertcat(VarAna3) ;
% % VarVar0 = max(max(Var)) - min(min(Var)) ;
% % Var3 = (Var3 - min(min(Var))) ./ VarVar0 ;
% % VarAna3 = (VarAna3 - min(min(Var))) ./ VarVar0 ;
% % 
% % Var = vertcat(Var1, VarAna1, Var2, VarAna2, Var3, VarAna3) ;
% % semilogx(dT, Var1,'oc', dT, Var2,'om', dT, Var3,'oy',...
% %     Ana_FittedX, VarAna1, 'c', Ana_FittedX, VarAna2, 'm', Ana_FittedX, VarAna3, 'y',...
% %     'LineWidth',DataLineWidth, 'MarkerSize', 10) ;

% Fig2 = figure;
% Fig2.OuterPosition = [500, 500, 800, 400] ;
% Var = vertcat(Ana_CorAutoTTT1(:,2), Ana_CorAutoTTT1_Fit(:,2), Ana_CorAutoTTT2_Fit(:,2)) ;
% semilogx(dT, Ana_CorAutoTTT1(:,2),'ok',...
%     Ana_FittedX,  Ana_CorAutoTTT1_Fit(:,2), 'r', Ana_FittedX,  Ana_CorAutoTTT2_Fit(:,2), 'b',...
%     'LineWidth',DataLineWidth, 'MarkerSize', 10) ;


title('Auto correation 2')
xlim([min(dT)*0.7 max(dT)*1.3]) ;
Var0 = max(max(Var)) - min(min(Var)) ;
Var1 = min(min(Var)) - Var0 *0.15 ;
Var2 = max(max(Var)) + Var0 *0.15 ;
ylim([Var1, Var2]) ;
xlabel('dT (s)','FontSize',FontSize,'FontWeight','bold','Color','k')
ylabel('Correlation','FontSize',FontSize,'FontWeight','bold','Color','k')
set(gca, 'tickdir', 'out') ; set(gca,'LineWidth',LineWidth); set(gca,'FontSize',FontSize); set(gca, 'ticklength', [0.01 0.01]) ;
legend(char(StrF_Legend_CorAuto),'Location','eastoutside') ;

pbaspect([25 10 1])
ax = gca; 
ax.TickLength = [0.02 0.02] ;
ax.XTick = [10^-4, 10^-3, 10^-2, 10^-1, 10^0, 10^1] ;
ax.XMinorTick = 'off' ;
ax.YMinorTick = 'off' ;
ax.XTickMode = 'manual';
ax.YTickMode = 'manual';
ax.ZTickMode = 'manual';
ax.XLimMode = 'manual';
ax.YLimMode = 'manual';
ax.ZLimMode = 'manual';
StrF = strcat('Cor_Auto2_WithFit','.pdf') ; print('-r1200','-dpdf',StrF);
StrF = strcat('Cor_Auto2_WithFit') ; saveas(gcf, StrF, 'png') ;
%close gcf ;

%%
Fig3 = figure;
Fig3.OuterPosition = [500, 500, 800, 400] ;
Var1 = Ana_CorCross(:,1) ;
VarAna = Ana_CorCross_Fit(:,1) ;
Var = vertcat(Var1, VarAna) ;
semilogx(dT, Var1,'ok',...
    Ana_FittedX, VarAna, 'r',...
    'LineWidth',DataLineWidth, 'MarkerSize', 10) ;

% VarAna = Ana_CorCross_Fit(:,1) ;
% VarAna1 = VarAna;  VarAna1(:) = Ana_CorCross_Fit(:,1) ;
% VarAna2 = VarAna;  VarAna2(:) = Ana_CorCross_Fit_free(:,1) ;
% Var = vertcat(Var1, VarAna1, VarAna2) ;
% semilogx(dT, Var1,'ok',...  
%     Ana_FittedX, VarAna1, 'r',...
%     Ana_FittedX, VarAna2, 'b',...
%     'LineWidth',DataLineWidth, 'MarkerSize', 10) ;

% VarAna = Ana_CorCross_Fit(:,1) ;
% VarAna1 = VarAna;  VarAna1(:) = Fit2stateFix_Ana_FittedY(1,2,:) ;
% VarAna2 = VarAna;  VarAna2(:) = Fit3stateFix_Ana_FittedY(1,2,:) ;
% VarAna3 = VarAna;  VarAna3(:) = Fit3state_Ana_FittedY(1,2,:) ;
% VarAna4 = VarAna;  VarAna4(:) = Fit3stateFixBG_Ana_FittedY(1,2,:) ;
% %VarAna5 = VarAna;  VarAna5(:) = Fit3stateFreeBG_Ana_FittedY(1,2,:) ;
% Var = vertcat(Var1, VarAna2, VarAna3, VarAna4) ;
% semilogx(dT, Var1,'ok',...
%     Ana_FittedX, VarAna1, 'b',...    
%     Ana_FittedX, VarAna2, 'r',...
%     Ana_FittedX, VarAna3, 'm',...
%     Ana_FittedX, VarAna4, 'g',...%Ana_FittedX, VarAna5, 'g',...
%     'LineWidth',DataLineWidth, 'MarkerSize', 10) ;

% Fig3 = figure;
% Fig3.OuterPosition = [500, 500, 800, 400] ;
% Var1 = Ana_CorCrossTTT1(:,1) ;
% VarAna1 = Ana_CorCrossTTT1_Fit(:,1) ;
% Var = vertcat(VarAna1) ;   %Var = vertcat(Var1, VarAna1, Ana_CorCrossTTT2(:,1),Ana_CorCrossTTT2_Fit(:,1)) ;
% VarVar0 = max(max(Var)) - min(min(Var)) ;
% Var1 = (Var1 - min(min(Var))) ./ VarVar0 ;
% VarAna1 = (VarAna1 - min(min(Var))) ./ VarVar0 ;
% 
% Var2 = Ana_CorCrossTTT2(:,1) ;
% VarAna2 = Ana_CorCrossTTT2_Fit(:,1) ;
% Var = vertcat(VarAna2) ;   %Var = vertcat(Ana_CorCrossTTT1(:,1),Ana_CorCrossTTT1_Fit(:,1), Var2, VarAna2) ;
% VarVar0 = max(max(Var)) - min(min(Var)) ;
% Var2 = (Var2 - min(min(Var))) ./ VarVar0 ;
% VarAna2 = (VarAna2 - min(min(Var))) ./ VarVar0 ;
% 
% Var = vertcat(Var1, VarAna1, Var2, VarAna2) ;
% semilogx(dT, Var1,'oc', dT, Var2,'om',...
%     Ana_FittedX, VarAna1, 'c', Ana_FittedX, VarAna2, 'm',...
%     'LineWidth',DataLineWidth, 'MarkerSize', 10) ;
% 
% % Var3 = Ana_CorCrossTTT3(:,1) ;
% % VarAna3 = Ana_CorCrossTTT3_Fit(:,1) ;
% % Var = vertcat(VarAna3) ;
% % VarVar0 = max(max(Var)) - min(min(Var)) ;
% % Var3 = (Var3 - min(min(Var))) ./ VarVar0 ;
% % VarAna3 = (VarAna3 - min(min(Var))) ./ VarVar0 ;
% % 
% % Var = vertcat(Var1, VarAna1, Var2, VarAna2, Var3, VarAna3) ;
% % semilogx(dT, Var1,'oc', dT, Var2,'om', dT, Var3,'oy',...
% %     Ana_FittedX, VarAna1, 'c', Ana_FittedX, VarAna2, 'm', Ana_FittedX, VarAna3, 'y',...
% %     'LineWidth',DataLineWidth, 'MarkerSize', 10) ;

% Fig3 = figure;
% Fig3.OuterPosition = [500, 500, 800, 400] ;
% Var = vertcat(Ana_CorCrossTTT1(:,1), Ana_CorCrossTTT1_Fit(:,1), Ana_CorCrossTTT2_Fit(:,1)) ;
% semilogx(dT, Ana_CorCrossTTT1(:,1),'ok',...
%     Ana_FittedX,  Ana_CorCrossTTT1_Fit(:,1), 'r', Ana_FittedX,  Ana_CorCrossTTT2_Fit(:,1), 'b',...
%     'LineWidth',DataLineWidth, 'MarkerSize', 10) ;


title('Cross correation 1-2')
xlim([min(dT)*0.7 max(dT)*1.3]) ;
Var0 = max(max(Var)) - min(min(Var)) ;
Var1 = min(min(Var)) - Var0 *0.15 ;
Var2 = max(max(Var)) + Var0 *0.15 ;
ylim([Var1, Var2]) ;
xlabel('dT (s)','FontSize',FontSize,'FontWeight','bold','Color','k')
ylabel('Correlation','FontSize',FontSize,'FontWeight','bold','Color','k')
set(gca, 'tickdir', 'out') ; set(gca,'LineWidth',LineWidth); set(gca,'FontSize',FontSize); set(gca, 'ticklength', [0.01 0.01]) ;
legend(char(StrF_Legend_CorAuto),'Location','eastoutside') ;

pbaspect([25 10 1])
ax = gca; 
ax.TickLength = [0.02 0.02] ;
ax.XTick = [10^-4, 10^-3, 10^-2, 10^-1, 10^0, 10^1] ;
ax.XMinorTick = 'off' ;
ax.YMinorTick = 'off' ;
ax.XTickMode = 'manual';
ax.YTickMode = 'manual';
ax.ZTickMode = 'manual';
ax.XLimMode = 'manual';
ax.YLimMode = 'manual';
ax.ZLimMode = 'manual';
StrF = strcat('Cor_Cross1_WithFit','.pdf') ; print('-r1200','-dpdf',StrF);
StrF = strcat('Cor_Cross1_WithFit') ; saveas(gcf, StrF, 'png') ;
%close gcf ;


%%

% Var = (NumOfState+1)*NumOfState/2 ;
% Ana_Correlation_each = zeros(length(dT), Var) ;
% Var = size(Ana_Correlation) ;
% Count = Count + 1 ;
% I = 0 ;
% while I < Var(1)
%     I = I + 1 ;
%     K = I-1 ;
%     while K < Var(2)
%         K = K + 1 ;
%         Count = Count + 1 ;
%         Ana_Correlation_each(:,Count) = Ana_Correlation ;
%     end
% end



%%
% TargetM = RisePoint_Result_1DFDC_Mat_M_Model_lin ;
% 
% RisePointIRF = 307 ;
% RisePointAveRange = 13 ;
% 
% 
% Var = RisePoint_Result_1DFDC_EstQ(:,2) ;
% [Var, Icenter] = min(abs(Var-RisePointIRF)) ;
% 
% Var = size(RisePoint_Result_1DFDC_Mat_M_Model_lin) ;
% Ave_1DFDC_Mat_M_Model_lin = zeros(Var(1),1) ;
% A = zeros(Var(1),2) ;
% 
% Count = 0 ;
% I = Icenter - round(RisePointAveRange/2) ;
% while Count < RisePointAveRange
%     Count = Count + 1 ;
%     I = I + 1 ;
%     Ave_1DFDC_Mat_M_Model_lin = Ave_1DFDC_Mat_M_Model_lin + RisePoint_Result_1DFDC_Mat_M_Model_lin(:,I) ;
% end
% Ave_1DFDC_Mat_M_Model_lin = Ave_1DFDC_Mat_M_Model_lin / Count ;
% 
% A(:,1) = Mat_1DFDC_lin ;
% A(:,2) = Ave_1DFDC_Mat_M_Model_lin ;







