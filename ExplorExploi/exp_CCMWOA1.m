% for i = 1:10
clear all
close all
clc
addpath(genpath(pwd));

%%parm setting
% option setting
SearchAgents_no=30; % Number of search agents
Div_algorithm_num=2;
NumofRecord=40;
dim=30;
MaxFEs=dim*5000*2;
algorithmName={'LRRIME', 'RIME'};
% Flod=30;
Flod=4;
% %% prepare xlsfile
% timestr=datestr(now,'mm-dd_HH_MM');
% dirname=['exp_result/CEC17-',algorithmName{1},'-',algorithmName{2},'-',timestr,'-Dim',num2str(dim)];
% % dirname=['exp_result/Functiontest-',timestr,'-Dim',num2str(dim)];
% mkdir(dirname);
% filename=[dirname,'/FES-',timestr];

%% perpare xlsfile
timestr=datestr(now,'mm-dd_HH_MM');
% dirname=['exp_result/11.18-',algorithmName{1},'-',algorithmName{2},'-',timestr,'-Dim',num2str(dim)];
dirname=['exp_result/', algorithmName{1},'-',timestr];
mkdir(dirname);
% Excel filename.
% filename=[dirname,'/jona-',timestr];
filename=[dirname,'/', algorithmName{1}];

te= {'F','Algrithm','max','min','mean','std'};
xlsfilename=[filename,'.xlsx'];
xlswrite(xlsfilename, te, 'overall')
startLineNum=2;% the startLineNum of overall sheet to write data


algrithmNum=size(algorithmName,2);
% generate the space of linestyles, MarkerEdgeColors,Markers
nLines = algrithmNum;
basic_linestyles = cellstr(char('-',':','-.','--'));
basic_Markers    = cellstr(char('o','x','+','*','s','d','v','^','<','>','p','h','.'));
MarkerEdgeColors = hsv(nLines);
linestyles       = repmat(basic_linestyles,ceil(nLines/numel(basic_linestyles)),1);
Markers          = repmat(basic_Markers,ceil(nLines/numel(basic_Markers)),1);
%% =======uncomment for specify function===============
% 'CEC14'
% Function_name_all={'F24','F25','F26','F27','F28','F29','F30','F31','F32','F33','F34','F35','F36','F37','F38','F39','F40','F41','F42','F43','F44','F45','F46','F47','F48','F49','F50','F51','F52','F53'};
% Function_name_all={'F1','F2','F3','F4','F5','F6','F7','F8','F9','F10','F11','F12','F13','F14','F15','F16','F17','F18','F19','F20','F21','F22','F23'};
% Function_name_all={'F107','F108','F109','F110','F111','F112','F113','F114','F115','F116','F117','F118','F119','F120','F121','F122','F123','F124','F125','F126','F127','F128','F129','F130','F131','F132','F133','F134','F135','F136'};
% Function_name_all={ 'F79','F80','F81','F82','F83','F84','F85','F86','F87','F88','F89','F90','F91','F92','F93','F94','F95','F96','F97','F98','F99','F100','F101','F102','F103','F104','F105','F106'};
% Function_name_all={ 'F54','F55','F56','F57','F58','F59','F60','F61','F62','F63','F64','F65','F66','F67','F68','F69','F70','F71','F72','F73','F74','F75','F76','F77','F78'};
% Function_name_all={'F1','F2','F3','F4','F5','F6','F7','F8','F9','F10','F11','F12','F13','F14','F15','F16','F17','F18','F19','F20','F21','F22','F23',...
%     'F24','F25','F26','F27','F28','F29','F30','F31','F32','F33','F34','F35','F36','F37','F38','F39','F40','F41','F42','F43','F44','F45','F46','F47','F48','F49','F50','F51','F52','F53',...
%     'F107','F108','F109','F110','F111','F112','F113','F114','F115','F116','F117','F118','F119','F120','F121','F122','F123','F124','F125','F126','F127','F128','F129','F130','F131','F132','F133','F134','F135','F136'};
% Function_name_all={'F1','F2','F3','F4','F5','F6','F7','F8','F9','F10','F11','F12','F13','F14','F15','F16','F17','F18','F19','F20','F21','F22','F23',...
%     };
% Function_name_all={'F24','F25','F26','F27','F28','F29','F30','F31','F32','F33','F34','F35','F36','F37','F38','F39','F40','F41','F42','F43','F44','F45','F46','F47','F48','F49','F50','F51','F52','F53',...
%     };

%% CEC2017
Function_name_all={'F107','F109','F110','F111','F112','F113','F114','F115','F116','F117','F118','F119','F120','F121','F122','F123','F124','F125','F126','F127','F128','F129','F130','F131','F132','F133','F134','F135','F136'};
% %%
% Function_name_all = {'F107'};
% Function_name_all={'F137','F138','F139','F140','F141','F142','F143','F144','F145','F146'};
% Function_name_all = {'F107', 'F111', 'F114','F127', 'F130', 'F132'};

for funcNum=1:size(Function_name_all,2)
    dim=30;
    Function_name=Function_name_all{funcNum};
    [lb,ub,dim,fhd]=Get_Functions(Function_name,dim);
    display(['----------------',Function_name,'--------------------']);
    Function_name=['F',num2str(funcNum)];
    %% =================uncomment for 1-23 function==========
    % for funcNum=109:146
    %     dim=30;
    %     Function_name=['F',num2str(funcNum)]; % Name of the test function that can be from F1 to F23
    %     [lb,ub,dim,fhd]=Get_Functions(Function_name,dim);
    %     display(['----------------',Function_name,'--------------------']);
    %%======================================================
    % benchmark function
    cg_curves=zeros(algrithmNum,Flod,NumofRecord);
%     %parfor cflod=1:Flod   [≥ı º∞Ê±æ]
    parfor cflod=1:Flod
%     for cflod = 1:Flod
        display(['flod',num2str(cflod)]);
        for cnum=1:algrithmNum
            alg_fhd=str2func(algorithmName{cnum});
            [average_distance,div,~,cg_curve]=alg_fhd(SearchAgents_no,MaxFEs,lb,ub,dim,fhd);
            cg_curves(cnum,cflod,:)=MySampling(cg_curve,NumofRecord);
            if cnum<=Div_algorithm_num
                Div(cnum,cflod,:)=MySampling(div,10000);
                Average_distance(cnum,cflod,:)=MySampling(average_distance,10000);
            end
        end
    end
    temp=zeros(Flod,10000);
    temp1=zeros(Flod,10000);
    Div_mean=zeros(Div_algorithm_num,10000);
    Average_distance_mean=zeros(Div_algorithm_num,10000);
    for i=1:Div_algorithm_num
        for j=1:Flod
            temp(j,:)=Div(i,j,:);
        end
        Div_mean(i,:)=mean(temp);
        for j=1:Flod
            temp1(j,:)=Average_distance(i,j,:);
        end
        Average_distance_mean(i,:)=mean(temp1);
    end
    %% write data to excel sheet
    algName_labels=cell(algrithmNum*Flod,1);
    all_curves=zeros(algrithmNum*Flod,NumofRecord);
    folds_labels=repmat(int32(1:Flod)',[algrithmNum,1]);
    for it=1:algrithmNum
        algName_labels((it-1)*Flod+1:(it-1)*Flod+Flod)=algorithmName(it);
        all_curves((it-1)*Flod+1:(it-1)*Flod+Flod,:)=cg_curves(it,:,:);
    end
    xlswrite(xlsfilename, algName_labels, Function_name, 'A1')
    xlswrite(xlsfilename, folds_labels,Function_name, 'B1')
    xlswrite(xlsfilename, all_curves, Function_name, 'C1')
    %write overall
    algName_labels2=algorithmName';
    statistic_values=zeros(algrithmNum,4);
    for it=1:algrithmNum
        statistic_values(it,:)=[max(cg_curves(it,:,end)),min(cg_curves(it,:,end)),mean(cg_curves(it,:,end)),std(cg_curves(it,:,end))];
    end
    funcNum_label=repmat({Function_name},algrithmNum,1);
    xlswrite(xlsfilename, funcNum_label, 'overall', ['A',num2str(startLineNum)])
    xlswrite(xlsfilename, algName_labels2, 'overall', ['B',num2str(startLineNum)])
    xlswrite(xlsfilename, statistic_values, 'overall', ['C',num2str(startLineNum)])
    startLineNum=startLineNum+algrithmNum;
    %% plot cg_curveline
    clf
    set(gcf,'Position',[0,0,1000,600])
    
    for it=1:algrithmNum
        yy(it,:)=mean(all_curves((it-1)*Flod+1:(it-1)*Flod+Flod,:));
    end
    xx=[1:NumofRecord]*(MaxFEs/NumofRecord);
    for it=1:algrithmNum
        % semilogy(xx,yy(it,:),[linestyles{it} Markers{it}],'LineWidth', 1.5,'Color',MarkerEdgeColors(it,:));
        semilogy(xx,yy(it,:),[linestyles{it}],'LineWidth', 1.5,'Color',MarkerEdgeColors(it,:));
        hold on;
    end
    hold off;
    title(Function_name);
    set(gcf,'color','white')

    xlabel('FEs');
    ylabel('Best score obtained so far');
    algrithmName1=strrep(algorithmName,'_','\_');
    
    legend(algrithmName1,'Location', 'NorthEast');
    legend('boxoff')
    
    a=findobj(gcf); % get the handles associated with the current figure
    allaxes=findall(a,'Type','axes');
    % alllines=findall(a,'Type','line');
    alltext=findall(a,'Type','text');
    set(allaxes,'FontName','Times','LineWidth',1,'FontSize',14,'FontWeight','bold');
    % set(alllines,'Linewidth',1);
    set(alltext,'FontName','Times','FontSize',14,'FontWeight','bold')
    %
    set(gcf, 'PaperUnits', 'inches');
    % set(gcf, 'PaperUnits', 'centimeters');
    krare=3.5;
    % x_width=krare*1.618 ;
    x_width=krare*5/3 ;
    %  x_width=3*1;
    y_width=krare*4/3;
    set(gcf, 'PaperPosition', [0 0 x_width y_width]); %
    % set(gcf,'PaperUnits','inches','PaperPosition',[0 0 5 4])
    set(gca, ...
        'Box'         , 'on'     , ...
        'TickDir'     , 'in'     , ...
        'TickLength'  , [.02 .02] , ...
        'XMinorTick'  , 'on'      , ...
        'YMinorTick'  , 'on'      , ...
        'YGrid'       , 'off'      , ...
        'XGrid'       , 'off'      , ...
        'XColor'      , [.3 .3 .3], ...
        'YColor'      , [.3 .3 .3], ...
        'LineWidth'   , 1         );
    axis tight
%     grid on
%     box on
    saveas(gcf,[filename,'-',Function_name,'-carve'],'fig')
    filename1= [filename,'-',Function_name,'-carve'];
    print(filename1,'-dtiff', '-r300'); %<-Save as PNG with 300 DPI
        %% plot balance
    
    for i=1:Div_algorithm_num
        clf
        yy1=smooth(Div_mean(i,:),'loess');
        Div_max=max(yy1);
        XPL=zeros(1,size(yy1,1));
        XPT=zeros(1,size(yy1,1));
        save=zeros(1,size(yy1,1)+1);
        Increment=zeros(1,size(yy1,1));
        D_value=0;
        for j =1:size(yy1,1)
            XPL(j)=(yy1(j)/Div_max)*100;
            XPT(j)=(abs(yy1(j)-Div_max)/Div_max)*100;
            D_value=D_value+XPL(j)-XPT(j);
            if D_value<0
                D_value=0;
            end
            save(j)=D_value;
        end
        save=[0 save];
        D_value_max=max(save);
        for j=1:size(save,2)
        Increment(j)=(save(j)/D_value_max)*100;
        end
         %  Data=[XPT;XPL;Increment];
        Color=hsv(3);
        XPL_mean=mean(XPL);
        XPT_mean=mean(XPT);
    
        hold on
        % a1=[rand rand rand]
        %'LineStyle','--',
        XPL(10000+1)=nan;
        patch(1:size(XPL,2),XPL,Color(1,:),'edgealpha',1,'facealpha',0,'LineWidth',2,'edgecolor',Color(1,:));
        XPT(10000+1)=nan;
        patch(1:size(XPT,2),XPT,Color(3,:),'edgealpha',1,'facealpha',0,'LineWidth',2,'edgecolor',Color(3,:));
        Increment(1001+1)=nan;
        patch(1:size(Increment,2),Increment,Color(2,:),'edgealpha',1,'facealpha',0,'LineWidth',2,'edgecolor',Color(2,:));
        
        Title=['Balance analysis of ',algorithmName{i}];
        title(Title,'interpreter','Tex',...
            'FontSize',15,...
            'FontWeight','bold',...
            'FontName','Times')

        xlabel({'Iteration'},...
            'FontUnits','points',...
            'interpreter','Tex',...
            'FontSize',15,...
            'FontWeight','bold',...
            'FontName','Times')
        ylabel({'Percentage'},...
            'FontUnits','points',...
            'interpreter','Tex',...
            'FontSize',15,...
            'FontWeight','bold',...
            'FontName','Times')
        
        Data_name={['Exploration(Avg, ',num2str(XPL_mean),'%)'],['Exploitation(Avg, ',num2str(XPT_mean),'%)'],'Incremental-Decremental'};
        legend(Data_name,'Location', 'NorthEast');
        legend('boxoff') 
        % xlabel('Iteration#')
        box on
        axis tight

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        a=findobj(gcf); % get the handles associated with the current figure

        allaxes=findall(a,'Type','axes');
        % alllines=findall(a,'Type','line');
        alltext=findall(a,'Type','text');

        set(allaxes,'FontName','Times','LineWidth',1,'FontSize',15,'FontWeight','bold');
        %  set(alllines,'Linewidth',1);
        set(alltext,'FontName','Times','FontSize',15,'FontWeight','bold')



        set(gcf, 'PaperUnits', 'inches');

        krare=3;
        % x_width=krare*1.618 ;
        x_width=krare*5/3 ;

        %  x_width=3*1;

        y_width=krare*4/3;
        set(gcf, 'PaperPosition', [0 0 x_width y_width]); %
        % set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3])

        set(gca, ...
            'Box'         , 'on'     , ...
            'TickDir'     , 'in'     , ...
            'TickLength'  , [.02 .02] , ...
            'XMinorTick'  , 'on'      , ...
            'YMinorTick'  , 'on'      , ...
            'YGrid'       , 'off'      , ...
            'XGrid'       , 'off'      , ...
            'XColor'      , [.3 .3 .3], ...
            'YColor'      , [.3 .3 .3], ...
            'LineWidth'   , 1         );
        axis tight
        grid on
        box on
        saveas(gcf,[filename,'-',Function_name,'-balance','-',num2str(i)],'fig')
        filename1= [filename,'-',Function_name,'-balance','-',num2str(i)];
        print(filename1,'-dtiff', '-r300'); %<-Save as PNG with 300 DPI
    end
    %% plot diversity
    clf
    Color=hsv(Div_algorithm_num);
    for i =1:Div_algorithm_num
        Average_distance_mean(i,10000+1)=nan;
        patch(1:size(Average_distance_mean(i,:),2),Average_distance_mean(i,:),Color(i,:),'edgealpha',1,'facealpha',0,'LineWidth',2,'edgecolor',Color(i,:));
    end
    title('Algorithm diversity analysis','interpreter','Tex',...
            'FontSize',15,...
            'FontWeight','bold',...
            'FontName','Times')

        xlabel({'Iteration'},...
            'FontUnits','points',...
            'interpreter','Tex',...
            'FontSize',15,...
            'FontWeight','bold',...
            'FontName','Times')
        ylabel({'Average distance between search agent'},...
            'FontUnits','points',...
            'interpreter','Tex',...
            'FontSize',15,...
            'FontWeight','bold',...
            'FontName','Times')
        legend(algorithmName{1:Div_algorithm_num},'Location', 'NorthWest');
        legend('boxoff') 
        % xlabel('Iteration#')
        box on
        axis tight

%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        a=findobj(gcf); % get the handles associated with the current figure

        allaxes=findall(a,'Type','axes');
        % alllines=findall(a,'Type','line');
        alltext=findall(a,'Type','text');

        set(allaxes,'FontName','Times','LineWidth',1,'FontSize',15,'FontWeight','bold');
        %  set(alllines,'Linewidth',1);
        set(alltext,'FontName','Times','FontSize',15,'FontWeight','bold')



        set(gcf, 'PaperUnits', 'inches');

        krare=3;
        % x_width=krare*1.618 ;
        x_width=krare*5/3 ;

        %  x_width=3*1;

        y_width=krare*4/3;
        set(gcf, 'PaperPosition', [0 0 x_width y_width]); %
        % set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3])

        set(gca, ...
            'Box'         , 'on'     , ...
            'TickDir'     , 'in'     , ...
            'TickLength'  , [.02 .02] , ...
            'XMinorTick'  , 'on'      , ...
            'YMinorTick'  , 'on'      , ...
            'YGrid'       , 'off'      , ...
            'XGrid'       , 'off'      , ...
            'XColor'      , [.3 .3 .3], ...
            'YColor'      , [.3 .3 .3], ...
            'LineWidth'   , 1         );
        axis tight
        grid on
        box on
        saveas(gcf,[filename,'-',Function_name,'-diversity'],'fig')
        filename1= [filename,'-',Function_name,'-diversity'];
        print(filename1,'-dtiff', '-r300'); %<-Save as PNG with 300 DPI
end

Orderhao(xlsfilename);
pValueToExcelhao(xlsfilename,Flod);
FridTest3( xlsfilename ,Flod)
FridTest4( xlsfilename ,Flod)

close all
% end