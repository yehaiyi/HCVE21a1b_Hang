%% compare escape time
%load('mean_escape_time_1b.mat')
%load('escape_time_1a_500.mat')

load('escape_time_1a_500.mat')
load("escape_residues_time_1b.mat")
load("escape_residues_time_3a.mat")
run startup.m

rng default

escape_residues_1a=[384   386   388   390   391   393   394   395   396   397   398   399   400 ...
                    401   402   403   404   405   407   408   410   415   416   417   422   424 ...
                    431   433   434   435   438   442   444   446   453   456   461   466   475 ...
                    482   501   524   528   531   533   538   557   558   560   580   608 ...
                    610   636   713]; %除去520，一共54个escape residues
escape_residues_time_1a=mean_escape_time(escape_residues_1a-383);

data = [escape_residues_time_1a escape_residues_time_1b  escape_residues_time_3a ];
G = [zeros(size(escape_residues_time_1a))-1 zeros(size(escape_residues_time_1b)) ones(size(escape_residues_time_3a))];


set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultTextFontName','Arial')
set(0,'DefaultAxesFontSize',8)
set(0,'DefaultTextFontSize',8)
box_lineWidth = 0.75;
box_widths_value = 0.5;
black = [0 0 0];
box_color = [black;black;black];
box_color_transparency = 0; %faceAlpha
median_lineWidth = 0.75;
median_color = 'k';
whisker_value = 1.5;
outlier_marker = '';
outlier_markerSize = 3.5;
outlier_marker_edgeWidth = 0.001;
outlier_marker_edgeColor = 'w';
label_xaxis_data = {'1a','1b',sprintf('3a')};
text_ylabel = {'Escape time of',' escape residues'};
text_xlabel = '';
text_title = '';%'E2-escape mutations [Keck2009],[Morin2012],[Bailey2015]';
label_orientation_choice = 'horizontal'; %'horizontal'
ylim_min = 0;
ylim_max = 0.2;
outlier_jitter_value = 0;
savefig = 0;
savefig_name = 'escape_time_compare';
fig_width_cm = 4;
fig_height_cm =5;
FIG=figure;
set(gcf,'renderer','Painters')
x0=0.8+0.4*(rand(length(escape_residues_time_1a),1));
x1=1.8+0.4*(rand(length(escape_residues_time_1b),1));
x2=2.8+0.4*(rand(length(escape_residues_time_3a),1));
size_marker = 20;
% green = [27 129 62]/255;
% lightgreen = (1-green)*0.6+green;
f0=scatter(x0,escape_residues_time_1a ,'o','MarkerEdgeColor','k','MarkerFaceColor',red,'SizeData',size_marker,'LineWidth',0.01);f0.MarkerFaceAlpha = 0.6;hold on 
f1=scatter(x1,escape_residues_time_1b ,'o','MarkerEdgeColor','k','MarkerFaceColor',purple,'SizeData',size_marker,'LineWidth',0.01);f1.MarkerFaceAlpha = 0.6;hold on 
f2=scatter(x2,escape_residues_time_3a,'o','MarkerEdgeColor','k','MarkerFaceColor',orange,'SizeData',size_marker,'LineWidth',0.01);f2.MarkerFaceAlpha = f1.MarkerFaceAlpha;hold on



figure_boxplot(data,G,...
    box_lineWidth,box_widths_value,box_color,box_color_transparency,...
    median_lineWidth,median_color,...
    whisker_value,...
    outlier_marker,outlier_markerSize,outlier_marker_edgeWidth,outlier_marker_edgeColor,outlier_jitter_value,...
    label_xaxis_data,text_ylabel,text_xlabel,text_title,label_orientation_choice,...
    ylim_min,ylim_max,...
    savefig,savefig_name,fig_width_cm,fig_height_cm);

% V=violinplot(data, [repmat("1a",1,length(mean_escape_time)) repmat("1b",1,length(all_mean_escape_time_1b))],'EdgeColor' ,[0 0 0],'BoxColor' ,[0 0 0]);
% SizeData=10;
% ylabel({'Escape time of',' exposed residues'});
% V(1, 1).ViolinColor = purple;
% V(1, 2).ViolinColor = orange;
% V(1, 1).EdgeColor = 'None';
% 
% V(1, 1).MedianPlot.SizeData  =SizeData+30;
% V(1, 2).MedianPlot.SizeData  =SizeData+30;
% V(1, 2).EdgeColor = 'None';
% V(1, 1).ScatterPlot.SizeData  =SizeData;
% V(1, 2).ScatterPlot.SizeData  =SizeData;


set(gca,'YTick',0:50:250)
yt = get(gca, 'YTick');
axis([xlim    0  600])
xt = get(gca, 'XTick');
hold on
%plot(xt([1 2]), [1 1]*550, '-k','LineWidth',0.3)
%plot(xt([2 3]), [1 1]*550, '-k','LineWidth',0.3)
% plot(xt([1 1]), [0.95 1]*550, '-k','LineWidth',0.3)
% plot(xt([2 2]), [0.95 1]*550, '-k','LineWidth',0.3)
P1 = ranksum(escape_residues_time_3a,escape_residues_time_1a,'tail','left')
P2 = ranksum(escape_residues_time_3a,escape_residues_time_1b,'tail','left')

ind1 = floor(log10(P1));
P1 = roundn(P1/10^ind1,-1);
t1 = ['$$ P1 = ',num2str(P1),' \times 10^{',num2str(ind1),'} $$'];

ind2 = floor(log10(P2));
P2 = roundn(P2/10^ind2,-1);
t2 = ['$$ P2 = ',num2str(P2),' \times 10^{',num2str(ind2),'} $$'];

%text(1.3,560,t1,'interpreter','latex','FontSize',15,'FontName','Arial')
%text(2.3,560,t2,'interpreter','latex','FontSize',15,'FontName','Arial')

FIG.Name = 'escape_time_compare_1b3a';
set(gca,'fontname','Arial')  % Set it to times

FIG.Units = 'centimeters';
% set(gcf,'Position',[10 10 8 8.58]);
set(gcf,'Position',[6.53 6.53 4 4.16]);
% set(gcf,'Position',[6.53 6.53 5.22 3.92]);
set(gca,'TickLength',[0.035, 0.01])
set(gca,'Position',[.37 .12 .6 .86]);  %调整 XLABLE和YLABLE不会被切掉
figure_FontSize=8;
set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
set(gca,'TickDir','out')
set(gca,'TickLength',[0.02, 0.01])
set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
set(get(gca,'YLabel'), 'Units', 'Normalized', 'Position', [-0.5, 0.5, 0]);
set(findobj('FontSize',10),'FontSize',figure_FontSize);
title("escape time for escape residues")
xlabel("subtype")
% set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);
% try print(['C:\Users\27909\Desktop\' FIG.Name],'-dpng','-r600'); catch      print(['C:\Users\hzhangbr\Desktop\' FIG.Name],'-dpng','-r600');     end
% 
% try print(['C:\Users\27909\Desktop\' FIG.Name],'-dpng','-r600'); catch      print(['C:\Users\hzhangbr\Desktop\' FIG.Name],'-dpdf');     end
%% Known Escape residues
load("known_escape.mat")

run startup.m
label_xaxis_data = {'1a','1b','3a'};

%懒得改变量名了
HVR12_1a=known_escape(1,:);
HVR12_1b=known_escape(2,:);HVR12_1b=HVR12_1b(HVR12_1b<400);
HVR12_3a=known_escape(3,:);

NonVRs_1a= remain_1a1b(1,:);
NonVRs_1b= remain_1a1b(2,:);
NonVRs_3a= remain_3a;
% escape 
subplot(1,2,1)
data = [HVR12_1a HVR12_1b  HVR12_3a];
G = [zeros(size(HVR12_1a))-1 zeros(size(HVR12_1b)) ones(size(HVR12_3a))];
%figure;
max_1a = prctile(HVR12_1a,75)+1.5*(prctile(HVR12_1a,75)-prctile(HVR12_1a,25));
max_1b = prctile(HVR12_1b,75)+1.5*(prctile(HVR12_1b,75)-prctile(HVR12_1b,25));
max_3a = prctile(HVR12_3a,75)+1.5*(prctile(HVR12_3a,75)-prctile(HVR12_3a,25));

boxplot(data,G,'whisker',1.5,'labels',label_xaxis_data);hold on
x0=0.8+0.4*(rand(length(HVR12_1a(HVR12_1a<=max_1a)),1));
x1=1.8+0.4*(rand(length(HVR12_1b(HVR12_1b<=max_1b)),1));
x2=2.8+0.4*(rand(length(HVR12_3a(HVR12_3a<=max_3a)),1));
size_marker = 10;
f0=scatter(x0,HVR12_1a(HVR12_1a<=max_1a) ,'o','MarkerEdgeColor','k','MarkerFaceColor',red,'SizeData',size_marker,'LineWidth',0.01);f0.MarkerFaceAlpha = 0.6;hold on
f1=scatter(x1,HVR12_1b(HVR12_1b<=max_1b) ,'o','MarkerEdgeColor','k','MarkerFaceColor',purple,'SizeData',size_marker,'LineWidth',0.01);f1.MarkerFaceAlpha = 0.6;hold on
f2=scatter(x2,HVR12_3a(HVR12_3a<=max_3a),'o','MarkerEdgeColor','k','MarkerFaceColor',orange,'SizeData',size_marker,'LineWidth',0.01);f2.MarkerFaceAlpha = f1.MarkerFaceAlpha;hold on
ylabel("Escape time")
xlabel("Subtype")
ylim([0,500])
title("Known escape residues")

% remaining (没完全算完)
subplot(1,2,2)
data = [NonVRs_1a NonVRs_1b  NonVRs_3a ];
G = [zeros(size(NonVRs_1a))-1 zeros(size(NonVRs_1b)) ones(size(NonVRs_3a))];
%figure;
max_1a = prctile(NonVRs_1a,75)+1.5*(prctile(NonVRs_1a,75)-prctile(NonVRs_1a,25));
max_1b = prctile(NonVRs_1b,75)+1.5*(prctile(NonVRs_1b,75)-prctile(NonVRs_1b,25));
max_3a = prctile(NonVRs_3a,75)+1.5*(prctile(NonVRs_3a,75)-prctile(NonVRs_3a,25));

boxplot(data,G,'whisker',1.5,'labels',label_xaxis_data);hold on
x0=0.8+0.4*(rand(length(NonVRs_1a(NonVRs_1a<=max_1a)),1));
x1=1.8+0.4*(rand(length(NonVRs_1b(NonVRs_1b<=max_1b)),1));
x2=2.8+0.4*(rand(length(NonVRs_3a(NonVRs_3a<=max_3a)),1));
size_marker = 10;
f0=scatter(x0,NonVRs_1a(NonVRs_1a<=max_1a) ,'o','MarkerEdgeColor','k','MarkerFaceColor',red,'SizeData',size_marker,'LineWidth',0.01);f0.MarkerFaceAlpha = 0.6;hold on
f1=scatter(x1,NonVRs_1b(NonVRs_1b<=max_1b) ,'o','MarkerEdgeColor','k','MarkerFaceColor',purple,'SizeData',size_marker,'LineWidth',0.01);f1.MarkerFaceAlpha = 0.6;hold on
f2=scatter(x2,NonVRs_3a(NonVRs_3a<=max_3a),'o','MarkerEdgeColor','k','MarkerFaceColor',orange,'SizeData',size_marker,'LineWidth',0.01);f2.MarkerFaceAlpha = f1.MarkerFaceAlpha;hold on
ylabel("Escape time")
xlabel("Subtype")
ylim([0,500])
title("Remaining")

%[p,h] = ranksum(escape_residues_time_3a,escape_residues_time_1a,"tail","left")
%[p,h] = ranksum(escape_residues_time_1b,escape_residues_time_3a,"tail","left")
%% 部分remaining 和 全部remaining的对比

load("known_escape.mat")

run startup.m
label_xaxis_data = {'1a','1b'};

%懒得改变量名了
HVR12_1a=remain_1a1b(1,:);
HVR12_1b=remain_1a1b(2,:);

NonVRs_1a= remain_all(1,:);
NonVRs_1b= remain_all(2,:);

% escape 
subplot(1,2,1)
data = [HVR12_1a HVR12_1b];
G = [zeros(size(HVR12_1a))-1 zeros(size(HVR12_1b))];
%figure;
max_1a = prctile(HVR12_1a,75)+1.5*(prctile(HVR12_1a,75)-prctile(HVR12_1a,25));
max_1b = prctile(HVR12_1b,75)+1.5*(prctile(HVR12_1b,75)-prctile(HVR12_1b,25));

boxplot(data,G,'whisker',1.5,'labels',label_xaxis_data);hold on
x0=0.8+0.4*(rand(length(HVR12_1a(HVR12_1a<=max_1a)),1));
x1=1.8+0.4*(rand(length(HVR12_1b(HVR12_1b<=max_1b)),1));
%x2=2.8+0.4*(rand(length(HVR12_3a(HVR12_3a<=max_3a)),1));
size_marker = 7;
f0=scatter(x0,HVR12_1a(HVR12_1a<=max_1a) ,'o','MarkerEdgeColor','k','MarkerFaceColor',red,'SizeData',size_marker,'LineWidth',0.01);f0.MarkerFaceAlpha = 0.6;hold on
f1=scatter(x1,HVR12_1b(HVR12_1b<=max_1b) ,'o','MarkerEdgeColor','k','MarkerFaceColor',purple,'SizeData',size_marker,'LineWidth',0.01);f1.MarkerFaceAlpha = 0.6;hold on
%f2=scatter(x2,HVR12_3a(HVR12_3a<=max_3a),'o','MarkerEdgeColor','k','MarkerFaceColor',orange,'SizeData',size_marker,'LineWidth',0.01);f2.MarkerFaceAlpha = f1.MarkerFaceAlpha;hold on
ylabel("Escape time")
%xlabel("Subtype")
ylim([0,500])
title("Part of the remaining residues used")

% remaining (没完全算完)
subplot(1,2,2)
data = [NonVRs_1a NonVRs_1b ];
G = [zeros(size(NonVRs_1a))-1 zeros(size(NonVRs_1b)) ];
%figure;
max_1a = prctile(NonVRs_1a,75)+1.5*(prctile(NonVRs_1a,75)-prctile(NonVRs_1a,25));
max_1b = prctile(NonVRs_1b,75)+1.5*(prctile(NonVRs_1b,75)-prctile(NonVRs_1b,25));
%max_3a = prctile(NonVRs_3a,75)+1.5*(prctile(NonVRs_3a,75)-prctile(NonVRs_3a,25));

boxplot(data,G,'whisker',1.5,'labels',label_xaxis_data);hold on
x0=0.8+0.4*(rand(length(NonVRs_1a(NonVRs_1a<=max_1a)),1));
x1=1.8+0.4*(rand(length(NonVRs_1b(NonVRs_1b<=max_1b)),1));
%x2=2.8+0.4*(rand(length(NonVRs_3a(NonVRs_3a<=max_3a)),1));
size_marker = 7;
f0=scatter(x0,NonVRs_1a(NonVRs_1a<=max_1a) ,'o','MarkerEdgeColor','k','MarkerFaceColor',red,'SizeData',size_marker,'LineWidth',0.01);f0.MarkerFaceAlpha = 0.6;hold on
f1=scatter(x1,NonVRs_1b(NonVRs_1b<=max_1b) ,'o','MarkerEdgeColor','k','MarkerFaceColor',purple,'SizeData',size_marker,'LineWidth',0.01);f1.MarkerFaceAlpha = 0.6;hold on
%f2=scatter(x2,NonVRs_3a(NonVRs_3a<=max_3a),'o','MarkerEdgeColor','k','MarkerFaceColor',orange,'SizeData',size_marker,'LineWidth',0.01);f2.MarkerFaceAlpha = f1.MarkerFaceAlpha;hold on
ylabel("Escape time")
%xlabel("Subtype")
ylim([0,500])
title("All of the remaining residues")


%% 2个Hypervariable regions [Isla 2009]
load('escape_time_1a_500.mat')
load('escape_time_HVR1_1b.mat')
load("escape_time_HVR1_3a.mat")
load("escape_time_HVR495_1b.mat")
load("escape_time_HVR495_3a.mat")

run startup.m
label_xaxis_data = {'1a','1b','3a'};

HVR1_residues_1a = 384:410;
HVR495_residues_1a = 495:501;

escape_time_HVR1_1a = mean_escape_time(HVR1_residues_1a-383);
escape_time_HVR495_1a = mean_escape_time(HVR495_residues_1a-383);

% HVR1
subplot(1,2,1)
data = [escape_time_HVR1_1a escape_time_HVR1_1b  escape_time_HVR1_3a ];
G = [zeros(size(escape_time_HVR1_1a))-1 zeros(size(escape_time_HVR1_1b)) ones(size(escape_time_HVR1_3a))];
%figure;
max_1a = prctile(escape_time_HVR1_1a,75)+1.5*(prctile(escape_time_HVR1_1a,75)-prctile(escape_time_HVR1_1a,25));
max_1b = prctile(escape_time_HVR1_1b,75)+1.5*(prctile(escape_time_HVR1_1b,75)-prctile(escape_time_HVR1_1b,25));
max_3a = prctile(escape_time_HVR1_3a,75)+1.5*(prctile(escape_time_HVR1_3a,75)-prctile(escape_time_HVR1_3a,25));

boxplot(data,G,'whisker',1.5,'labels',label_xaxis_data);hold on
x0=0.8+0.4*(rand(length(escape_time_HVR1_1a(escape_time_HVR1_1a<=max_1a)),1));
x1=1.8+0.4*(rand(length(escape_time_HVR1_1b(escape_time_HVR1_1b<=max_1b)),1));
x2=2.8+0.4*(rand(length(escape_time_HVR1_3a(escape_time_HVR1_3a<=max_3a)),1));
size_marker = 15;
f0=scatter(x0,escape_time_HVR1_1a(escape_time_HVR1_1a<=max_1a) ,'o','MarkerEdgeColor','k','MarkerFaceColor',red,'SizeData',size_marker,'LineWidth',0.01);f0.MarkerFaceAlpha = 0.6;hold on
f1=scatter(x1,escape_time_HVR1_1b(escape_time_HVR1_1b<=max_1b) ,'o','MarkerEdgeColor','k','MarkerFaceColor',purple,'SizeData',size_marker,'LineWidth',0.01);f1.MarkerFaceAlpha = 0.6;hold on
f2=scatter(x2,escape_time_HVR1_3a(escape_time_HVR1_3a<=max_3a),'o','MarkerEdgeColor','k','MarkerFaceColor',orange,'SizeData',size_marker,'LineWidth',0.01);f2.MarkerFaceAlpha = f1.MarkerFaceAlpha;hold on
ylabel("Escape time")
xlabel("Subtype")
ylim([0,420])
title("Residues in HVR1")

% HVR495
subplot(1,2,2)
data = [escape_time_HVR495_1a escape_time_HVR495_1b  escape_time_HVR495_3a ];
G = [zeros(size(escape_time_HVR495_1a))-1 zeros(size(escape_time_HVR495_1b)) ones(size(escape_time_HVR495_3a))];
%figure;
max_1a = prctile(escape_time_HVR495_1a,75)+1.5*(prctile(escape_time_HVR495_1a,75)-prctile(escape_time_HVR495_1a,25));
max_1b = prctile(escape_time_HVR495_1b,75)+1.5*(prctile(escape_time_HVR495_1b,75)-prctile(escape_time_HVR495_1b,25));
max_3a = prctile(escape_time_HVR495_3a,75)+1.5*(prctile(escape_time_HVR495_3a,75)-prctile(escape_time_HVR495_3a,25));

boxplot(data,G,'whisker',1.5,'labels',label_xaxis_data);hold on
x0=0.8+0.4*(rand(length(escape_time_HVR495_1a(escape_time_HVR495_1a<=max_1a)),1));
x1=1.8+0.4*(rand(length(escape_time_HVR495_1b(escape_time_HVR495_1b<=max_1b)),1));
x2=2.8+0.4*(rand(length(escape_time_HVR495_3a(escape_time_HVR495_3a<=max_3a)),1));
size_marker = 15;
f0=scatter(x0,escape_time_HVR495_1a(escape_time_HVR495_1a<=max_1a) ,'o','MarkerEdgeColor','k','MarkerFaceColor',red,'SizeData',size_marker,'LineWidth',0.01);f0.MarkerFaceAlpha = 0.6;hold on
f1=scatter(x1,escape_time_HVR495_1b(escape_time_HVR495_1b<=max_1b) ,'o','MarkerEdgeColor','k','MarkerFaceColor',purple,'SizeData',size_marker,'LineWidth',0.01);f1.MarkerFaceAlpha = 0.6;hold on
f2=scatter(x2,escape_time_HVR495_3a(escape_time_HVR495_3a<=max_3a),'o','MarkerEdgeColor','k','MarkerFaceColor',orange,'SizeData',size_marker,'LineWidth',0.01);f2.MarkerFaceAlpha = f1.MarkerFaceAlpha;hold on
ylabel("Escape time")
xlabel("Subtype")
ylim([0,420])
title("Residues in HVR495")

[p,h] = ranksum(escape_time_HVR495_3a(find(escape_time_HVR495_3a<150)),escape_time_HVR495_1b(find(escape_time_HVR495_1b<150)),"tail","left");

%% VRs (HVR1+HVR2+igVR) vs 非VRs
load("VRs_NonVRs.mat")
run startup.m

label_xaxis_data = {'1a','1b','3a'};

% 去掉values of conserved residues (这里把最大的escape time都去掉)
%HVR12_1a = HVR12_1a(HVR12_1a <473);
%HVR12_1b = HVR12_1b(HVR12_1b <441);
%HVR12_3a = HVR12_3a(HVR12_3a <396);

%NonVRs_1a = NonVRs_1a(NonVRs_1a <473);
%NonVRs_1b = NonVRs_1b(NonVRs_1b <441);
%NonVRs_3a = NonVRs_3a(NonVRs_3a <396);


% HVR1+HVR2 (igVR的escape time还没算到）
subplot(1,2,1)
data = [HVR12_1a HVR12_1b  HVR12_3a];
G = [zeros(size(HVR12_1a))-1 zeros(size(HVR12_1b)) ones(size(HVR12_3a))];
%figure;
max_1a = prctile(HVR12_1a,75)+1.5*(prctile(HVR12_1a,75)-prctile(HVR12_1a,25));
max_1b = prctile(HVR12_1b,75)+1.5*(prctile(HVR12_1b,75)-prctile(HVR12_1b,25));
max_3a = prctile(HVR12_3a,75)+1.5*(prctile(HVR12_3a,75)-prctile(HVR12_3a,25));

boxplot(data,G,'whisker',1.5,'labels',label_xaxis_data);hold on
x0=0.8+0.4*(rand(length(HVR12_1a(HVR12_1a<=max_1a)),1));
x1=1.8+0.4*(rand(length(HVR12_1b(HVR12_1b<=max_1b)),1));
x2=2.8+0.4*(rand(length(HVR12_3a(HVR12_3a<=max_3a)),1));
size_marker = 10;
f0=scatter(x0,HVR12_1a(HVR12_1a<=max_1a) ,'o','MarkerEdgeColor','k','MarkerFaceColor',red,'SizeData',size_marker,'LineWidth',0.01);f0.MarkerFaceAlpha = 0.6;hold on
f1=scatter(x1,HVR12_1b(HVR12_1b<=max_1b) ,'o','MarkerEdgeColor','k','MarkerFaceColor',purple,'SizeData',size_marker,'LineWidth',0.01);f1.MarkerFaceAlpha = 0.6;hold on
f2=scatter(x2,HVR12_3a(HVR12_3a<=max_3a),'o','MarkerEdgeColor','k','MarkerFaceColor',orange,'SizeData',size_marker,'LineWidth',0.01);f2.MarkerFaceAlpha = f1.MarkerFaceAlpha;hold on
ylabel("Escape time")
xlabel("Subtype")
ylim([0,500])
title("Residues in variable regions(HVR1+HVR2+[igVR])")

% Non variable regions (目前是101个residues)
subplot(1,2,2)
data = [NonVRs_1a NonVRs_1b  NonVRs_3a ];
G = [zeros(size(NonVRs_1a))-1 zeros(size(NonVRs_1b)) ones(size(NonVRs_3a))];
%figure;
max_1a = prctile(NonVRs_1a,75)+1.5*(prctile(NonVRs_1a,75)-prctile(NonVRs_1a,25));
max_1b = prctile(NonVRs_1b,75)+1.5*(prctile(NonVRs_1b,75)-prctile(NonVRs_1b,25));
max_3a = prctile(NonVRs_3a,75)+1.5*(prctile(NonVRs_3a,75)-prctile(NonVRs_3a,25));

boxplot(data,G,'whisker',1.5,'labels',label_xaxis_data);hold on
x0=0.8+0.4*(rand(length(NonVRs_1a(NonVRs_1a<=max_1a)),1));
x1=1.8+0.4*(rand(length(NonVRs_1b(NonVRs_1b<=max_1b)),1));
x2=2.8+0.4*(rand(length(NonVRs_3a(NonVRs_3a<=max_3a)),1));
size_marker = 10;
f0=scatter(x0,NonVRs_1a(NonVRs_1a<=max_1a) ,'o','MarkerEdgeColor','k','MarkerFaceColor',red,'SizeData',size_marker,'LineWidth',0.01);f0.MarkerFaceAlpha = 0.6;hold on
f1=scatter(x1,NonVRs_1b(NonVRs_1b<=max_1b) ,'o','MarkerEdgeColor','k','MarkerFaceColor',purple,'SizeData',size_marker,'LineWidth',0.01);f1.MarkerFaceAlpha = 0.6;hold on
f2=scatter(x2,NonVRs_3a(NonVRs_3a<=max_3a),'o','MarkerEdgeColor','k','MarkerFaceColor',orange,'SizeData',size_marker,'LineWidth',0.01);f2.MarkerFaceAlpha = f1.MarkerFaceAlpha;hold on
ylabel("Escape time")
xlabel("Subtype")
ylim([0,500])
title("Other residues")

%% 3a VRs vs remaining
clear
load("Isla2009.mat")
run startup.m

label_xaxis_data = {'HVR1+HVR495+HVR575','Remaining'};

data = [VRs_3a NonVRs_3a];
G = [zeros(size(VRs_3a))-1 zeros(size(NonVRs_3a))];
%figure;
max = prctile(VRs_3a,75)+1.5*(prctile(VRs_3a,75)-prctile(VRs_3a,25));
max_3a = prctile(NonVRs_3a,75)+1.5*(prctile(NonVRs_3a,75)-prctile(NonVRs_3a,25));

boxplot(data,G,'whisker',1.5,'labels',label_xaxis_data);hold on
x0=0.8+0.4*(rand(length(VRs_3a(VRs_3a<=max)),1));
x1=1.8+0.4*(rand(length(NonVRs_3a(NonVRs_3a<=max_3a)),1));
size_marker = 10;
f0=scatter(x0,VRs_3a(VRs_3a<=max) ,'o','MarkerEdgeColor','k','MarkerFaceColor',red,'SizeData',size_marker,'LineWidth',0.01);f0.MarkerFaceAlpha = 0.6;hold on
f1=scatter(x1,NonVRs_3a(NonVRs_3a<=max_3a) ,'o','MarkerEdgeColor','k','MarkerFaceColor',purple,'SizeData',size_marker,'LineWidth',0.01);f1.MarkerFaceAlpha = 0.6;hold on
ylabel("Escape time")
xlabel("Type of residues")
ylim([0,420])
title("Escape time comparison for 3a")


%% figure 4 b can be generated by the respective .pse file

%% Box plot of minumum escape times of binding sites of each antibody
load('escape_time_1a_500.mat', 'mean_escape_time')
all_E=mean_escape_time;
HmAbs_all=[];

HmAb_Pierce2016 = [];

%CBH-4D
HmAb_Pierce2016{1} = [494 497 502 504 506:509 511 537 539 540 542:545 547 549:552 554 556 559 561 562 564 565 584 585 592 594 598 600 602 603 607:611 614 617:619 621 623 624 626 627 629:633 638 640 642:644];
%CBH-4G
HmAb_Pierce2016{2} = [494 497 503 505:509 511 517 537 539 540 542:545 547 549:552 554 556 559 561 564 565 584 592:594 597 598 600 602 603 607:611 614 617:619 621 623 624 627 629 631:633 638 640 642:644];
%CBH-4B
HmAb_Pierce2016{3} = [494 497 503 505:509 537 539 540 550:552 554 559 561 564 565 600 602 603 607:611 614 617:619 621 624 627 629 631:633 638 640 642:644];
%CBH-20
% HmAb_Pierce2016{4} = [494 497 503 505:509 537 539 540 542:545 547 549:552 554 559 561 564 565 569 581 583:585 592 594 597 598 600 602 603 607 608 611 614 617:619 621 623 624 627 631 632 638 640 642:644];
%included the residues with close to threshold (RB = 21,22) as well
HmAb_Pierce2016{4} = [494 497 503 505:509 537 539 540 542:545 547 549:552 554 559 561 564 565 569 581 583:585 592 594 597 598 600 602 603 607 608 610 611 614 617:619 621 623 624 627 631 632 638 640 642:644]; 
%CBH-21
HmAb_Pierce2016{5} = [452 494 497 503 505:509 537 539 540 542:545 547 549:552 554 559 561 564 565 569 581 583:585 592:594 597 598 600 602 603 607:611 614 617:619 621 623 624 627 631:633 638 640 642:644];
%CBH-22
HmAb_Pierce2016{6} = [452 494 497 503 505:509 537 539 540 542:545 547 549:552 554 559 561 564 565 569 581 583:585 592:594 597 598 600 602 603 607:611 614 617:619 621 623 624 627 631:633 638 640 642:644];

%HC-1
HmAb_Pierce2016{7} = [429 494 503:506 508 509 529 530 535 537 539 552 554 559 564 607 611 614 617 644];
%HC-11
HmAb_Pierce2016{8} = [425 428 429 436:438 442 443 494 497 502:504 506:509 511 520 530 535 537 539 550:552 554 556 559 564 565 602 603 607 608 611 614 617:621 624 640 643 644];
%A27
HmAb_Pierce2016{9} = [424:429 437 438 494 497 499 502:504 506:509 511 520 529 530 535 537 539 540 550:552 554 556 559 564 565 602 603 607:611 614 616:619 621 623 624 638 640 642:644];

%CBH-23
HmAb_Pierce2016{10} = [494 508 509 537 539 549 552 554 564 611 614 644]; %same as CBH7
%CBH-7
HmAb_Pierce2016{11} = [494 506 508 509 537 539 549 552 554 564 611 614 617 621 644];

%HC84-20
HmAb_Pierce2016{12} = [429 441 494 497 502:509 511 537 539 552 554 559 564 607 608 611 613 614 616:619 621 640 643 644];
%HC84-24
HmAb_Pierce2016{13} = [429 442 443 494 497 502:509 511 537 539 552 554 559 564 607 608 611 614 617:619 621 643 644];
%HC84-26
HmAb_Pierce2016{14} = [441 442 494 497 502:509 511 537 539 552 554 559 564 603 607 608 611 614 617:619 621 640 643 644];

%HC33-1
HmAb_Pierce2016{15} = [413 418 420];
%HC33-4
HmAb_Pierce2016{16} = [408 413 420];
% %CD81_bs
% HmAb_Pierce2016{17} = [420 421 424 427 430 436:438 440:443 523 527 529 530 535 540 613 614 616:618];



%HCV1
HmAbs_Gopal2017{1} = [413 415 418 420 422 624];
%AR1A
HmAbs_Gopal2017{2} = [417 485 494 497 502 504 506:511 514 516 519 537:540 542 544 545 547 549:552 554 559 564 565 602 603 607:611 614 617:619 621 623 624 632 638 640 642 643 644]; %559 is RB=21%
%AR1B
HmAbs_Gopal2017{3} = [494 497 504 506:509 511 537 539 544 545 547:552 554 564 607:608 610 611 614 617:619 621 640 644]; %610 is RB=21%
%AR2A
HmAbs_Gopal2017{4} = [552 607:611 614 617:619 621 623:625 628 638 640 643 644];
%AR3A
HmAbs_Gopal2017{5} = [424 425 428 429 436:438 440:442 485 494 496 497 502:509 511 516 523 525 529 530 535 537 539 552 554 559 564 602 603 607:611 614 616:619 621 623 624 625 640 642:644];
%AR3B
HmAbs_Gopal2017{6} = [424 425 427:429 432 436 437 440:442 494 496 497 502:509 511 516 518 520 523 529 530 535 537 539 540 552 554 555 556 558 559 564 602 603 607:611 614 616:619 621 623 624 625 630 638 640:644];
%AR3C
HmAbs_Gopal2017{7} = [424 425 428:429 437 438 441:443 494 496 497 502:509 511 516 525 529 530 535 537 539 540 552 554 559 564 602 603 607:611 614 616:619 621 623 624 625 630 638 640 641 643 644];
%AR3D
HmAbs_Gopal2017{8} = [424 425 427:430 436:438 440:442 459 494 496 497 499 502:509 511 516 518 520 523 529 530 535 537 539 540 550:552 554 555 556 558 559 564 565 600 602 603 607:612 614 616:619 621 623 624 625 638 640:644];
 
 

%CBH-5
HmAbs_Keck2019{1} = [424:429 436:438 441:443 535 638];
%212.1.1
HmAbs_Keck2019{2} = [425:429 433 434 529 530 535];
%212.1.10
HmAbs_Keck2019{3} = [426 428 429 442];
%212.15
HmAbs_Keck2019{4} = [544 547 549 637 638 639];
%212.25
HmAbs_Keck2019{5} = [549 636 638 639];

 
% Bailey2017

%HEPC3
HmAbs_Bailey2017{1} = [425 427 428 437 499 520 530 535];
%HEPC43
HmAbs_Bailey2017{2} = [425 427 428 432 436 437 438 442 443 499 517 520 527 529 530 535 616];
%HEPC74
HmAbs_Bailey2017{3} = [425 428 436 437 530 535];
%HEPC46
HmAbs_Bailey2017{4} = [541:546 548 549 594 598 633];
%HEPC50
HmAbs_Bailey2017{5} = [543 544 545 549 594 597 598];
%HEPC98
HmAbs_Bailey2017{6} = [402 405 408];

HmAbs_all=[HmAbs_all  HmAbs_Gopal2017 HmAbs_Keck2019 HmAbs_Bailey2017 HmAb_Pierce2016];


for kk = 1:length(HmAbs_all)
    data_mean_HmAb_new{kk} = all_E(HmAbs_all{kk}-383); 
    min_data_mean_HmAb_1a(kk) =min(data_mean_HmAb_new{kk});
end

load('mean_escape_time_1b.mat')
all_E=all_mean_escape_time_1b;

for kk = 1:length(HmAbs_all)
    data_mean_HmAb_new{kk} = all_E(HmAbs_all{kk}-383); 
    min_data_mean_HmAb_1b(kk) =min(data_mean_HmAb_new{kk});
end


data = [min_data_mean_HmAb_1a  min_data_mean_HmAb_1b ];
G = [zeros(size(min_data_mean_HmAb_1a )) ones(size(min_data_mean_HmAb_1b))];



set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultTextFontName','Arial')
set(0,'DefaultAxesFontSize',8)
set(0,'DefaultTextFontSize',8)
box_lineWidth = 0.75;
box_widths_value = 0.5;
black = [0 0 0];
box_color = [black;black];
box_color_transparency = 0; %faceAlpha
median_lineWidth = 1.5;
median_color = 'k';
whisker_value = 1.5;
outlier_marker = '';
outlier_markerSize = 3.5;
outlier_marker_edgeWidth = 0.001;
outlier_marker_edgeColor = 'w';label_xaxis_data = {'1a',sprintf('1b')};
text_ylabel = {'Minumum escape time of',' binding resiudes of HmAbs'};
text_xlabel = '';
text_title = '';%'E2-escape mutations [Keck2009],[Morin2012],[Bailey2015]';
label_orientation_choice = 'horizontal'; %'horizontal'
ylim_min = 0;
ylim_max = 0.2;
savefig = 0;
savefig_name = 'escape_time_compare';
fig_width_cm = 4;
fig_height_cm = 5;
FIG=figure;
set(gcf,'renderer','Painters')

x1=0.8+0.4*(rand(length(min_data_mean_HmAb_1a),1));
x2=1.8+0.4*(rand(length(min_data_mean_HmAb_1b),1));
size_marker = 10;
% green = [27 129 62]/255;
% lightgreen = (1-green)*0.6+green;
f1=scatter(x1,min_data_mean_HmAb_1a ,'o','MarkerEdgeColor','w','MarkerFaceColor',purple,'SizeData',size_marker,'LineWidth',0.01);f1.MarkerFaceAlpha = 0.6;hold on 
f2=scatter(x2,min_data_mean_HmAb_1b,'o','MarkerEdgeColor','w','MarkerFaceColor',orange,'SizeData',size_marker,'LineWidth',0.01);f2.MarkerFaceAlpha = f1.MarkerFaceAlpha;hold on





figure_boxplot(data,G,...
    box_lineWidth,box_widths_value,box_color,box_color_transparency,...
    median_lineWidth,median_color,...
    whisker_value,...
    outlier_marker,outlier_markerSize,outlier_marker_edgeWidth,outlier_marker_edgeColor,outlier_jitter_value,...
    label_xaxis_data,text_ylabel,text_xlabel,text_title,label_orientation_choice,...
    ylim_min,ylim_max,...
    savefig,savefig_name,fig_width_cm,fig_height_cm);

% V=violinplot(data, [repmat("1a",1,length(min_data_mean_HmAb_1a )) repmat("1b",1,length(min_data_mean_HmAb_1b))],'EdgeColor' ,[0 0 0],'BoxColor' ,[0 0 0]);
% SizeData=10;
% ylabel({'Minimum escape time of',' binding residues of HmAbs'});
% V(1, 1).ViolinColor = purple;
% V(1, 2).ViolinColor = orange;
% V(1, 1).EdgeColor = 'None';
% 
% V(1, 1).MedianPlot.SizeData  =SizeData+30;
% V(1, 2).MedianPlot.SizeData  =SizeData+30;
% V(1, 2).EdgeColor = 'None';
% V(1, 1).ScatterPlot.SizeData  =SizeData;
% V(1, 2).ScatterPlot.SizeData  =SizeData;

set(gca,'YTick',0:200:400)
yt = get(gca, 'YTick');
axis([xlim    0  480])
xt = get(gca, 'XTick');
hold on
plot(xt([1 2]), [1 1]*440, '-k','LineWidth',0.3)
% plot(xt([1 1]), [0.95 1]*440, '-k','LineWidth',0.3)
% plot(xt([2 2]), [0.95 1]*440, '-k','LineWidth',0.3)
P2 = ranksum(min_data_mean_HmAb_1b,min_data_mean_HmAb_1a,'tail','left')

ind = floor(log10(P2));
P2 = roundn(P2/10^ind,-1);
t2 = ['$$ P = ',num2str(P2),' \times 10^{',num2str(ind),'} $$'];

% text(0.9,480,t,'interpreter','latex','FontSize',8)


FIG.Name = 'escape_time_compare';


FIG.Units = 'centimeters';
% set(gcf,'Position',[10 10 8 8.58]);
set(gcf,'Position',[6.53 6.53 4 6.27]);
% set(gcf,'Position',[6.53 6.53 5.22 3.92]);
set(gca,'TickLength',[0.035, 0.01])
set(gca,'Position',[.37 .1 .6 .88]);  %调整 XLABLE和YLABLE不会被切掉
figure_FontSize=8;
set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
set(gca,'TickDir','out')
set(gca,'TickLength',[0.02, 0.01])
set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
set(get(gca,'YLabel'), 'Units', 'Normalized', 'Position', [-0.5, 0.5, 0]);
set(findobj('FontSize',10),'FontSize',figure_FontSize);
% try print(['C:\Users\27909\Desktop\' FIG.Name],'-dpng','-r600'); catch      print(['C:\Users\hzhangbr\Desktop\' FIG.Name],'-dpng','-r600');     end
% try print(['C:\Users\27909\Desktop\' FIG.Name],'-dpng','-r600'); catch      print(['C:\Users\hzhangbr\Desktop\' FIG.Name],'-dpdf');     end
