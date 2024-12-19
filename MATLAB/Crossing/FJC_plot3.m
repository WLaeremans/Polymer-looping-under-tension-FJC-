clear;
clc;
close all;
purple = [160/256, 156/256, 252/256];

%% READ DATA

%read data N 19
N19time1 = readtable("N19_f1.txt");
N19time1 = table2array(N19time1);
N19time1 = 0.005.*(N19time1 - 500);
N19t1 = mean(N19time1);

N19time1_25 = readtable("N19_f1_25.txt");
N19time1_25 = table2array(N19time1_25);
N19time1_25 = 0.005.*(N19time1_25 - 500);
N19t1_25 = mean(N19time1_25);

N19time1_5 = readtable("N19_f1_5.txt");
N19time1_5 = table2array(N19time1_5);
N19time1_5 = 0.005.*(N19time1_5 - 500);
N19t1_5 = mean(N19time1_5);

N19time1_75 = readtable("N19_f1_75.txt");
N19time1_75 = table2array(N19time1_75);
N19time1_75 = 0.005.*(N19time1_75 - 500);
N19t1_75 = mean(N19time1_75);

N19time2 = readtable("N19_f2.txt");
N19time2 = table2array(N19time2);
N19time2 = 0.005.*(N19time2 - 500);
N19t2 = mean(N19time2);

times19 = [N19time1_25, N19time1_5, N19time1_75, N19time2];
TT19 = [N19t1_25, N19t1_5, N19t1_75, N19t2];        %f = 1 contained zeros
FF19 = [1.25, 1.5, 1.75, 2];

%read data N 49
N49time05 = readtable("N49_f05.txt");
N49time05 = table2array(N49time05);
N49time05 = 0.005.*(N49time05 - 500);
N49t05 = mean(N49time05);

N49time075 = readtable("N49_f075.txt");
N49time075 = table2array(N49time075);
N49time075 = 0.005.*(N49time075 - 500);
N49t075 = mean(N49time075);

N49time0875 = readtable("N49_f0875.txt");
N49time0875 = table2array(N49time0875);
N49time0875 = 0.005.*(N49time0875 - 500);
N49t0875 = mean(N49time0875);

N49time1 = readtable("N49_f1.txt");
N49time1 = table2array(N49time1);
N49time1 = 0.005.*(N49time1 - 500);
N49t1 = mean(N49time1);

N49time1_125 = readtable("N49_f1_125.txt");
N49time1_125 = table2array(N49time1_125);
N49time1_125 = 0.005.*(N49time1_125 - 500);
N49t1_125 = mean(N49time1_125);

N49time1_25 = readtable("N49_f1_25.txt");
N49time1_25 = table2array(N49time1_25);
N49time1_25 = 0.005.*(N49time1_25 - 500);
N49t1_25 = mean(N49time1_25);

times49 = [N49time075, N49time0875, N49time1, N49time1_125, N49time1_25];
TT49 = [N49t075, N49t0875, N49t1, N49t1_125, N49t1_25];     %t05 contained zeros
FF49 = [0.75, 0.875, 1, 1.125, 1.25];

%read data N 99
N99time045 = readtable("N99_f045.txt");
N99time045 = table2array(N99time045);
N99time045 = 0.005.*(N99time045 - 500);
N99t045 = mean(N99time045);

N99time055 = readtable("N99_f055.txt");
N99time055 = table2array(N99time055);
N99time055 = 0.005.*(N99time055 - 500);
N99t055 = mean(N99time055);

N99time065 = readtable("N99_f065.txt");
N99time065 = table2array(N99time065);
N99time065 = 0.005.*(N99time065 - 500);
N99t065 = mean(N99time065);

N99time075 = readtable("N99_f075.txt");
N99time075 = table2array(N99time075);
N99time075 = 0.005.*(N99time075 - 500);
N99t075 = mean(N99time075);

N99time085 = readtable("N99_f085.txt");
N99time085 = table2array(N99time085);
N99time085 = 0.005.*(N99time085 - 500);
N99t085 = mean(N99time085);

times99 = [N99time045, N99time055, N99time065, N99time075, N99time085];
TT99 = [N99t045, N99t055, N99t065, N99t075, N99t085];
FF99 = [0.45, 0.55, 0.65, 0.75, 0.85];

%% PLOT SCALING WITH 1/P_L
N = 19;
g = @(x) log((sinh(x)/x)^N/(   (sinh(FF19(1))/FF19(1))^N   ) * TT19(1)  );
fplot(g, 'Color', purple, 'LineStyle','-', 'LineWidth', 2)
hold on
SEM = std(times19)./10;
STDs_log19 = SEM./TT19;
errorbar(FF19(1:end),log(TT19(1:end)),STDs_log19,"LineStyle","none", 'Marker','o',"MarkerEdgeColor","k","MarkerFaceColor","k", "MarkerSize", 2, 'LineWidth',2,'LineStyle','none','Color',[0 0 0])
    legend('MD simulation', 'theory', '', 'Position',[0.161183560873863 0.775317475627404 0.38560909428516 0.117380954924084],...
    'Interpreter','latex',...
    'FontSize',14);
N = 49;
g = @(x) log((sinh(x)/x)^N/(   (sinh(FF49(1))/FF49(1))^N   ) * TT49(1)  );
fplot(g, 'Color', purple, 'LineStyle','-', 'LineWidth', 2)
SEM = std(times49)./10;
STDs_log49 = SEM./TT49;
errorbar(FF49(1:end),log(TT49(1:end)),STDs_log49,"LineStyle","none", 'Marker','o',"MarkerEdgeColor","k","MarkerFaceColor","k", "MarkerSize", 2, 'LineWidth',2,'LineStyle','none','Color',[0 0 0])
    legend('MD simulation', 'theory', '', 'Position',[0.161183560873863 0.775317475627404 0.38560909428516 0.117380954924084],...
    'Interpreter','latex',...
    'FontSize',14);
N = 99;
SEM = std(times99)./10;
STDs_log99 = SEM./TT99;
errorbar(FF99(1:end),log(TT99(1:end)),STDs_log99,"LineStyle","none", 'Marker','o',"MarkerEdgeColor","k","MarkerFaceColor","none", "MarkerSize", 2, 'LineWidth',2,'LineStyle','none','Color',[0 0 0])
    legend('MD simulation', 'theory', '', 'Position',[0.161183560873863 0.775317475627404 0.38560909428516 0.117380954924084],...
    'Interpreter','latex',...
    'FontSize',14);
g = @(x) log((sinh(x)/x)^N/(   (sinh(FF99(1))/FF99(1))^N   ) * TT99(1)  );
fplot(g, 'Color', purple, 'LineStyle','-', 'LineWidth', 2)
errorbar(FF99(1:end),log(TT99(1:end)),STDs_log99,"LineStyle","none", 'Marker','o',"MarkerEdgeColor","k","MarkerFaceColor","k", "MarkerSize", 2, 'LineWidth',2,'LineStyle','none','Color',[0 0 0])
    legend('MD simulation', 'theory', '', 'Position',[0.161183560873863 0.775317475627404 0.38560909428516 0.117380954924084],...
    'Interpreter','latex',...
    'FontSize',14);

xlim([0.2 2.4])
ylim([3.6 16.7])
lgnd = legend('', 'Simulation', '', '',  '', 'Theory', '','Location','southeast','Interpreter','latex','FontSize',14);
box on
ax = gca;
ax.FontSize = 12; 
xlabel('$f$', 'Interpreter','latex', 'FontSize',18)
ylabel('$\ln[\tau_{L,x}(N,f)]$', 'Interpreter','latex', 'FontSize',18)
set(lgnd,'color','none');


annotation('textarrow',[0.675357142857143 0.618214285714285],...
    [0.455714285714287 0.495714285714288],'String',{'$N = 19$'},...
    'Interpreter','latex',...
    'FontSize',14);

annotation('textarrow',[0.544642857142856 0.487499999999998],...
    [0.727142857142859 0.767142857142861],'String',{'$N = 49$'},...
    'Interpreter','latex',...
    'FontSize',14);

annotation('textarrow',[0.268928571428571 0.308928571428571],...
    [0.807142857142858 0.767142857142858],'String',{'$N = 99$'},...
    'Interpreter','latex',...
    'FontSize',14);

