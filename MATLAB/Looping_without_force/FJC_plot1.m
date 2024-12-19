clear;
clc;
close all;

bbblue = [0, 0.76, 1];
bbpink = [1, 0.42, 1];

%% Read data

%%%%%%%%%%%%%% rc = 1 *************
time14 = readtable("N14_f0.txt");
time14 = table2array(time14);
time14 = 0.005.*(time14 - 500);
t14 = mean(time14);

time19 = readtable("N19_f0.txt");
time19 = table2array(time19);
time19 = 0.005.*(time19 - 500);
t19 = mean(time19);

time29 = readtable("N29_f0.txt");
time29 = table2array(time29);
time29 = 0.005.*(time29 - 500);
t29 = mean(time29);

time39 = readtable("N39_f0.txt");
time39 = table2array(time39);
time39 = 0.005.*(time39 - 500);
t39 = mean(time39);

time49 = readtable("N49_f0.txt");
time49 = table2array(time49);
time49 = 0.005.*(time49 - 500);
t49 = mean(time49);

timesrc1 = [time14, time19, time29, time39, time49];
TTrc1 = [t14, t19, t29, t39, t49];
FF = [14, 19:10:49];

%%%%%%%%%%%%%% rc = 0.5 *************
time14 = readtable("N14_f0_rc05.txt");
time14 = table2array(time14);
time14 = 0.005.*(time14 - 500);
t14 = mean(time14);

time19 = readtable("N19_f0_rc05.txt");
time19 = table2array(time19);
time19 = 0.005.*(time19 - 500);
t19 = mean(time19);

time29 = readtable("N29_f0_rc05.txt");
time29 = table2array(time29);
time29 = 0.005.*(time29 - 500);
t29 = mean(time29);

time39 = readtable("N39_f0_rc05.txt");
time39 = table2array(time39);
time39 = 0.005.*(time39 - 500);
t39 = mean(time39);

time49 = readtable("N49_f0_rc05.txt");
time49 = table2array(time49);
time49 = 0.005.*(time49 - 500);
t49 = mean(time49);

timesrc05 = [time14, time19, time29, time39, time49];
TTrc05 = [t14, t19, t29, t39, t49];

%%%%%%%%%%%%%% rc = 1.5 *************
time14 = readtable("N14_f0_rc1_5.txt");
time14 = table2array(time14);
time14 = 0.005.*(time14);
t14 = mean(time14);

time19 = readtable("N19_f0_rc1_5.txt");
time19 = table2array(time19);
time19 = 0.005.*(time19);
t19 = mean(time19);

time29 = readtable("N29_f0_rc1_5.txt");
time29 = table2array(time29);
time29 = 0.005.*(time29);
t29 = mean(time29);

time39 = readtable("N39_f0_rc1_5.txt");
time39 = table2array(time39);
time39 = 0.005.*(time39);
t39 = mean(time39);

time49 = readtable("N49_f0_rc1_5.txt");
time49 = table2array(time49);
time49 = 0.005.*(time49);
t49 = mean(time49);

timesrc1_5 = [time14, time19, time29, time39, time49];
TTrc1_5 = [t14, t19, t29, t39, t49];

%% ERRORS

%%%%%%%%%%%%%% rc = 0.5 *************
SEM = std(timesrc05)./10;
STDs_log = SEM./TTrc05;

%error bars
errorbar(log(FF(1:end)),log(TTrc05(1:end)),STDs_log,"LineStyle","none", 'Marker','^','LineWidth',2,'LineStyle','none','Color',bbpink)
    legend('MD simulation', 'theory', '', 'Position',[0.161183560873863 0.775317475627404 0.38560909428516 0.117380954924084],...
    'Interpreter','latex',...
    'FontSize',14);
hold on

%linear fit
f = @(x) 1.5*x + 0.8302;
hold on
fplot(f,'MarkerSize',6,'Color',[0.67843137254902 0.67843137254902 1],'LineWidth',2)

%error bars again for legend
errorbar(log(FF(1:end)),log(TTrc05(1:end)),STDs_log,"LineStyle","none", 'Marker','^','LineWidth',2,'LineStyle','none','Color',bbpink)

%%%%%%%%%%%%%% rc = 1 *************
SEM = std(timesrc1)./10;
STDs_log = SEM./TTrc1;

%error bars
errorbar(log(FF(1:end)),log(TTrc1(1:end)),STDs_log,"LineStyle","none", 'Marker','o','LineWidth',2,'LineStyle','none','Color',[0 0 0])
    legend('MD simulation', 'theory', '', 'Position',[0.161183560873863 0.775317475627404 0.38560909428516 0.117380954924084],...
    'Interpreter','latex',...
    'FontSize',14);
hold on

%linear fit
f = @(x) 1.5*x - 0.5357;
hold on
fplot(f,'MarkerSize',6,'Color',[0.67843137254902 0.67843137254902 1],'LineWidth',2)

%error bars again for legend
errorbar(log(FF(1:end)),log(TTrc1(1:end)),STDs_log,"LineStyle","none", 'Marker','o','LineWidth',2,'LineStyle','none','Color',[0 0 0])

%%%%%%%%%%%%%% rc = 1.5 *************
SEM = std(timesrc1_5)./10;
STDs_log = SEM./TTrc1_5;

%error bars
errorbar(log(FF(1:end)),log(TTrc1_5(1:end)),STDs_log,"LineStyle","none", 'Marker','square','LineWidth',2,'LineStyle','none','Color',bbblue)
    legend('MD simulation', 'theory', '', 'Position',[0.161183560873863 0.775317475627404 0.38560909428516 0.117380954924084],...
    'Interpreter','latex',...
    'FontSize',14);
hold on

%linear fit
f = @(x) 1.5*x - 1.3219;
hold on
fplot(f,'MarkerSize',6,'Color',[0.67843137254902 0.67843137254902 1],'LineWidth',2)

%error bars again for legend
errorbar(log(FF(1:end)),log(TTrc1_5(1:end)),STDs_log,"LineStyle","none", 'Marker','square','LineWidth',2,'LineStyle','none','Color',bbblue)

%legend
lgnd = legend('Sim. $r_c = 0.5$', '', '','Sim. $r_c = 1$', '', '', 'Sim. $r_c = 1.5$', 'Theory $\sim N^{3/2}$', 'Location','northwest','Interpreter','latex','FontSize',14);

%general figure settings
box on
ax = gca;
ax.FontSize = 12; 
xlabel('$\ln\left[ N \right]$', 'Interpreter','latex', 'FontSize',18)
ylabel('$\ln\left[ \tau_L(N,f=0) \right]$', 'Interpreter','latex', 'FontSize',18)
set(lgnd,'color','none');
xlim([2.5 4])
ylim([2.3 7.1])
