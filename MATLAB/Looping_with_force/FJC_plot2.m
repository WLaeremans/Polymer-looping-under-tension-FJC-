clear;
clc;
close all;

bbblue = [0, 0.76, 1];
bbpink = [1, 0.42, 1];

%% read data N19
N = 19;

time0 = readtable("N19_f0.txt");
time0 = table2array(time0);
time0 = 0.005.*(time0 - 500);
t0 = mean(time0);

time05 = readtable("N19_f05.txt");
time05 = table2array(time05);
time05 = 0.005.*(time05 - 500);
t05 = mean(time05);

time1 = readtable("N19_f1.txt");
time1 = table2array(time1);
time1 = 0.005.*(time1 - 500);
t1 = mean(time1);

time1_5 = readtable("N19_f1_5.txt");
time1_5 = table2array(time1_5);
time1_5 = 0.005.*(time1_5 - 500);
t1_5 = mean(time1_5);

time2 = readtable("N19_f2.txt");
time2 = table2array(time2);
time2 = 0.005.*(time2 - 500);
t2 = mean(time2);

time3 = readtable("N19_f3.txt");
time3 = table2array(time3);
time3 = 0.005.*(time3 - 500);
t3 = mean(time3);

% rc = 0.5
N = 19;

time0rc05 = readtable("N19f0_rc05.txt");
time0rc05 = table2array(time0rc05);
time0rc05 = 0.005.*(time0rc05);
t0rc05 = mean(time0rc05);

time05rc05 = readtable("N19f05_rc05.txt");
time05rc05 = table2array(time05rc05);
time05rc05 = 0.005.*(time05rc05);
t05rc05 = mean(time05rc05);

time1rc05 = readtable("N19f1_rc05.txt");
time1rc05 = table2array(time1rc05);
time1rc05 = 0.005.*(time1rc05);
t1rc05 = mean(time1rc05);

time1_5rc05 = readtable("N19f1_5_rc05.txt");
time1_5rc05 = table2array(time1_5rc05);
time1_5rc05 = 0.005.*(time1_5rc05);
t1_5rc05 = mean(time1_5rc05);

time2rc05 = readtable("N19f2_rc05.txt");
time2rc05 = table2array(time2rc05);
time2rc05 = 0.005.*(time2rc05);
t2rc05 = mean(time2rc05);

TTrc05 = [t0rc05, t05rc05, t1rc05, t1_5rc05, t2rc05];
FFrc05 = [0:0.5:2];

timesrc05 = [time0rc05, time05rc05, time1rc05, time1_5rc05, time2rc05];

% rc = 1.5
N = 19;

time0rc1_5 = readtable("N19f0_rc1_5.txt");
time0rc1_5 = table2array(time0rc1_5);
time0rc1_5 = 0.005.*(time0rc1_5);
t0rc1_5 = mean(time0rc1_5);

time075rc1_5 = readtable("N19f075_rc1_5.txt");
time075rc1_5 = table2array(time075rc1_5);
time075rc1_5 = 0.005.*(time075rc1_5);
t075rc1_5 = mean(time075rc1_5);

time1_25rc1_5 = readtable("N19f1_25_rc1_5.txt");
time1_25rc1_5 = table2array(time1_25rc1_5);
time1_25rc1_5 = 0.005.*(time1_25rc1_5);
t1_25rc1_5 = mean(time1_25rc1_5);

time1_75rc1_5 = readtable("N19f1_75_rc1_5.txt");
time1_75rc1_5 = table2array(time1_75rc1_5);
time1_75rc1_5 = 0.005.*(time1_75rc1_5);
t1_75rc1_5 = mean(time1_75rc1_5);

TTrc1_5 = [t0rc1_5, t075rc1_5, t1_25rc1_5, t1_75rc1_5];
FFrc1_5 = [0, 0.75, 1.25, 1.75];

timesrc1_5 = [time0rc1_5, time075rc1_5, time1_25rc1_5, time1_75rc1_5];

SEMrc1_5 = std(timesrc1_5)./10;
STDs_logrc1_5 = SEMrc1_5./TTrc1_5;

%% ERROR N = 19
times = [time0, time05, time1, time1_5, time2];
TT = [t0, t05, t1, t1_5, t2];
FF = [0:0.5:2];

SEM = std(times)./10;
STDs_log = SEM./TT;

SEMrc05 = std(timesrc05)./10;
STDs_logrc05 = SEMrc05./TTrc05;


%% Plot N = 19

errorbar(FFrc05(1:end),log(TTrc05(1:end)./TTrc05(1)),STDs_logrc05,"LineStyle","none", 'Marker','^','LineWidth',2,'LineStyle','none','Color',bbpink)
hold on
errorbar(FF(1:end),log(TT(1:end)./TT(1)),STDs_log,"LineStyle","none", 'Marker','o','LineWidth',2,'LineStyle','none','Color',[0 0 0])
errorbar(FFrc1_5(1:end),log(TTrc1_5(1:end)./TTrc1_5(1)),STDs_logrc1_5,"LineStyle","none", 'Marker','square', 'MarkerSize', 8,'LineWidth',2,'LineStyle','none','Color',bbblue);
Theory = @(f) N*log(sinh(f)/f);
fplot(Theory,'MarkerSize',6,'Color',[0.67843137254902 0.67843137254902 1],'LineWidth',2)

% exact integral
% Define parameters
N = 19; % Example value for N
rc = 1; % Example value for rc
alpha = 2*N/3; % Calculate alpha

% Define the function P(f)
extra = @(f) exp(-((rc * (rc + f * alpha)) / alpha)) .* (-2 * (-1 + exp(2 * f * rc)) * alpha + exp((2 * rc + f * alpha).^2 / (4 * alpha)) .* f .* sqrt(pi) .* alpha^(3/2) .* ( erf((2 * rc - f * alpha) ./ (2 * sqrt(alpha))) + erf((2 * rc + f * alpha) ./ (2 * sqrt(alpha))) ) ) ;
P_f = @(f)  (sinh(f) ./ f).^(-N) .* ( (1 ./ (4 * f)) .* exp(-((rc * (rc + f * alpha)) / alpha)) .* (-2 * (-1 + exp(2 * f * rc)) * alpha + exp((2 * rc + f * alpha).^2 / (4 * alpha)) .* f .* sqrt(pi) .* alpha^(3/2) .* ( erf((2 * rc - f * alpha) ./ (2 * sqrt(alpha))) + erf((2 * rc + f * alpha) ./ (2 * sqrt(alpha))) ) ) );

% Define the range of f
f_values = linspace(0.01, 2, 100); % Avoid f=0 to prevent division by zero

% Compute P(0)
P_0 = P_f(0.01); % Approximate P(0) using a very small f to avoid division by zero

% Compute the ratio and its logarithm
log_ratio = log(P_0 ./ P_f(f_values));

% Plot
plot(f_values, log_ratio,'LineWidth',1,'LineStyle', '--', 'Color', 'k');

errorbar(FFrc05(1:end),log(TTrc05(1:end)./TTrc05(1)),STDs_logrc05,"LineStyle","none", 'Marker','^','LineWidth',2,'LineStyle','none','Color',bbpink)

errorbar(FF(1:end),log(TT(1:end)./TT(1)),STDs_log,"LineStyle","none", 'Marker','o','LineWidth',2,'LineStyle','none','Color',[0 0 0])

errorbar(FFrc1_5(1:end),log(TTrc1_5(1:end)./TTrc1_5(1)),STDs_logrc1_5,"LineStyle","none", 'Marker','square', 'MarkerSize', 8,'LineWidth',2,'LineStyle','none','Color',bbblue);

%% read N = 49
N = 49;

time0 = readtable("N49_f0.txt");
time0 = table2array(time0);
time0 = 0.005.*(time0);
t0 = mean(time0);

time05 = readtable("N49_f05.txt");
time05 = table2array(time05);
time05 = 0.005.*(time05);
t05 = mean(time05);

time075 = readtable("N49_f075.txt");
time075 = table2array(time075);
time075 = 0.005.*(time075);
t075 = mean(time075);

time1 = readtable("N49_f1.txt");
time1 = table2array(time1);
time1 = 0.005.*(time1);
t1 = mean(time1);

times = [time0, time05, time075, time1];
TT = [t0, t05, t075, t1];
FF = [0, 0.5, 0.75, 1];

SEM = std(times)./10;
STDs_log = SEM./TT;


time0rc05 = readtable("N49f0_rc05.txt");
time0rc05 = table2array(time0rc05);
time0rc05 = 0.005.*(time0rc05);
t0rc05 = mean(time0rc05);

time05rc05 = readtable("N49f05_rc05.txt");
time05rc05 = table2array(time05rc05);
time05rc05 = 0.005.*(time05rc05);
t05rc05 = mean(time05rc05);

time075rc05 = readtable("N49f075_rc05.txt");
time075rc05 = table2array(time075rc05);
time075rc05 = 0.005.*(time075rc05);
t075rc05 = mean(time075rc05);

time1rc05 = readtable("N49f1_rc05.txt");
time1rc05 = table2array(time1rc05);
time1rc05 = 0.005.*(time1rc05);
t1rc05 = mean(time1rc05);

TTrc05 = [t0rc05, t05rc05, t075rc05, t1rc05];
FFrc05 = [0, 0.5, 0.75, 1];

timesrc05 = [time0rc05, time05rc05, time075rc05, time1rc05];

SEMrc05 = std(timesrc05)./10;
STDs_logrc05 = SEMrc05./TTrc05;

time0rc1_5 = readtable("N49f0_rc1_5.txt");
time0rc1_5 = table2array(time0rc1_5);
time0rc1_5 = 0.005.*(time0rc1_5);
t0rc1_5 = mean(time0rc1_5);

time064rc1_5 = readtable("N49f064_rc1_5.txt");
time064rc1_5 = table2array(time064rc1_5);
time064rc1_5 = 0.005.*(time064rc1_5);
t064rc1_5 = mean(time064rc1_5);

time09rc1_5 = readtable("N49f09_rc1_5.txt");
time09rc1_5 = table2array(time09rc1_5);
time09rc1_5 = 0.005.*(time09rc1_5);
t09rc1_5 = mean(time09rc1_5);

time117rc1_5 = readtable("N49f117_rc1_5.txt");
time117rc1_5 = table2array(time117rc1_5);
time117rc1_5 = 0.005.*(time117rc1_5);
t117rc1_5 = mean(time117rc1_5);

TTrc1_5 = [t0rc1_5, t064rc1_5, t09rc1_5, t117rc1_5];
FFrc1_5 = [0, 0.64, 0.9, 1.17];

timesrc1_5 = [time0rc1_5, time064rc1_5, time09rc1_5, time117rc1_5];

SEMrc1_5 = std(timesrc1_5)./10;
STDs_logrc1_5 = SEMrc1_5./TTrc1_5;

%% Plot N = 49

errorbar(FFrc05(1:end),log(TTrc05(1:end)./TTrc05(1)),STDs_logrc05,"LineStyle","none", 'Marker','^','LineWidth',2,'LineStyle','none','Color',bbpink)
hold on
errorbar(FF(1:end),log(TT(1:end)./TT(1)),STDs_log,"LineStyle","none", 'Marker','o','LineWidth',2,'LineStyle','none','Color',[0 0 0])
errorbar(FFrc1_5(1:end),log(TTrc1_5(1:end)./TTrc1_5(1)),STDs_logrc1_5,"LineStyle","none", 'Marker','square', 'MarkerSize', 8,'LineWidth',2,'LineStyle','none','Color',bbblue);

hold on
Theory = @(f) N*log(sinh(f)/f);
fplot(Theory,'MarkerSize',6,'Color',[0.67843137254902 0.67843137254902 1],'LineWidth',2)

% exact integral
% Define parameters
N = 49; % Example value for N
rc = 1; % Example value for rc
alpha = 2*N/3; % Calculate alpha

% Define the function P(f)
extra = @(f) exp(-((rc * (rc + f * alpha)) / alpha)) .* (-2 * (-1 + exp(2 * f * rc)) * alpha + exp((2 * rc + f * alpha).^2 / (4 * alpha)) .* f .* sqrt(pi) .* alpha^(3/2) .* ( erf((2 * rc - f * alpha) ./ (2 * sqrt(alpha))) + erf((2 * rc + f * alpha) ./ (2 * sqrt(alpha))) ) ) ;
P_f = @(f)  (sinh(f) ./ f).^(-N) .* ( (1 ./ (4 * f)) .* exp(-((rc * (rc + f * alpha)) / alpha)) .* (-2 * (-1 + exp(2 * f * rc)) * alpha + exp((2 * rc + f * alpha).^2 / (4 * alpha)) .* f .* sqrt(pi) .* alpha^(3/2) .* ( erf((2 * rc - f * alpha) ./ (2 * sqrt(alpha))) + erf((2 * rc + f * alpha) ./ (2 * sqrt(alpha))) ) ) );

% Define the range of f
f_values = linspace(0.01, 2, 100); % Avoid f=0 to prevent division by zero

% Compute P(0)
P_0 = P_f(0.01); % Approximate P(0) using a very small f to avoid division by zero

% Compute the ratio and its logarithm
log_ratio = log(P_0 ./ P_f(f_values));

% Plot
plot(f_values, log_ratio,'LineWidth',1,'LineStyle', '--', 'Color', 'k');
errorbar(FFrc05(1:end),log(TTrc05(1:end)./TTrc05(1)),STDs_logrc05,"LineStyle","none", 'Marker','^','LineWidth',2,'LineStyle','none','Color',bbpink)
errorbar(FF(1:end),log(TT(1:end)./TT(1)),STDs_log,"LineStyle","none", 'Marker','o','LineWidth',2,'LineStyle','none','Color',[0 0 0])
errorbar(FFrc1_5(1:end),log(TTrc1_5(1:end)./TTrc1_5(1)),STDs_logrc1_5,"LineStyle","none", 'Marker','square', 'MarkerSize', 8,'LineWidth',2,'LineStyle','none','Color',bbblue);


annotation('textarrow',[0.640357142857143 0.583214285714285],...
    [0.744285714285714 0.784285714285716],'String',{'$N = 49$'},...
    'Interpreter','latex',...
    'FontSize',14);

annotation('textarrow',[0.679642857142857 0.622499999999999],...
    [0.394761904761906 0.434761904761907],'String',{'$N = 19$'},...
    'Interpreter','latex',...
    'FontSize',14);

lgnd = legend('Sim. $r_c = 0.5$', 'Sim. $r_c = 1$', 'Sim. $r_c = 1.5$', 'Thr. appr.', 'Thr. ex. $r_c = 1$', 'Location','northwest','Interpreter','latex','FontSize',14);
box on
ax = gca;
ax.FontSize = 12; 
xlabel('$\beta f b$', 'Interpreter','latex', 'FontSize',18)
ylabel('$\ln\left[ \tau_L(N,f)/\tau_L(N,f=0) \right]$', 'Interpreter','latex', 'FontSize',18)
set(lgnd,'color','none');
xlim([0 2])
ylim([-0.5 12.1])


