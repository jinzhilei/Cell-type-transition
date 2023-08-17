clear; clc;

FIG=figure(1);
outpath='output/figdat/';


colors=[0 1 0; 0  0  1;0  1  1; 0.3  0.3  0.3];

A=load(strcat(outpath,'sys.dat'));

%1-4th rows,Mean cell number
%5-8th rows, Cell number standard deviation
%9-11th rows£¬The average ratio of the number of cells
%12-14th rows£¬The standard deviation of the ratio of cell numbers

subplot(3,2,1)

x=1:5;

h1=plot(x,A(2,:),'-','LineWidth',1,'Color',colors(1,:));
hold on;
h2=plot(x,A(3,:),'-','LineWidth',1,'Color',colors(2,:));
hold on;
plot(x,A(4,:),'-','LineWidth',1,'Color',colors(3,:))
hold on;
plot(x,A(1,:),'-','LineWidth',1,'Color',colors(4,:)) 
hold on;

ylim([0 12000])

set(gca,'Ytick',[0 2000 4000 6000 8000 10000 12000])
set(gca,'YTickLabel',{'0','2000','4000','6000','8000','10000','12000'});

errorbar(x,A(2,:),A(6,:),'-','LineWidth',1,'Color',colors(1,:));
hold on;
errorbar(x,A(3,:),A(7,:),'-','LineWidth',1,'Color',colors(2,:));
hold on;
errorbar(x,A(4,:),A(8,:),'-','LineWidth',1,'Color',colors(3,:));
hold on;
errorbar(x,A(1,:),A(5,:),':','LineWidth',1,'Color',colors(4,:))
hold off;

h=legend('SC','TA1','TA2','Total','Location','north');

h.Box='off';
set(gca,'Xtick',[1 2 3 4 5])
set(gca,'XTickLabel',{'0.800','0.825','0.850','0.975','0.900'});


h.NumColumns = 2;

xlabel('m_2');

ylabel('Cell Number');

title('A');

hold off

subplot(3,2,2)

r1=[9 10 11];
r2=[12 13 14];
groupnames={'0.800','0.825','0.850','0.875','0.900'};
Title = 'B';
Xlabel = 'm_2';
Ylabel='Cell Ratio(%)';
barweb(A(r1,:)',A(r2,:)',1,groupnames,Title,Xlabel,Ylabel,jet,'none',[],2,'plot');
ylim([0 85])
% 
h=legend('SC','TA1','TA2','Location','north');
h.NumColumns = 3;
h.Box='off';




A=load(strcat(outpath,'tran.dat'));


%1-9th rows represent the mean transition probability
%10-18th rows represents the standard deviation of the transition probability

% 1-th row 1:SC-SC; 5-th row : TA1-TA1;  9-th row: TA2-TA2;
% 2-th row :SC-TA1; 3-th row: SC-TA2;
% 4-th row  :TA1-SC;  7-th row: TA2-SC;
% 6-th row :TA1-TA2;  8-th row: TA2-TA1;


subplot(3,2,3)

r=[1 5 9];
groupnames={'0.800','0.825','0.850','0.875','0.900'};
Title = 'C';
Xlabel = 'm_2';
Ylabel='Transition Probability(%)';
barweb(A(r,:)',A(r+9,:)',1,groupnames,Title,Xlabel,Ylabel,jet,'none',[],2,'plot');

ylim([0 160])
h=legend('SC-SC','TA1-TA1','TA2-TA2','Location','north');
h.NumColumns = 2;
h.Box='off';

subplot(3,2,4)

r=[2 3];
groupnames={'0.800','0.825','0.850','0.875','0.900'};
Title = 'D';
Xlabel = 'm_2';
Ylabel='Transition Probability(%)';
barweb(A(r,:)',A(r+9,:)',1,groupnames,Title,Xlabel,Ylabel,jet,'none',[],2,'plot');

ylim([0 20])
h=legend('SC-TA1','SC-TA2','Location','north');

h.NumColumns = 3;
h.Box='off';


subplot(3,2,5)

r=[4 7];
groupnames={'0.800','0.825','0.850','0.875','0.900'};
Title = 'E';
Xlabel = 'm_2';
Ylabel='Transition Probability(%)';
barweb(A(r,:)',A(r+9,:)',1,groupnames,Title,Xlabel,Ylabel,jet,'none',[],2,'plot');



ylim([0 2])
set(gca,'Ytick',[0 0.5 1 1.5 2])
set(gca,'YTickLabel',{'0','0.5','1.0','1.5','2.0'});

h=legend('TA1-SC','TA2-SC','Location','north');

h.NumColumns = 3;

h.Box='off';


subplot(3,2,6)

r=[6 8];
groupnames={'0.800','0.825','0.850','0.875','0.900'};
Title = 'F';
Xlabel = 'm_2';
Ylabel='Transition Probability(%)';
barweb(A(r,:)',A(r+9,:)',1,groupnames,Title,Xlabel,Ylabel,jet,'none',[],2,'plot');

ylim([0 0.15])
set(gca,'Ytick',[0 0.05 0.10 0.15])
set(gca,'YTickLabel',{'0','0.05','0.10','0.15'});
h=legend('TA1-TA2','TA2-TA1','Location','north');

h.NumColumns = 3;

h.Box='off';

exportfig(FIG,'output/fig/res.eps', 'FontMode' ,'fixed','FontSize',8, 'color', 'cmyk' );

clear; clc;