close all; clear all; clc
days=[1,2,3,4,5];
temperature=[0.9074,0.9179,0.9244,0.9189,0.9123];
concentration=[0,2,6,4,2];
fun=@(x,y) bar(x,y,0.4);
[hAxes,hBar,hLine]=plotyy(days,temperature,days,concentration,fun,'plot');

% for i = 1:length(days)
% text(days(i),temperature(i)+0.5,num2str(temperature));
% end 

grid on;
set(gca,'XTickLabel',{'5×5×5','5×5×5+3×3×3','5×5×5+5×5×5','5×5×5+7×7×7','5×5×5+9×9×9'}); 
set(hLine,'color',[0.4,0.4,0.4],'LineWidth',2,'Marker','.','MarkerSize',25,...
    'MarkerFace','r');
    title('')
xlabel('卷积核组规模')
ylabel(hAxes(1),'平均准确率')
ylabel(hAxes(2),'最优次数')

set(hBar,'Facecolor',[0.5,0.7,1],'EdgeColor',[1,1,1])
set(gca,'Ylim',[0.90 0.93],'ytick',[0.90:0.01:0.93])
