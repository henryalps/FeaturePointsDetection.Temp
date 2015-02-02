function [dicNotchs,dicPeaks] = tmpFunc(data,peaks,notchs)
%%预定义
dicNotchs=peaks;
dicPeaks=peaks;
plotOrNot = 1;
complex =0;
calcDxOrNot=1;%作小波变换前是否对原始信号做差分
diffOrNot = 0;%是否对小波系数求差分
scale = 4;
%%求尺度1:6的bior6.8小波变换并绘图
% wl = cwt(data,[1:6],'bior6.8');
% wl = [data(:)';wl];
% subplotNWayFig(wl);
%%作一阶差分
if calcDxOrNot==1
    dataBackup=data;
    data=(lowPassFilter(diff(data)));
end
%%求尺度4的bior1.3小波变换
if complex == 1;
    wl4 = cwt(data,scale,'cgau2');%bior1.3 cgau2
    wl4 = abs(wl4).^2;
else
    wl4 = cwt(data,scale,'bior1.3');%gaus1
end
wl4(1:100) = wl4(1:100)*0;
wl4(end-100:end) = wl4(end-100:end)*0;
if diffOrNot == 1 || complex == 1
    wl4 = diff(wl4);
end
wl4=wl4/(max(abs(wl4)));
%%如果求过一阶差分,那么需要恢复原data
if calcDxOrNot==1
    data = dataBackup;
end
if plotOrNot == 1
    close all
    plot(data);
    hold on
    if calcDxOrNot==1
        plot((lowPassFilter(diff(data))/max(abs(lowPassFilter(diff(data))))),'k');
    end
    plot(wl4,'r');
    legend('脉搏波','脉搏波差分','bior1.5小波3尺度序列');
    title('降中峡与重博波检测');
    line([1,length(data)],[0,0],'color','k');
end
%%求小波变换的过零点
pos = findPassZeroPointPos(wl4);
%%如果没找到过零点，则返回空
if isempty(pos)
    dicNotchs = [];
    dicPeaks = [];
    return
end
%%绘制所有的过零点
if plotOrNot == 1
    for i=1:length(pos)
        line([pos(i),pos(i)],[-1,1],'color','k');
    end
end
%%找到波峰后的第1-4个过零点

function [pos] = findPassZeroPointPos(data)
%%找一组离散数据的近似过零点
%  输入:
%  data: N维离散数据(N>2)
%  输出:
%  pos : 离散数据中近似过零点的位置
% 预处理数据
if min(data)>1 || length(data)<3
    pos = [];
    return
end
%% 求数据移位相乘的结果
plusResult = data(1:end-1).*data(2:end);
%% 找到小于0的数据点的位置
pos = find(double(plusResult<=0));
for i=1:length(pos)
    if abs(data(pos(i))) > abs(data(pos(i) + 1))
        pos(i) = pos(i)+1;
    end
end