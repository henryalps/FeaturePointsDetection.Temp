function [dicNotchs,dicPeaks] = tmpFunc(data,peaks,notchs)
%%预定义
dicNotchs=peaks;
dicPeaks=peaks;
plotOrNot = 1;
complex =1;
n=10;%差分阶数
calcDxOrNot=0;%作小波变换前是否对原始信号做差分
diffOrNot = 0;%是否对小波系数求差分
resampleOrNot = 1;%是否减采样
resampleInterval = 10;%减采样周期
scale = 4;
shiftFactor = 0.0;
data = data - shiftFactor;
%% 减采样
if resampleOrNot == 1 && length(data)>resampleInterval
    rData = downsample(data,resampleInterval);
end
%% 先滤波 - 求尺度5/10的meyer小波变换并绘图
% wl = wavedec(rData,3,'dmey');
% wl = [rData(:)';wl];
% subplotNWayFig(wl);
% close all
% plot(wl(2,:)-wl(6,:));
%% 作n偏移差分
if calcDxOrNot==1
    rDataBackup=rData;
    rData =(nShiftDerivation(rData,n));
    if isempty(rData)
        return
    end
end
%% 如果使用复小波，则求高斯二阶导小波。否则求尺度4的bior1.3小波变换
if complex == 1;
    wl4 = cwt(rData,scale,'cgau2');%bior1.3 cgau2
    wl4 = abs(wl4).^2;
else
    wl4 = cwt(rData,scale,'bior1.3');%gaus1
end
wl4(1:100) = wl4(1:100)*0;
wl4(end-100:end) = wl4(end-100:end)*0;
if diffOrNot == 1 || complex == 1
    wl4 = diff(wl4);
end
wl4=wl4/(max(abs(wl4)));
%% 如果求过n偏移差分,那么需要恢复原rData
if calcDxOrNot==1
    swap(rData,rDataBackup);
end
if plotOrNot == 1
    close all
    plot(rData);
    hold on
    if calcDxOrNot==1
        plot(rDataBackup/max(abs(rDataBackup)),'k');
    end
    plot(wl4,'r');
    %legend('脉搏波','脉搏波差分','bior1.5小波3尺度序列');
    title('降中峡与重博波检测');
    line([1,length(rData)],[0,0],'color','k');
end
%% 求小波变换的过零点
pos = findPassZeroPointPos(wl4);
%% 如果没找到过零点，则返回空
if isempty(pos)
    dicNotchs = [];
    dicPeaks = [];
    return
end
%% 绘制所有的过零点
if plotOrNot == 1
    for i=1:length(pos)
        line([pos(i),pos(i)],[-1,1],'color','k');
    end
end
%% 找到波峰后第一个0点的第1-3个小波过零点
% 还原pos到减采样前的data序列位置
if plotOrNot == 1
    close all
    plot(data)
    hold on
end
pos = pos*resampleInterval;
for i=1:length(peaks(:,1))
    % 找到peaks后第一个data与x轴的交点
    pointPos=find(data(peaks(i,1):end)<0,1);
    % 找到该交点后离这个交点最近的3个小波过0点及其对应的位置
    tmp=pos-(peaks(i,1)+pointPos);
    ppos = find(tmp>0,1);
    if ppos<=length(pos)-2
        dicNotchs(i,:) = [pos(ppos),data(pos(ppos))];
        dicPeaks(i,:) = [pos(ppos+2),data(pos(ppos+2))];
    end
    if plotOrNot == 1
        plot(dicNotchs(i,1),dicNotchs(i,2),'ro');
        plot(dicPeaks(i,1),dicPeaks(i,2),'go');
    end
end
    
