function [dicNotchs,dicPeaks] = tmpFunc(data,peaks,onsets)
%%预定义
peaks=peaks(peaks(:,1)~=0,:);
dicNotchs=peaks;
dicPeaks=peaks;
plotOrNot = 0;
complex =0;
n=10;%差分阶数
calcDxOrNot=0;%作小波变换前是否对原始信号做差分
diffOrNot = 0;%是否对小波系数求差分
resampleOrNot = 1;%是否减采样
resampleInterval = 10;%减采样周期
scale = 5;
shiftFactor = 0.17;

range = [200,0,0];%用距离法检测重博波的起始变化范围
numOfZeroPassPoints = zeros(1,length(peaks(:,1))-1);%对每个心搏周期存放一个半周期小波变换过零点数量
confactor = 3;%置信区间因子

% %% 减采样
% if resampleOrNot == 1 && length(data)>resampleInterval
%     rData = downsample(data,resampleInterval);
% end
% rData = rData - shiftFactor;
%% 先滤波 - 求尺度5/10的meyer小波变换并绘图
% wl = wavedec(rData,3,'dmey');
% wl = [rData(:)';wl];
% subplotNWayFig(wl);
% close all
% plot(wl(2,:)-wl(6,:));
%% 作n偏移差分
% if calcDxOrNot==1
%     rDataBackup=rData;
%     rData =(nShiftDerivation(rData,n));
%     if isempty(rData)
%         return
%     end
% end
% %% 如果使用复小波，则求高斯二阶导小波。否则求尺度4的bior1.3小波变换
% if complex == 1;
%     wl4 = cwt(rData,scale,'cgau2');%bior1.3 cgau2
%     wl4 = abs(wl4).^2;
% else
%     wl4 = cwt(rData,scale,'bior1.3');%gaus1
% end
% wl4(1:100) = wl4(1:100)*0;
% wl4(end-100:end) = wl4(end-100:end)*0;
% if diffOrNot == 1 || complex == 1
%     wl4 = diff(wl4);
% end
% wl4=wl4/(max(abs(wl4)));
% %% 如果求过n偏移差分,那么需要恢复原rData
% if calcDxOrNot==1
%     swap(rData,rDataBackup);
% end
% if plotOrNot == 1
%     close all
%     plot(rData);
%     hold on
%     if calcDxOrNot==1
%         plot(rDataBackup/max(abs(rDataBackup)),'k');
%     end
%     plot(wl4,'r');
%     %legend('脉搏波','脉搏波差分','bior1.5小波3尺度序列');
%     title('降中峡与重博波检测');
%     line([1,length(rData)],[0,0],'color','k');
% end
% %% 求小波变换的过零点
% pos = findPassZeroPointPos(wl4);
% %% 如果没找到过零点，则返回空
% if isempty(pos)
%     dicNotchs = [];
%     dicPeaks = [];
%     return
% end
% %% 绘制所有的过零点
% if plotOrNot == 1
%     for i=1:length(pos)
%         line([pos(i),pos(i)],[-1,1],'color','k');
%     end
% end
%% 找到波峰后第一个0点的第1-3个小波过零点
% 还原pos到减采样前的data序列位置
% if plotOrNot == 1
%     close all
%     plot(data)
%     hold on
% end
close all
plot(data)
hold on
% pos = pos*resampleInterval;
for i=1:length(peaks(:,1))-1
%     % 找到peaks后第一个data与x轴的交点
%     pointPos=find(data(peaks(i,1):end)<0,1);
%     % 找到该交点后离这个交点最近的3个小波过0点及其对应的位置
%     tmp=pos-(peaks(i,1)+pointPos);
%     ppos = find(tmp>0,1);
%     if ppos<=length(pos)-2
%         dicNotchs(i,:) = [pos(ppos),data(pos(ppos))];
%         dicPeaks(i,:) = [pos(ppos+2),data(pos(ppos+2))];
%     end
%   先拆分信号，再做小波变换
    dataPart = data(peaks(i,1):peaks(i+1,1));
%   信号段减采样与归一化 
    dataPart = downsample(dataPart,resampleInterval);
    dataPart = dataPart/max(abs(dataPart));
%   求小波变换
    wlp = waveletMethodB(dataPart);
%   求小波变换的第4个过零点并增采样，作为降中峡位置
    pos = findPassZeroPointPos(wlp(1:floor(end/2))) ;    
    numOfZeroPassPoints(i) = length(pos);
    pos =  pos(4) * resampleInterval + peaks(i);
%   将降中峡位置增采样，并在[-resampleInterval,resampleInterval]内寻找二阶导最大值点，作为降中峡原始波形位置
    if pos + resampleInterval<=length(data) && pos - resampleInterval>0
            [~,shift] = max(abs(diff(data(pos- resampleInterval:pos + resampleInterval),2)));
            pos = pos + shift - resampleInterval; 
    end
%   在[pos+1:pos+range-1]范围内寻找到包含[pos,data(pos)]与[pos+range,data(pos+range)]这两点的直线
%       的距离最近的点，然后在起点到这点之间找波峰。如果找到了波峰，则将其位置作为重博波位置，否则使用终点位置
    tmp = data(pos+1:pos+range(1)-1);
    [~,maxpos,~] = poinToLineDistance([pos+1:pos+range(1)-1;tmp(:)']',...
        [pos,data(pos)],[pos+range(1),data(pos+range(1))],1);
    tmp = findPeakShiftInData(data(pos:pos+maxpos));
    dicNotchs(i,:) = [pos,data(pos)];
    if tmp ~= -1
        dicPeaks(i,:) = [pos+tmp-1,data(pos+tmp-1)];
    else
        dicPeaks(i,:) = [pos+maxpos,data(pos+maxpos)];
    end
    %% 衰减法计算新的range
    range = seqShifter(range,dicPeaks(i,1) - dicNotchs(i,1)+100,i);

%%  将第4个过零点作为降中峡并绘图
%     close all
%     plot(data(peaks(i,1):peaks(i+1,1)));
%     hold on
%     hold on
%     plot(wlp,'r');
%     plot(pos(4),dataPart(pos(4)),'go');
    
%     if plotOrNot == 1
%         plot(dicNotchs(i,1),dicNotchs(i,2),'ro');
%         plot(dicPeaks(i,1),dicPeaks(i,2),'go');
%     end
end
%% 除去过零点数量在置信区间外的心搏周期notchs与peaks
    numOfZeroPassPoints = numOfZeroPassPoints - mean(numOfZeroPassPoints);
    var = std(numOfZeroPassPoints);
    numOfZeroPassPoints = abs(numOfZeroPassPoints) < confactor*var;
    dicNotchs = dicNotchs(numOfZeroPassPoints,:);
    dicPeaks = dicPeaks(numOfZeroPassPoints,:);
    for i = 1:length(dicNotchs(:,1))        
        plot(dicNotchs(i,1),dicNotchs(i,2),'ro');
        plot(dicPeaks(i,1),dicPeaks(i,2),'go');
    end
end

    
