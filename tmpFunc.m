function [dicNotchs,dicPeaks] = tmpFunc(data,peaks,notchs)
%%预定义
plotOrNot = 1;
complex =0;
diffOrNot = 0;
scale = 4;
%%求尺度4的bior1.3小波变换
if complex == 1;
    wl4 = cwt(data,scale,'bior1.5');%bior1.3 cgau2
    wl4 = abs(wl4).^2;
else
    wl4 = cwt(data,scale,'gaus1');
end
wl4(1:100) = wl4(1:100)*0;
wl4(end-100:end) = wl4(end-100:end)*0;
if diffOrNot == 1 || complex == 1
    wl4 = diff(wl4);
end
wl4=wl4/(abs(max(wl4)));
if plotOrNot == 1
    close all
    plot(data);
    hold on
    plot(wl4,'r');
    legend('脉搏波','bior1.5小波3尺度序列');
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
%%绘制数据以及所有的过零点
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