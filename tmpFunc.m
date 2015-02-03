function [dicNotchs,dicPeaks] = tmpFunc(data,peaks,notchs)
%%Ԥ����
peaks=peaks(peaks(:,1)~=0,:);
dicNotchs=peaks;
dicPeaks=peaks;
plotOrNot = 1;
complex =1;
n=10;%��ֽ���
calcDxOrNot=0;%��С���任ǰ�Ƿ��ԭʼ�ź������
diffOrNot = 0;%�Ƿ��С��ϵ������
resampleOrNot = 1;%�Ƿ������
resampleInterval = 10;%����������
scale = 4;
shiftFactor = 0.1;
data = data - shiftFactor;
%% ������
if resampleOrNot == 1 && length(data)>resampleInterval
    rData = downsample(data,resampleInterval);
end
%% ���˲� - ��߶�5/10��meyerС���任����ͼ
% wl = wavedec(rData,3,'dmey');
% wl = [rData(:)';wl];
% subplotNWayFig(wl);
% close all
% plot(wl(2,:)-wl(6,:));
%% ��nƫ�Ʋ��
if calcDxOrNot==1
    rDataBackup=rData;
    rData =(nShiftDerivation(rData,n));
    if isempty(rData)
        return
    end
end
%% ���ʹ�ø�С���������˹���׵�С����������߶�4��bior1.3С���任
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
%% ������nƫ�Ʋ��,��ô��Ҫ�ָ�ԭrData
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
    %legend('������','���������','bior1.5С��3�߶�����');
    title('����Ͽ���ز������');
    line([1,length(rData)],[0,0],'color','k');
end
%% ��С���任�Ĺ����
pos = findPassZeroPointPos(wl4);
%% ���û�ҵ�����㣬�򷵻ؿ�
if isempty(pos)
    dicNotchs = [];
    dicPeaks = [];
    return
end
%% �������еĹ����
if plotOrNot == 1
    for i=1:length(pos)
        line([pos(i),pos(i)],[-1,1],'color','k');
    end
end
%% �ҵ�������һ��0��ĵ�1-3��С�������
% ��ԭpos��������ǰ��data����λ��
if plotOrNot == 1
    close all
    plot(data)
    hold on
end
pos = pos*resampleInterval;
for i=1:length(peaks(:,1))
    % �ҵ�peaks���һ��data��x��Ľ���
    pointPos=find(data(peaks(i,1):end)<0,1);
    % �ҵ��ý������������������3��С����0�㼰���Ӧ��λ��
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
    
