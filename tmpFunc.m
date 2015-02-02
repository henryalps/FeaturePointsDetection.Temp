function [dicNotchs,dicPeaks] = tmpFunc(data,peaks,notchs)
%%Ԥ����
dicNotchs=peaks;
dicPeaks=peaks;
plotOrNot = 1;
complex =0;
calcDxOrNot=1;%��С���任ǰ�Ƿ��ԭʼ�ź������
diffOrNot = 0;%�Ƿ��С��ϵ������
scale = 4;
%%��߶�1:6��bior6.8С���任����ͼ
% wl = cwt(data,[1:6],'bior6.8');
% wl = [data(:)';wl];
% subplotNWayFig(wl);
%%��һ�ײ��
if calcDxOrNot==1
    dataBackup=data;
    data=(lowPassFilter(diff(data)));
end
%%��߶�4��bior1.3С���任
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
%%������һ�ײ��,��ô��Ҫ�ָ�ԭdata
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
    legend('������','���������','bior1.5С��3�߶�����');
    title('����Ͽ���ز������');
    line([1,length(data)],[0,0],'color','k');
end
%%��С���任�Ĺ����
pos = findPassZeroPointPos(wl4);
%%���û�ҵ�����㣬�򷵻ؿ�
if isempty(pos)
    dicNotchs = [];
    dicPeaks = [];
    return
end
%%�������еĹ����
if plotOrNot == 1
    for i=1:length(pos)
        line([pos(i),pos(i)],[-1,1],'color','k');
    end
end
%%�ҵ������ĵ�1-4�������

function [pos] = findPassZeroPointPos(data)
%%��һ����ɢ���ݵĽ��ƹ����
%  ����:
%  data: Nά��ɢ����(N>2)
%  ���:
%  pos : ��ɢ�����н��ƹ�����λ��
% Ԥ��������
if min(data)>1 || length(data)<3
    pos = [];
    return
end
%% ��������λ��˵Ľ��
plusResult = data(1:end-1).*data(2:end);
%% �ҵ�С��0�����ݵ��λ��
pos = find(double(plusResult<=0));
for i=1:length(pos)
    if abs(data(pos(i))) > abs(data(pos(i) + 1))
        pos(i) = pos(i)+1;
    end
end