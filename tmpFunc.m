function [dicNotchs,dicPeaks] = tmpFunc(data,peaks,notchs)
%%Ԥ����
plotOrNot = 1;
complex =0;
diffOrNot = 0;
scale = 4;
%%��߶�4��bior1.3С���任
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
    legend('������','bior1.5С��3�߶�����');
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
%%���������Լ����еĹ����
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