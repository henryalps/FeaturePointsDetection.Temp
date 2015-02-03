function [dicNotchs,dicPeaks] = tmpFunc(data,peaks,onsets)
%%Ԥ����
peaks=peaks(peaks(:,1)~=0,:);
dicNotchs=peaks;
dicPeaks=peaks;
plotOrNot = 0;
complex =0;
n=10;%��ֽ���
calcDxOrNot=0;%��С���任ǰ�Ƿ��ԭʼ�ź������
diffOrNot = 0;%�Ƿ��С��ϵ������
resampleOrNot = 1;%�Ƿ������
resampleInterval = 10;%����������
scale = 5;
shiftFactor = 0.17;

range = [200,0,0];%�þ��뷨����ز�������ʼ�仯��Χ
numOfZeroPassPoints = zeros(1,length(peaks(:,1))-1);%��ÿ���Ĳ����ڴ��һ��������С���任���������
confactor = 3;%������������

% %% ������
% if resampleOrNot == 1 && length(data)>resampleInterval
%     rData = downsample(data,resampleInterval);
% end
% rData = rData - shiftFactor;
%% ���˲� - ��߶�5/10��meyerС���任����ͼ
% wl = wavedec(rData,3,'dmey');
% wl = [rData(:)';wl];
% subplotNWayFig(wl);
% close all
% plot(wl(2,:)-wl(6,:));
%% ��nƫ�Ʋ��
% if calcDxOrNot==1
%     rDataBackup=rData;
%     rData =(nShiftDerivation(rData,n));
%     if isempty(rData)
%         return
%     end
% end
% %% ���ʹ�ø�С���������˹���׵�С����������߶�4��bior1.3С���任
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
% %% ������nƫ�Ʋ��,��ô��Ҫ�ָ�ԭrData
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
%     %legend('������','���������','bior1.5С��3�߶�����');
%     title('����Ͽ���ز������');
%     line([1,length(rData)],[0,0],'color','k');
% end
% %% ��С���任�Ĺ����
% pos = findPassZeroPointPos(wl4);
% %% ���û�ҵ�����㣬�򷵻ؿ�
% if isempty(pos)
%     dicNotchs = [];
%     dicPeaks = [];
%     return
% end
% %% �������еĹ����
% if plotOrNot == 1
%     for i=1:length(pos)
%         line([pos(i),pos(i)],[-1,1],'color','k');
%     end
% end
%% �ҵ�������һ��0��ĵ�1-3��С�������
% ��ԭpos��������ǰ��data����λ��
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
%     % �ҵ�peaks���һ��data��x��Ľ���
%     pointPos=find(data(peaks(i,1):end)<0,1);
%     % �ҵ��ý������������������3��С����0�㼰���Ӧ��λ��
%     tmp=pos-(peaks(i,1)+pointPos);
%     ppos = find(tmp>0,1);
%     if ppos<=length(pos)-2
%         dicNotchs(i,:) = [pos(ppos),data(pos(ppos))];
%         dicPeaks(i,:) = [pos(ppos+2),data(pos(ppos+2))];
%     end
%   �Ȳ���źţ�����С���任
    dataPart = data(peaks(i,1):peaks(i+1,1));
%   �źŶμ��������һ�� 
    dataPart = downsample(dataPart,resampleInterval);
    dataPart = dataPart/max(abs(dataPart));
%   ��С���任
    wlp = waveletMethodB(dataPart);
%   ��С���任�ĵ�4������㲢����������Ϊ����Ͽλ��
    pos = findPassZeroPointPos(wlp(1:floor(end/2))) ;    
    numOfZeroPassPoints(i) = length(pos);
    pos =  pos(4) * resampleInterval + peaks(i);
%   ������Ͽλ��������������[-resampleInterval,resampleInterval]��Ѱ�Ҷ��׵����ֵ�㣬��Ϊ����Ͽԭʼ����λ��
    if pos + resampleInterval<=length(data) && pos - resampleInterval>0
            [~,shift] = max(abs(diff(data(pos- resampleInterval:pos + resampleInterval),2)));
            pos = pos + shift - resampleInterval; 
    end
%   ��[pos+1:pos+range-1]��Χ��Ѱ�ҵ�����[pos,data(pos)]��[pos+range,data(pos+range)]�������ֱ��
%       �ľ�������ĵ㣬Ȼ������㵽���֮���Ҳ��塣����ҵ��˲��壬����λ����Ϊ�ز���λ�ã�����ʹ���յ�λ��
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
    %% ˥���������µ�range
    range = seqShifter(range,dicPeaks(i,1) - dicNotchs(i,1)+100,i);

%%  ����4���������Ϊ����Ͽ����ͼ
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
%% ��ȥ�����������������������Ĳ�����notchs��peaks
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

    
