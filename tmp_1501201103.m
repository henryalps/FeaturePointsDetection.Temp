close all
tmp = bpWrist(1:100000);
tmp = sigfilter(tmp);
tmp=tmp(1000:end);
dx1 = (diff(tmp));
dx2 = nPointsAverage((diff(tmp,2)),7);

%subplotNWayFig({tmp,dx1,dx2});
figure
hold on,plot(tmp/max(tmp)),plot(dx1/max(dx1),'r'),plot(dx2/max(dx2),'k');

data=tmp/max(tmp);
len = length(data);
maxLenRet = ceil(len / 300);
peak = zeros(maxLenRet, 2);
valley = peak;
key = peak;
rise = peak;
dicNotch = peak;
dicPeak = peak;
num = 0;

%% ����1�������?����С������ܴ��ڵķ�Χ
wing50 = wingFunc(data, 50);
thres = mean(wing50(wing50 > 0));
idx = wing50 > thres;
wing10 = wingFunc(data, 10);
idx = idx .* (wing10 > 0);
wing3 = wingFunc(data, 3);
idx = idx .* (wing3 > 0);

%% ����2������һ�׵�����С��Χ
data4diff1 = [data(1); data];
d1 = diff(data4diff1);
idx = idx .* (d1 > 0);

%% ����3��������ܵ㣬��������������Ϊ���壺
%       1��ǰ70��1�׵���Ϊ���100��һ�׵���Ϊ����
%       2��Ϊǰ��300���ڵ����ֵ
idxs = find(idx);
idxs = idxs(idxs > 500);
idxs = idxs(idxs + 300 < length(data));
for i = 1 : length(idxs)
    idx = idxs(i);
    if diff(data(idx - 50:idx)) >= 0
        if diff(data(idx : idx + 70)) <= 0
            if data(idx) == max(data(idx -300 : idx + 300))
                num = num + 1;
                peak(num, :) = [idx, data(idx)];
            end
        end
    end
end
%% ȥ���ֵ��λ�����зǷ�ֵ
peak=peak(peak>1);
dicnotchs=ones(length(peak),1)*[-1,0];
dicwaves=dicnotchs;
findNotchThed=100;
%% һ�׵������λ��������ȡ
tmp=dx1(1:end-1).*dx1(2:end);
position=tmp<=0;
positions=find(position);
%% ���׵������λ��������ȡ
tmp=dx2(1:end-1).*dx2(2:end);
position1=tmp<=0;
positions1=find(position1);
%% ������׵����������:�������������С����ֵ�ĵ㶪��
tmp=positions1;
ts = 60;
for i=2:length(positions1)-1
    if (positions1(i)-positions1(i-1)<=ts)||(positions1(i+1)-positions1(i)<=ts)
        tmp(i)=-1;
    end
end
positions1save = positions1;
positions1=positions1(tmp~=-1);

%x=(1:length(dx2));
plot(positions1,dx2(positions1),'go');
%% 7/16λ������
posEstimate = floor(diff(peak)*7/16)+peak(1:end-1);
for i=1:length(peak(:,1))-1
    plot(posEstimate(i),data(posEstimate(i)),'bo');
    tmp=positions1 - posEstimate(i);
    [~,pos]=max(tmp(tmp<0)); 
    if (dx2(positions1(pos)+1)>= dx2(positions1(pos))) && (pos>1)
        [~,tmp]=min(dx2(positions1(pos-1):positions1(pos)));
        if isempty(find(dx1(positions1(pos-1):positions1(pos))>0,1))
            tmp=positions1(pos-1)+tmp-1;
            dicwaves(i,:)=[tmp data(tmp)];
            line([tmp,tmp],[-1,1],[0,0],'color','c');

%             METHOD I            
            tmpsave=tmp;
            tmp=positions1save-tmp;
            [~,pos]=max(tmp(tmp<0));
            dicnotchs(i,:)=[positions1save(pos),data(positions1save(pos))];
            tmp=tmpsave;

%             METHOD II
             x=(tmp - findNotchThed:tmp);
            [~, pos] = poinToLineDistance([x(:)'; data(x(:))']',...
                [x(1),data(x(1))],[x(end),data(x(end))]);
            pos = tmp + pos - 1 - findNotchThed;
%             dicnotchs(i,:)=[pos, data(pos)];

%            METHOD III
            if dicnotchs(i,1)>pos
                x=(pos:dicnotchs(i,1));
                [~, tmp] = poinToLineDistance([x(:)'; data(x(:))']',...
                     [x(1),data(x(1))],[x(end),data(x(end))]);
                 pos=pos+tmp+6;
                 dicnotchs(i,:)=[pos,data(pos)];
            end
                
            plot(dicnotchs(i,1), dicnotchs(i,2), 'ko');
            continue;
        end
        pos=tmp + positions1(pos-1) - 1;
    elseif (dx2(positions1(pos)+1)< dx2(positions1(pos))) && (pos<length(positions1))
        [~,tmp]=min(dx2(positions1(pos):positions1(pos+1)));
        if isempty(find(dx1(positions1(pos):positions1(pos+1))>0,1))
            tmp=positions1(pos)+tmp-1;
            dicwaves(i,:)=[tmp data(tmp)];              
            line([tmp,tmp],[-1,1],[0,0],'color','c');           
           
            tmpsave=tmp;
            tmp=positions1save-tmp;
            [~,pos]=max(tmp(tmp<0));
            dicnotchs(i,:)=[positions1save(pos),data(positions1save(pos))];   
            tmp=tmpsave;

            x=(tmp - findNotchThed:tmp);
            [~, pos] = poinToLineDistance([x(:)'; data(x(:))']',...
                [x(1),data(x(1))],[x(end),data(x(end))]);
            pos = tmp + pos -1 - findNotchThed;
%             dicnotchs(i,:)=[pos, data(pos)];

%            METHOD III
            if dicnotchs(i,1)>pos
                x=(pos:dicnotchs(i,1));
                [~, tmp] = poinToLineDistance([x(:)'; data(x(:))']',...
                     [x(1),data(x(1))],[x(end),data(x(end))]);
                 pos=pos+tmp+6;
                 dicnotchs(i,:)=[pos,data(pos)];
            end
            
            plot(dicnotchs(i,1), dicnotchs(i,2), 'ko');   
            
            continue;
        end
        pos=tmp + positions1(pos) - 1;
    else
            continue;
    end
    [~,pos]=min(abs(positions-pos));
    if positions(pos)>posEstimate(i) || positions(pos)<=peak(i,1)
        continue;
    end
    dicwaves(i,:)=[positions(pos) data(positions(pos))];
    plot(dicwaves(i,1), dicwaves(i,2), 'ko');    
    line([positions(pos),positions(pos)],[-1,1],[0,0],'color','r');
    
    posave=pos;
    tmp=positions1save-positions(pos);
    [~,pos]=max(tmp(tmp<0));
    dicnotchs(i,:)=[positions1save(pos),data(positions1save(pos))];
    pos=posave;
    
    tmp=positions(pos);
    x=(tmp - findNotchThed:tmp);
    [~, pos] = poinToLineDistance([x(:)'; data(x(:))']',...
        [x(1),data(x(1))],[x(end),data(x(end))]);
    pos = tmp + pos -1 - findNotchThed;
%     dicnotchs(i,:)=[pos, data(pos)]; 

%            METHOD III
            if dicnotchs(i,1)>pos
                x=(pos:dicnotchs(i,1));
                [~, tmp] = poinToLineDistance([x(:)'; data(x(:))']',...
                     [x(1),data(x(1))],[x(end),data(x(end))]);
                 pos=pos+tmp+6;
                 dicnotchs(i,:)=[pos,data(pos)];
            end
    
    plot(dicnotchs(i,1), dicnotchs(i,2), 'ko');   
end