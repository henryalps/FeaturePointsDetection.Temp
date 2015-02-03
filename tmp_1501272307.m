close all
% tmp = bpElbow(1:100000);
% tmp = sigfilter(tmp);
% tmp=tmp(1000:end);
% 
% data=tmp/max(tmp);
data=bpElbow;
len = length(data);
maxLenRet = ceil(len / 300);
peak = zeros(maxLenRet, 2);
valley = peak;
key = peak;
rise = peak;
dicNotch = peak;
dicPeak = peak;
num = 0;

%% ����1������������С������ܴ��ڵķ�Χ
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

%% ����3���������ܵ㣬��������������Ϊ���壺
%       1��ǰ70��1�׵���Ϊ������100��һ�׵���Ϊ����
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

%%����4��ʹ��С���任�����ز����뽵��Ͽ��λ��
%
[dicpeaks,dicnotchs] = tmpFunc(data,peak,[]);