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

%% 步骤1：利用翼函数缩小波峰可能存在的范围
wing50 = wingFunc(data, 50);
thres = mean(wing50(wing50 > 0));
idx = wing50 > thres;
wing10 = wingFunc(data, 10);
idx = idx .* (wing10 > 0);
wing3 = wingFunc(data, 3);
idx = idx .* (wing3 > 0);

%% 步骤2：利用一阶导数缩小范围
data4diff1 = [data(1); data];
d1 = diff(data4diff1);
idx = idx .* (d1 > 0);

%% 步骤3：遍历可能点，满足两个条件则为波峰：
%       1，前70点1阶导数为正，后100点一阶导数为负；
%       2，为前后300点内的最大值
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

%%步骤4：使用小波变换估计重博波与降中峡的位置
%
[dicpeaks,dicnotchs] = tmpFunc(data,peak,[]);