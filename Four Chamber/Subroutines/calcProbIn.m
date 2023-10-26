function [bin_inside_rim_sliding_SEM,bin_inside_rim_sliding_mean,inside_rim,framesinBin]...
    = calcProbIn(self,lab,border,startPts)
inside_rimTmp = self.(['r' lab])'<border;
time_to_end = size(inside_rimTmp,1)-min(startPts)+1;
inside_rim = nan(time_to_end,self.nFly);
for fly = 1:self.nFly
    tmp = inside_rimTmp(startPts(fly):end,fly);
    inside_rim(1:numel(tmp),fly) = tmp;
end
inside_rim_mean = nanmean(inside_rim,2);
framesinBin = ceil(0.2.*self.fs); %200 ms

%==========================================================================
%sliding bins with 6 frames/bin
bin_inside_rim_sliding = zeros(length(inside_rim_mean),self.nFly);

for fly=1:self.nFly
    for i = 1:length(inside_rim_mean)
        if i< 3 % get the average of frames 1~i
            bin_inside_rim_sliding(i,fly) = mean(inside_rim(1:i,fly));
        elseif i > (length(inside_rim_mean) -3)
            bin_inside_rim_sliding(i,fly) = nanmean(inside_rim(i-2:end,fly));
        else
            bin_inside_rim_sliding(i,fly) = nanmean(inside_rim(i-2:i+3,fly));
        end
    end
end

%mean of binned data
bin_inside_rim_sliding_mean =  nanmean(bin_inside_rim_sliding,2);

%std of binned data
bin_inside_rim_sliding_std = nanstd(bin_inside_rim_sliding,0,2);

%SEM
bin_inside_rim_sliding_SEM = bin_inside_rim_sliding_std./(sqrt(self.nFly));

end