function [d,bgvid,bgvid2] = generateBackground(headerinfo,frames,readframe)
sum1 = zeros(headerinfo.max_width,headerinfo.max_height);
sum2 = zeros(headerinfo.max_width,headerinfo.max_height);
% repetitive monte carlo sampling because of cpu memory restrictions
for j = 1:10
    tcount=1;
    v = ceil([rand(1,10).*(frames./2-1) rand(1,10).*(frames./2)+frames./2]); %increments through the frames in sets of 50
    v = sort(v);
    for k=v
        %frames that will be averaged to get the background that will be subtracted
        d.backgroundvideo(:,:,tcount) = readframe(k);
        tcount=tcount+1;
    end
    
    d.background_double = im2double(d.backgroundvideo); % changing format of the video to be able to do math on it
    nframes = size(d.background_double,3);
    d.background_double_before= d.background_double(:,:,1:floor(nframes/2));
    d.background_double_during= d.background_double(:,:,floor(nframes/2+1):nframes);
    sum1 = sum(d.background_double_before,3)+sum1;% sum across all the frames
    sum2 = sum(d.background_double_during,3)+sum2;% sum across all the frames
end

nframes = length(v)*j/2;
d.average = sum1/(nframes);% average frame
d.average1 = sum2/(nframes);% average frame
bgvid = im2uint8(squeeze(d.average)); %?????????
bgvid2 = im2uint8(squeeze(d.average1)); %?????????

end