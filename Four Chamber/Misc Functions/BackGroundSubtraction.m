function [] = BackGroundSubtraction()
addpath([pwd '\matlab']);
addpath([pwd '\util']);
clear
s=dir('Videos/*video.ufmf');
filelist={s.name};

for i = 1:length(filelist)
    filename = filelist{i};
    [readframe,frames,~,headerinfo] = get_readframe_fcn(filename);
    
    % create background
    newvidname = ['Videos/' filename(1:end-5) 'bgs.avi'];
    
    clearvars -except i filename readframe frames headerinfo newvidname filelist
    
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
        nframes = size(d.background_double,4);
        d.background_double_before= d.background_double(:,:,1:floor(nframes/2));
        d.background_double_during= d.background_double(:,:,floor(nframes/2+1):nframes);
        sum1 = sum(d.background_double_before,3)+sum1;% sum across all the frames
        sum2 = sum(d.background_double_during,3)+sum2;% sum across all the frames
    end

    nframes = length(v)*j/2;
    d.average = sum1/(nframes);% average frame
    d.average1 = sum2/(nframes);% average frame
    d.average2 = im2uint8(squeeze(d.average)); %?????????
    d.average3 = im2uint8(squeeze(d.average1)); %?????????
    % clear large variables to make room for memory
    d.average = [];%??????
    d.average1 = [];
    d.background_double = [];
    d.background_double_before = [];
    d.background_double_during = [];

    fragmentsize = 25; %sets size of fragments in which parallel processing occurs
    startframe = 1; %sets start
    endframe = startframe+fragmentsize-1; %sets end
    vidObj = VideoWriter(newvidname,'Motion JPEG AVI');%creates new video objects
    open(vidObj); %opens the created video
    bgvid = repmat(d.average2,[1,1,fragmentsize]);%sets average image to the new video
    bgvid2 = repmat(d.average3,[1,1,fragmentsize]);%sets average image to the new video
    
    clear('d.average2');
    while (endframe<=(frames)) %loop from start-->endframe
        for j = 1:fragmentsize
            d.behaviorvideo(:,:,j)=readframe(startframe+j-1);
        end
        if startframe>(frames/2)
            bgvid = bgvid2;
        end
        d.subtra_video(:,:,1,1:fragmentsize) = imsubtract(squeeze(d.behaviorvideo),bgvid); %when object is brighter than background
        d.subtra_video(:,:,2,1:fragmentsize) = imsubtract(squeeze(d.behaviorvideo),bgvid); %when object is brighter than background
        d.subtra_video(:,:,3,1:fragmentsize) = imsubtract(squeeze(d.behaviorvideo),bgvid); %when object is brighter than background
        writeVideo(vidObj,mat2gray(d.subtra_video)); %creates subtracted form of video
        clear('d.subtra_video');
        clear('d.behavior_video');
        startframe = endframe+1;
        if startframe>(frames/2-1)
            bgvid = repmat(d.average3,[1,1,1,fragmentsize]);%sets average image to the new video
        end
        display(startframe);
        endframe = startframe+fragmentsize-1;
    end
    d.behaviorvideo = uint8([]);d.subtra_video = uint8([]);
    if startframe-1 ~= frames%????????????????????
        endframe = frames; %sets endframe = frames
        display('entering end loop'); %displays string
        for j = 1:endframe-startframe+1
            d.behaviorvideo(:,:,j)=readframe(startframe+j-1);
        end
        tmp = bgvid(:,:,1:size(d.behaviorvideo,3));
        endpiece(:,:,1,:) = imsubtract(squeeze(d.behaviorvideo),bgvid(:,:,1:size(d.behaviorvideo,3))); %when object is brighter than background
        endpiece(:,:,1,:) = imsubtract(squeeze(d.behaviorvideo),bgvid(:,:,1:size(d.behaviorvideo,3))); %when object is brighter than background
        endpiece(:,:,1,:) = imsubtract(squeeze(d.behaviorvideo),bgvid(:,:,1:size(d.behaviorvideo,3))); %when object is brighter than background
        writeVideo(vidObj,endpiece);
        close(vidObj);
        clear('d.subtra_video');
    end
    close(vidObj);

end