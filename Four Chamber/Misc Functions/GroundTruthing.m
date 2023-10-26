function [] = GroundTruthing()
clear
addpath(genpath([fileparts(which(mfilename)) '\matlab']));
addpath(genpath([fileparts(which(mfilename)) '\util']));

s=dir('Videos/*video.ufmf');
s2=dir('Data Files/*.mat');
filelist={s2.name};
fName = cell(length(filelist),1);
ndx = zeros(length(filelist),1);
for i = 1:length(filelist)
    C = strsplit(filelist{i},'.');
    fName{i} = C{1}(1:end-1);
    ndx(i,1) = str2double(C{1}(end));
end

[videos,~,ndx(:,2)] = unique(cellfun(@num2str,fName,'uni',0));

nFrames = 20;
nFiles = length(filelist);
us.xH = zeros(nFiles,nFrames);ip.xH = zeros(nFiles,nFrames);
us.yH = zeros(nFiles,nFrames);ip.yH = zeros(nFiles,nFrames);
us.x = zeros(nFiles,nFrames);ip.x = zeros(nFiles,nFrames);
us.y = zeros(nFiles,nFrames);ip.y = zeros(nFiles,nFrames);
us.ang = zeros(nFiles,nFrames);ip.ang = zeros(nFiles,nFrames);

for i = 2:length(videos)
    vidName = ['Videos/' videos{i} 'video.ufmf'];
    
    [readframe,~,~,~] = get_readframe_fcn(vidName);
    load([vidName(1:end-5) '_Inside.mat'],'Inside','vertices')
    load([vidName(1:end-5) '_Inside.mat'],'Inside','vertices')
    
    validArenas = ndx(ndx(:,2)==i,1);
    
    v.randFrames = randi(10700,nFrames,1)+99;
    l = length(v.randFrames);
    for k = 1:l
        v.frames(:,:,k) = readframe(v.randFrames(k));
    end
    
    n = find(ndx(:,2)==i,1);
    for k = 1:length(validArenas)
        datName = ['Data Files/' filelist{n}];
        [ip,s] = getImageProcessingData(ip,datName,v.randFrames,n);
        
        [xrange,yrange] = getLimits(validArenas,vertices,k);
        for kk = 1:l
            [us] = selectHeadBody(us,v,s,xrange,yrange,n,kk);
        end
        n = n+1;
    end
end

save('tmpGroundTruth.mat','ip','us')

end

function [ip,s] = getImageProcessingData(ip,datName,randomFrames,n)
load(datName,'s');
ip.xH(n,:) = s.Head.x(randomFrames);
ip.yH(n,:) = s.Head.y(randomFrames);
ip.x(n,:) = s.Center.xUS(randomFrames);
ip.y(n,:) = s.Center.yUS(randomFrames);
ip.ang(n,:) = myatan(ip.xH(n,:)-ip.x(n,:),ip.yH(n,:)-ip.y(n,:),'degrees',2);
end

function [xrange,yrange] = getLimits(validArenas,vertices,k)

m = validArenas(k);
xCent = mean(vertices{m}(:,1));
yCent = mean(vertices{m}(:,2));

if xCent<1024 && yCent<1024
    xrange = [1 1024];yrange = [1 1024];
elseif xCent<1024 && yCent>1024
    xrange = [1 1024];yrange = [1024 2048];
elseif xCent>1024 && yCent<1024
    xrange = [1024 2048];yrange = [1 1024];
else
    xrange = [1024 2048];yrange = [1024 2048];
end

end

function [us] = selectHeadBody(us,v,s,xrange,yrange,n,kk)

currFrame = v.randFrames(kk);

figure(1);set(gcf,'Position',[2 42 958 954])
imagesc(v.frames(:,:,kk));hold on;
%plot(s.Center.x(currFrame-91:15:currFrame-1),s.Center.y(currFrame-91:15:currFrame-1),'*k')
priorFrames = currFrame-91:15:currFrame-1;
for i = 1:7
    plot(s.Center.x(priorFrames(i)),s.Center.y(priorFrames(i)),'*','Color',1-i/7*ones(1,3))
end
axis([xrange yrange])
title('Zoom into fly; Press any key to continue')
pause

title('Select Head position. Press any key to continue')
[xHTmp,yHTmp] = getpts;
us.xH(n,kk) = xHTmp(end);
us.yH(n,kk) = yHTmp(end);
title('Select Centroid position. Press any key to continue')
[xTmp,yTmp] = getpts;
us.x(n,kk) = xTmp(end);
us.y(n,kk) = yTmp(end);
us.ang(n,kk) = myatan(us.xH(n,kk)-us.x(n,kk),us.yH(n,kk)-us.y(n,kk),'degrees',2);

close all
end
