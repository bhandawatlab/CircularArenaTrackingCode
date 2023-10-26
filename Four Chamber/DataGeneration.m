%function [] = DataGeneration()
%% This code works correctly as of 230511, please do not change anything (SPW)

clear all; close all;
%--------------------------------------------------------------------------
addpath(genpath([fileparts(which(mfilename)) '\matlab']));
addpath(genpath([fileparts(which(mfilename)) '\util']));
addpath('C:\Users\lt532\Desktop\New Code (4-rig)_ data acquisition\ufmfMatFiles\filehandling')
addpath('C:\Users\lt532\Desktop\New Code (4-rig)_ data acquisition\ufmfMatFiles\misc')
addpath('C:\Users\lt532\Desktop\4Rig_ImageProcessing\Util')
addpath('D:\Image Processing')
addpath('C:\Users\lt532\Desktop\ORN-Optogenetics-main\Util\ParforProgMon-master')
addpath('C:\Users\sbm94\Desktop\DrosoRT-master')
addpath('C:\Users\sbm94\Desktop\ORN-Optogenetics-mainFinalMarch2022\PlottingFunctions')
addpath('C:\Users\spw43\Desktop\ORN-Optogenetics-main\Util')
%%
tic
%%

dateVal = "20-Oct-2023";
[fileslist2, badFlies] = getBadFlies(dateVal);

%%

s=dir('Videos/*video.ufmf');
filelist={s.name};


%%
fileNames = cell(numel(filelist), 3);
for ww = 1:numel(filelist)
    tempName = filelist{ww}; underscoreIdx = find(tempName == '_');
   fileNames{ww, 1} = filelist{ww}(1:end-5);
   fileNames{ww, 2} = str2num(string(tempName(underscoreIdx(2)+1:underscoreIdx(3)-1)));
   fileNames{ww, 3} = ww;    
end % ww
fileNamesIdx = cell2mat(fileNames(:, 2)); fileNamesIdx = sort(fileNamesIdx);


%%
errorIdx = [];
for i = 1:length(filelist)
    try
    [readframe,~,~,~] = get_readframe_fcn(['Videos/' filelist{i}]);
    catch
        errorIdx = [];
    end
end
% close all
% clear

if ~exist([pwd '\Analysis'], 'dir')
    mkdir([pwd '\Analysis'])
end
if ~exist([pwd '\Figures'], 'dir')
    mkdir([pwd '\Figures'])
end
if ~exist([pwd '\Data Files'], 'dir')
    mkdir([pwd '\Data Files'])
end

s=dir('Videos/*video.ufmf');
filelist={s.name};

% Add a check chere for any video_inside.mat files and if so, skip this
% section and pre-load everything instead. 

sRing=dir('Videos/*Inside.mat');
filelistRing={sRing.name};
if isempty(filelistRing)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
InAll = cell(1,length(filelist)); vertAll = cell(1,length(filelist));
i = 1;
[InAll{i},vertAll{i}] = getArenaBounds([],[],['Videos/' filelist{i}]);

for i = 1:length(filelist)
    try
        InAll{i} = InAll{1};
        vertAll{i} = vertAll{1};
        [InAll{i},vertAll{i}] = getArenaBounds(InAll{1},vertAll{1},['Videos/' filelist{i}]);
    catch
    end
end

else
 % Do nothing    
 sRing=dir('Videos/*Inside.mat');
filelistRing={sRing.name};
InAll2 = cell(1,length(filelist));
vertAll = cell(1,length(filelist));
for i = 1:length(filelist)
    load(['Videos/' filelistRing{i}]);
    InAll2{i} = Inside;vertAll{i} = vertices;
end
InAll = InAll2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sRing=dir('Videos/*Inside.mat');
filelistRing={sRing.name};
InAll2 = cell(1,length(filelist));
vertAll = cell(1,length(filelist));
for i = 1:length(filelist)
    load(['Videos/' filelistRing{i}]);
    InAll2{i} = Inside;vertAll{i} = vertices;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


delete(gcp('nocreate'))
%parpool('local',6) old, might be neccessary
parpool('local',4)

exp = cell(1,length(filelist));
ppm = ParforProgMon('Processing: ', length(filelist) , 1);

aa = 1; bb = length(filelist);
parfor i = aa:bb
    filename = filelist{i};
    [readframe,frames,~,headerinfo] = get_readframe_fcn(['Videos/' filename]);

    % create background
    [~,bgvid,bgvid2] = generateBackground(headerinfo,frames,readframe);

    Inside = InAll{i};
    vertices = vertAll{i};
    badFlies = [];
    % do tracking with background subtraction
    [fly] = Tracking2(Inside,vertices,readframe,bgvid,bgvid2,frames,i); % binarize with 0.25

    cF = zeros(1,4);cent = zeros(4,2);rad = zeros(1,4);
    for k = 1:4
        I = Inside(:,:,k);
        stats = regionprops(I, 'EquivDiameter', 'Centroid');
        cF(k) = 80./stats.EquivDiameter;
        rad(k) = stats.EquivDiameter./2;
        cent(k,:) = stats.Centroid;
    end
    cF = mean(cF);

    exp{i}.fly = fly;
    exp{i}.Inside = Inside;
    exp{i}.vertices = vertices;
    exp{i}.filename = filename;
    exp{i}.filenamebgs = filename;
    exp{i}.badFlies = badFlies;
    exp{i}.sArena.cF = cF;
    exp{i}.sArena.rad = rad;
    exp{i}.sArena.arenaCent = cent';
    ppm.increment();

end

save('temp.mat','-v7.3');

% --------------------------------------------------------------------------
load('temp.mat');
% If the day value is a single digit (1-9), then it messes up here. Fix
% eventually. low priority.

if isempty(badFlies)
   newDateVal = tempName(1:6);
   oldDateVal = char(dateVal);
   hyphenIdx = find(oldDateVal == '-');
   if newDateVal(end-1) == '0'
   dateVal = string([newDateVal(end) oldDateVal(hyphenIdx(1):end)]);
   else
       dateVal = string([newDateVal(end-1:end) oldDateVal(hyphenIdx(1):end)]);
   end
   [fileslist2, badFlies] = getBadFlies(dateVal);
    save('temp.mat','-v7.3');
   
else
    % do nothing
end



tempDateVal = char(dateVal);
hyphenVal = find(tempDateVal == '-'); tempDateVal = string(tempDateVal(1:hyphenVal(1)-1));
% tempDateVal = string(tempDateVal(1:2)); 
if numel(char(tempDateVal)) == 1
    tempDateVal = string(['0' char(tempDateVal)]);
    singleDigitTrigger = 1;
else
    singleDigitTrigger = 0;
end
currentDir = pwd; tempIdx = find(currentDir == '\'); dateValCheck = string([currentDir(tempIdx(2)+5:tempIdx(2)+6)]);
if dateValCheck == tempDateVal
    if singleDigitTrigger == 1
        tempDateVal = char(dateVal);
        dateVal = string(['0' tempDateVal]);
    else
    end
    % do nothing
else
    tempDateVal = char(dateVal);
    dateVal = string([char(dateValCheck) tempDateVal(3:end)]);
end


figure;
for i = aa:bb
    for k = 1:4
        if isfield(exp{i}.fly{k},'Center')
            half = ceil(length(exp{i}.fly{k}.Center.x)./2);
            
            x = exp{i}.sArena.arenaCent(1,k);
            y = exp{i}.sArena.arenaCent(2,k);
            r = ((max(vertAll{i}{k}(:,1))-min(vertAll{i}{k}(:,1)))./2+(max(vertAll{i}{k}(:,2))-min(vertAll{i}{k}(:,2)))./2)./2;
            [xLight,yLight] = circle(x,y,r.*1.25/4);
            subplot(2,2,k);hold off
            plot(vertAll{i}{k}(:,1),vertAll{i}{k}(:,2),'k');
            hold on
            plot(exp{i}.fly{k}.Center.x(1:half),exp{i}.fly{k}.Center.y(1:half),'g');
            plot(exp{i}.fly{k}.Center.x(half:end),exp{i}.fly{k}.Center.y(half:end),'r');
            plot(xLight,yLight,'k');hold off
            title(['Fly ' num2str(k)])
            axis tight
        else
            subplot(2,2,k);hold off
            plot(0,0);
        end
    end
end % break point here to debug outliers**********

% dateVal = "15-Jun-2023";
[fileslist2, badFlies] = getBadFlies(dateVal);
itxIdx = intersect(string(fileNames(:, 1)), string(fileslist2)');

badFlies_old = badFlies;
for i = 1:numel(filelist)
    ndx = find(strcmp(fileslist2,filelist{i}(1:end-5)));
    badFlies(i) = badFlies_old(ndx);
end


% conduct post hoc corrections
flyProc = cell(1,4);
for i = aa:bb
    vertices = vertAll{i};
    [flyProc{i}] = PostHocAnalysis(exp{i}.fly,exp{i}.sArena,vertices,badFlies{i});
    
    for k = 1:length(flyProc{i})
        if ~isempty(flyProc{i}{k})
            disp(["K value:" num2str(k)])
            s = [];
            s.Kinematics.thrust = flyProc{1,i}{1,k}.thrust;
            s.Kinematics.slip = flyProc{1,i}{1,k}.slip;
            s.Kinematics.yaw = flyProc{1,i}{1,k}.yaw;
            s.Distances.DistanceR = flyProc{1,i}{1,k}.DistanceR;
            s.Distances.xDist = flyProc{1,i}{1,k}.xDist;
            s.Distances.yDist = flyProc{1,i}{1,k}.yDist;
            s.MjrAxs = flyProc{1,i}{1,k}.MjrAxs;
            s.MinAxs = flyProc{1,i}{1,k}.MinAxs;
            s.AngVec = flyProc{1,i}{1,k}.ang;
            s.Center = flyProc{1,i}{1,k}.Center;
            s.Head = flyProc{1,i}{1,k}.Head;
            s.Status = flyProc{1,i}{1,k}.Status;
            s.LightOn = 5400;%exp{i}.fly{1,k}.LightOn;
            s.Loss = flyProc{1,i}{1,k}.Loss;
            
            if numel(s.Center.x) == 10799
                data(1,:) = s.Center.x; %%
                data(2,:) = s.Center.y;
                data(3,:) = s.AngVec;
                data(4,1:end-1) = s.Kinematics.thrust;
                data(5,1:end-1) = s.Kinematics.slip;
                data(6,1:end-1) = s.Kinematics.yaw;
            else
                
            end
            sArena.cF = exp{i}.sArena.cF;
            sArena.rad = exp{i}.sArena.rad(k);
            sArena.arenaCent = exp{i}.sArena.arenaCent(:,k);
            
            orginalVid = filelist{i};
            
            foldername = ['Analysis/' filelist{i}(1:end-11) '_' num2str(k)];
            save(['Data Files/' filelist{i}(1:end-11) '_' num2str(k)],'s','sArena','data','orginalVid')
            try
            plottingFunctions(s,sArena,length(s.Center.x),30,foldername)
            catch
            end
            %if size(s.Center.x) == 10799
           
            %else
            % do nothing
            %end
            
        end
    end
    close all
end
close all

%% Generate the XY tracks from raw data:
currentDir = pwd; tempIdx = find(currentDir == '\'); 
dateVal2 = [currentDir(tempIdx(2)+1:tempIdx(2)+6)];
plotXYTracks_fromRaw(dateVal2, 1);
cd(currentDir)

%%
toc
