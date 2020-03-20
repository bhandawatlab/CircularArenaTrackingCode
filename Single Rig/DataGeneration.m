function [] = DataGeneration()
addpath(genpath([fileparts(which(mfilename)) '\matlab']));
addpath(genpath([fileparts(which(mfilename)) '\util']));

s=dir('Videos/*video.ufmf');
filelist={s.name};
for i = 1:length(filelist)
    [readframe,~,~,~] = get_readframe_fcn(['Videos/' filelist{i}]);
end
close all
clear

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
regions = 1;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
InAll = cell(1,length(filelist));vertAll = cell(1,length(filelist));
for i = 1:length(filelist)
    [InAll{i},vertAll{i}] = getArenaBounds([],[],['Videos/' filelist{i}],regions);
end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sRing=dir('Videos/*Inside.mat');
% filelistRing={sRing.name};
% InAll = cell(1,length(filelist));vertAll = cell(1,length(filelist));
% for i = 1:length(filelist)
%     load(['Videos/' filelistRing{i}]);
%     InAll{i} = Inside;vertAll{i} = vertices;
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


delete(gcp('nocreate'))
parpool('local',6)

exp = cell(1,length(filelist));
ppm = ParforProgMon('Processing: ', length(filelist) , 1);

parfor i = 1:length(filelist)
    filename = filelist{i};
    [readframe,frames,~,headerinfo] = get_readframe_fcn(['Videos/' filename]);
    
    % create background
    [~,bgvid,bgvid2] = generateBackground(headerinfo,frames,readframe);
    
    Inside = InAll{i};
    vertices = vertAll{i};
    badFlies = [];
    % do tracking with background subtraction
    [fly] = Tracking2(Inside,vertices,readframe,bgvid,bgvid2,frames,i,regions); % binarize with 0.25
    
    cF = zeros(1,regions);cent = zeros(regions,2);rad = zeros(1,regions);
    for k = 1:regions
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
% load('temp.mat');

% prompt users to eliminate bad flies
figure;
badFlies = cell(1,length(filelist));
for i = 1:length(filelist)
    for k = 1:regions
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
    badFlies{i} = input(['What Flies to Remove? (list 1-' num2str(regions) ')']);
end

% conduct post hoc corrections
close all
flyProc = cell(1,regions);
for i = 1:length(filelist)
    vertices = vertAll{i};
    [flyProc{i}] = PostHocAnalysis(exp{i}.fly,exp{i}.sArena,vertices,badFlies{i},regions);
    
    for k = 1:length(flyProc{i})
        if any(k == badFlies{i})
        else
            if ~isfield(exp{i}.fly{k},'Center')
                badFlies{i} = [badFlies{i}, k];
            else
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

                data(1,:) = s.Center.x;
                data(2,:) = s.Center.y;
                data(3,:) = s.AngVec;
                data(4,1:end-1) = s.Kinematics.thrust;
                data(5,1:end-1) = s.Kinematics.slip;
                data(6,1:end-1) = s.Kinematics.yaw;

                sArena.cF = exp{i}.sArena.cF;
                sArena.rad = exp{i}.sArena.rad(k);
                sArena.arenaCent = exp{i}.sArena.arenaCent(:,k);

                orginalVid = filelist{i};

                foldername = ['Analysis/' filelist{i}(1:end-11) '_' num2str(k)];
                plottingFunctions(s,sArena,length(s.Center.x),30,foldername)
                save(['Data Files/' filelist{i}(1:end-11) '_' num2str(k)],'s','sArena','data','orginalVid')
            end
        end
    end
    close all
end
close all

tic
parfor i = 1:length(filelist)
    filenamebgs = ['Videos/' filelist{i}];
    Inside = InAll{i};
    replayVideo(flyProc{1,i},Inside,filenamebgs,badFlies{i},regions)
end
toc

end