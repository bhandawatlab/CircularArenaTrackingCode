clc;
clear all;
close all;
%%
dateVal = "08-May-2023";
dateVal2 = '230508';
cd(['D:\4-rig data\' dateVal2 'Data Not Done'])
[filelist, badFlies] = getBadFlies(dateVal);
s=dir('Videos/*video.ufmf');
filelist={s.name};

load('temp.mat');
% prompt users to eliminate bad flies
figure;
%[fileslist2, badFlies] = getBadFlies(dateVal)
%badFlies = cell(1,length(filelist));
%aa = 1; bb = length(filelist);
%for i = 1:length(filelist)
%dateVal
%[fileslist2, badFlies] = getBadFlies(dateVal)
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
    filelist(i)
    %badFlies{i} = input('What Flies to Remove?');
end % break point here to debug outliers********** 


% conduct post hoc corrections
close all
flyProc = cell(1,4);
% for i = 1:length(filelist)
% disp("Here...")
for i = aa:bb
    vertices = vertAll{i};
    [flyProc{i}] = PostHocAnalysis(exp{i}.fly,exp{i}.sArena,vertices,badFlies{i});
    
    %disp("Working Till Here")
    for k = 1:length(flyProc{i})
        %if any(k == badFlies{i})
        %else
            %if ~isfield(exp{i}.fly{k},'Center')
            %    badFlies{i} = [badFlies{i}, k];
            %else
            %disp("here")
            if isempty(flyProc{i}{k})
                %do nothing
            else
                disp("K value:")
                disp(k)
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
                
                %if size(s.Center.x) > 10000
                data(1,:) = s.Center.x; %%
                data(2,:) = s.Center.y;
                data(3,:) = s.AngVec;
                data(4,1:end-1) = s.Kinematics.thrust;
                data(5,1:end-1) = s.Kinematics.slip;
                data(6,1:end-1) = s.Kinematics.yaw;
                %else
                    
                %end
                sArena.cF = exp{i}.sArena.cF;
                sArena.rad = exp{i}.sArena.rad(k);
                sArena.arenaCent = exp{i}.sArena.arenaCent(:,k);

                orginalVid = filelist{i};

                foldername = ['Analysis/' filelist{i}(1:8) '_' num2str(i) '_' num2str(k)];
                plottingFunctions(s,sArena,length(s.Center.x),30,foldername)
                %if size(s.Center.x) == 10799
                
                save(['Data Files/' foldername(10:end)],'s','sArena','data','orginalVid')
%                 save(['Data Files/' filelist{i}(1:10) '_' num2str(k)],'s','sArena','data','orginalVid')
                disp(['saved: ' ['Data Files/' foldername(10:end) '_' num2str(k)]])
                clear data
                %else
                   % do nothing 
                %end
        end
    end
    close all
end
close all

plotXYTracks_fromRaw(dateVal2, 1);


