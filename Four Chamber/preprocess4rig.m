addpath(genpath(pwd))
foldName = 'D:\4 Chamber Data\210701\Data Files\';
behavior_list = dir([foldName '*.mat']);
tracks = {behavior_list.name};
numflies = length(tracks);
track_names = tracks';
m1=0;
fig_count=1;

currentFolder = pwd;
k = strfind(currentFolder, '\');
start = max(k) + 1;
r = 1.25/4;
attNdx = nan(numflies,2);

for ii=1:numflies % do this for each fly
    % get the tracks
    figure(fig_count)
    set(gcf,'Position',[434 58 1200 600],'color','white')
    set(gca,'Box','off','Xtick',[],'Ytick',[],'XColor',[1,1,1],'YColor',[1,1,1]);
    axis tight;
    
    
    filename_behavior = tracks{ii};
    load([foldName filename_behavior]);
    LightOn=s.LightOn;
    
    in = s.Distances.DistanceR<r;
    during = true(1,numel(in));
    during(1:LightOn) = false;
    
    firstEntry = find(in & during,1,'first');
    
    if ~isempty(firstEntry)
        attNdx(ii,1) = sum(in(1:firstEntry-1))./(firstEntry-1);                 % before
        attNdx(ii,2) = sum(in(firstEntry:end))./(numel(in)-firstEntry+1);       % after
    end
    
    j=rem(ii,6);
    if j==0
        j=6;
        fig_count=fig_count+1;
    end
    
    subplot(2,3,j)
    plot(s.Distances.xDist(1:LightOn),s.Distances.yDist(1:LightOn),'b');hold on;
    plot(s.Distances.xDist(LightOn:end),s.Distances.yDist(LightOn:end),'r');
    axis off;
    
    plotCircle([0,0],r,100,'k');
    plotCircle([0,0],1,100,'k');
    
    xlim([-1 1]); ylim([-1 1]);
    
    text(0,0,num2str(ii))
    title(sprintf ('%c', filename_behavior(1:end)), 'Interpreter', 'none')
    axis off;
    
    set(gcf,'NextPlot','add');
    axes;
    % h = title(sprintf ('%c', currentFolder(start:end)));
    h = title('LexAop-Ch;UAS-GtACR1;Orco-LexA (Control)');
    set(gca,'Visible','off');
    set(h,'Visible','on');
    set(h,'Position',[0.5,0.5]);
    clear rim
end


attNdx(all(isnan(attNdx),2),:) = [];
attNdxLin = reshape(attNdx,[],1);
id = {'Before','During'};

if size(attNdx,2)>=10
    % figure;suptitle(currentFolder(start:end))
    figure;suptitle('LexAop-Ch;UAS-GtACR1;Orco-LexA (Control)')
    [ss,ssPMC,moes] = dabest2(attNdx,id,'Y');
    ylabel('delta Attraction Ndx')
end


