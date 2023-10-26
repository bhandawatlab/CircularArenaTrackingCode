function [] = plotXYTracks_fromRaw(dateVal, saveState)
% dateVal = '230421';
% saveState = 1;
currentDir = pwd;

cd(['D:\4-rig data\' dateVal 'Data Not Done\Data Files'])
s=dir();
filelist={s.name};
filelist = filelist(3:end);
s2 = cell(numel(filelist), 2);
for ww = 1:numel(filelist)
    load(filelist{ww})
    s2{ww, 1} = filelist{ww};
    s2{ww, 2} = s;
%     s2{ww, 3} = data;
end % ww

xCent = sArena.arenaCent(1, 1); 
yCent = sArena.arenaCent(2, 1); 


subplotRow = 4;
subplotNum = round(numel(filelist)/subplotRow);
if subplotNum*subplotRow < numel(filelist)
    subplotRow = 5;
else
    % do nothing
end
border = 0.985;%1./1.02;
figure;set(gcf,'Position',[2 42 838 924]); hold on;
for ww = 1:numel(filelist)
    lightOn = s2{ww, 2}.LightOn;
%     x = s2{ww, 2}.Center.x; y = s2{ww, 2}.Center.y;
    x = s2{ww, 2}.Center.x-xCent; y = s2{ww, 2}.Center.y;    
%     x = (s2{ww, 3}(1, :));
%     y = (s2{ww, 3}(2, :));
    
    
%% Normalize:
% %     zi = (xi – min(x)) / (max(x) – min(x))
      xy = zeros(2, numel(x));
      for yy = 1:numel(x)
          xy(1, yy) = (x(yy)-min(x))/(max(x)-min(x))*7.8;
          xy(2, yy) = (y(yy)-min(y))/(max(y)-min(y))*7.8;
      end % yy
      x = xy(1, :)-4; y = xy(2, :)-4;

    
    subplot(subplotNum, subplotRow, ww); hold on;
    tempName = s2{ww, 1}; underscoreIdx = find(tempName == '_');
    tempName(underscoreIdx) = ' ';
%     subtitle([tempName])
    plot(x(1:lightOn),y(1:lightOn),'g');hold on; % pre-stim
    plot(x(lightOn+1:end),y(lightOn+1:end),'r'); % stim
    plotCircle([0 0],4,100,'k');
    plotCircle([0 0],border,100,'b');
    xlim([-4 4]);ylim([-4 4])
    axis square
    title([tempName])
%     k = k+1;   
end

if saveState == 1
    fName = [dateVal ' Fly Tracks.pdf'];
    cd(['D:\4-rig data\' dateVal 'Data Not Done\Figures'])
    print('-painters','-dpdf',[fName])
else
    cd(currentDir)
end


end
