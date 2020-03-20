function [fly] = Tracking2(Inside,vertices,readframe,bgvid,bgvid2,frames,fileNum,regions)

fly = cell(1,regions);maxVals = cell(1,regions);minVals = cell(1,regions);
for k = 1:regions
    fly{k}.ang = zeros(1,frames);
    maxVals{k} = ceil(max(vertices{k}));
    minVals{k} = floor(min(vertices{k}));
end
se = strel('disk',6);
%tic
for i = 1:frames %loop from start-->endframe
    IBase = readframe(i);
    if i<floor(frames./2)
        IBase = imsubtract(IBase,bgvid);
    else
        IBase = imsubtract(IBase,bgvid2);
    end
    IBase = IBase./max(IBase(:));
    %tic
    for k = 1:regions
        I = IBase.*Inside(:,:,k);
        I2 = I(minVals{k}(2):maxVals{k}(2),minVals{k}(1):maxVals{k}(1));
        idx = [];thresh = 0.25;
        while isempty(idx) && thresh>0
            if thresh>0
                J = imbinarize(I2,thresh);
                temp = imdilate(J,se);
                K= bwareaopen(temp,8);
                L = imfill(K,'holes');
                %compute largest connected group of pixels to keep as fly
                CC = bwconncomp(L);
                numPixels = cellfun(@numel,CC.PixelIdxList);
                [~,idx] = max(numPixels);
                thresh = thresh-0.05;
            end
        end
        if ~isempty(idx)
            img = zeros(size(L));
            img(CC.PixelIdxList{idx}) = 1;
            
            %calculate centroid, major, minor, and orientation of fly
            stats = regionprops(img,'Centroid',...
                'MajorAxisLength','MinorAxisLength','Orientation');
            stats.MajorAxisLength = stats.MajorAxisLength-6;
            stats.MinorAxisLength = stats.MinorAxisLength-6;
            
            %centroid (take mean)
            cent = cat(1, stats.Centroid);
            index = find(isnan(cent(:,1)));
            cent(isnan(cent(:,1)),:) = [];
            centroid = mean(cent,1);
            fly{k}.Center.x(i) = centroid(1)+minVals{k}(1);
            fly{k}.Center.y(i) = centroid(2)+minVals{k}(2);
            
            %major and minor axis (take mean)
            mjrAL = cat(1, stats.MajorAxisLength);
            mjrAL(index) = [];
            mjrAL = mean(mjrAL)/2;      %radius
            fly{k}.MjrAxs(1,i) = mjrAL*2;    %diameter
            minAL = cat(1, stats.MinorAxisLength);
            minAL(index) = [];
            minAL = mean(minAL)/2;      %radius
            fly{k}.MinAxs(1,i) = minAL*2;    %diameter
            
            %orientation (take orientation associated with largest major axis value)
            ort = cat(1, stats.Orientation);
            idx = cat(1, stats.MajorAxisLength) == max(cat(1,stats.MajorAxisLength));
            ang = ort(idx);
            
            %Format where orientation is the angle is between the negative
            %horizontal axis and the major axis.
            ang = 180-ang;
            if ang<0
                ang = 360+ang;
            end
            angVec = ang+180;
            angVec(angVec>360) = angVec(angVec>360)-360;
            
            
            % If the general oval shape is more circular, then the head
            % position may be in either the major or minor axis
            if minAL/mjrAL > 0.8 && (i>1) && length(fly{k}.ang)>=(i-1)
                oldAng = fly{k}.ang(i-1);
                
                angCases = [ang,ang+90,ang+180,ang+270];
                angCases(angCases>360) = angCases(angCases>360)-360;
                
                [~,ndx] = min(min(abs([oldAng-angCases;oldAng-angCases+360])));
                angVec = angCases(ndx);
            else
                % choose the head position based on closeness of major axis
                % end points to the front of the location of the brightest
                % pixel
                angVec = findHeadByIntensity(I2,centroid,angVec,mjrAL);
            end
            
            angVec = [angVec,angVec+180];
            
            %calculate major and minor axis endpoints
            xMajor1 = centroid(1)+minVals{k}(1) + mjrAL * cosd(angVec(1));
            yMajor1 = centroid(2)+minVals{k}(2) + mjrAL * sind(angVec(1));
            %calculate major and minor axis endpoint candidates
            xMajor2 = centroid(1)+minVals{k}(1) + mjrAL * cosd(angVec(2));
            yMajor2 = centroid(2)+minVals{k}(2) + mjrAL * sind(angVec(2));
            
            % hungarian algorithm
            if i > 1
                currPt = [xMajor1,yMajor1;xMajor2,yMajor2];
                lastpt = [fly{k}.Head.y(i-1) fly{k}.Head.y(i-1)];
                D1 = pdist2(currPt,lastpt,'euclidean');
                [newNdx,finCost] = munkres(D1);
                newNdx = boolean(newNdx);
                
                % assign head positions
                fly{k}.Head.x(i) = currPt(newNdx,1);
                fly{k}.Head.y(i) = currPt(newNdx,2);
                
                % assign angle
                fly{k}.ang(i) = angVec(newNdx);
            end
            
            %calculate head
            fly{k}.Head.x(i) = xMajor1;
            fly{k}.Head.y(i) = yMajor1;
            
            fly{k}.ang(i) = angVec;
            
        end
        
    end
    %end
    
    if mod(i,1000)==1
        display([num2str(fileNum) ': ' num2str(i)]);%toc
    end
    
end

end

function [ang] = findHeadByIntensity(I,centroid,ang,mjrAL)

%calculate possible head positions as ends of major axis
xMajor1 = centroid(1) + mjrAL * cosd(ang);
xMajor2 = centroid(1) - mjrAL * cosd(ang);
yMajor1 = centroid(2) + mjrAL * sind(ang);
yMajor2 = centroid(2) - mjrAL * sind(ang);

% calculate the maximum intensity pixel in the general
% viscinity of the fly
x = round(centroid(2)-50:centroid(2)+50);
y = round(centroid(1)-50:centroid(1)+50);
x(x<1) = 1;x(x>size(I,1)) = size(I,1);
y(y<1) = 1;y(y>size(I,2)) = size(I,2);
[r,c] = find(I==max(max(I(x,y))));

% calculate the distance of the each possible head position
% from the maximum intensity pixel
hPos1 = norm([xMajor1,yMajor1]-[c,r]);
hPos2 = norm([xMajor2,yMajor2]-[c,r]);

% the head is likely at the position closest to the maximum
% intensity pixel
if hPos2<hPos1
    ang = ang+180;
end
ang(ang>360) = ang(ang>360)-360;

end


