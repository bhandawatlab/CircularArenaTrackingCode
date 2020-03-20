function [fly2] = PostHocAnalysis(fly,sArena,vertices,badFlies,regions)
addpath('C:\Nicholas Data\Image Processing New Rig\Util')
f = 30;
cF = sArena.cF;

for k = 1:regions
    if any(k == badFlies) || ~isfield(fly{k},'Center')
        fly{k} = [];
    else
        ang = fly{k}.ang;
        xTmp = fly{k}.Center.x;
        x = [xTmp xTmp(end).*ones(1,length(ang)-length(xTmp))];
        yTmp = fly{k}.Center.y;
        y = [yTmp yTmp(end).*ones(1,length(ang)-length(yTmp))];
        mjrALTmp = fly{k}.MjrAxs./2;
        mjrAL = [mjrALTmp mean(mjrALTmp).*ones(1,length(ang)-length(mjrALTmp))];
        minALTmp = fly{k}.MinAxs./2;
        minAL = [minALTmp mean(minALTmp).*ones(1,length(ang)-length(minALTmp))];
        
        x = interp1q(find(x~=0)',x(x~=0)',[1:1:length(x)]')';
        y = interp1q(find(y~=0)',y(y~=0)',[1:1:length(y)]')';
        fly{k}.Center.x = x;
        fly{k}.Center.y = y;
        
        mjrAL = interp1q(find(mjrAL~=0)',mjrAL(mjrAL~=0)',[1:1:length(mjrAL)]')';
        minAL = interp1q(find(minAL~=0)',minAL(minAL~=0)',[1:1:length(minAL)]')';
        fly{k}.MjrAxs = mjrAL.*2;
        fly{k}.MinAxs = minAL.*2;
        
        len = length(x);

        %calculating displacement, thrust, slip
        dspmt = sqrt(diff(x).^2 + diff(y).^2);
        mvmtAngle = ang(1:len-1) - myatan(diff(x),diff(y),'degrees',2);
        speed = dspmt.*f*cF;
        thrustV = dspmt.*cosd(mvmtAngle)*f*cF;

        tmp = find(thrustV<-5);
        endNdx = [find(diff(tmp)>3) length(tmp)];
        startNdx = [1 endNdx(1:end-1)+1];

        angVecCor = ang;
        for i = 1:length(endNdx)
            angVecCor(tmp(startNdx(i)):tmp(endNdx(i))) = ang(tmp(startNdx(i)):tmp(endNdx(i)))+180;
        end
        angVecCor(angVecCor>360) = angVecCor(angVecCor>360)-360;

        for i = 1:30
            %correct for locations where large yall is found twice very close
            % correct by thrust and change in yaw
            [angVecCor] = correctForOpposites(angVecCor,30);
            [angVecCor] = correctForSame(angVecCor,30);
        end
        for i = 1:11
            [angVecCor] = correctForOpposites(angVecCor,3000);
            [angVecCor] = correctForSame(angVecCor,3000);
        end

        yaw = diff(angVecCor);
        yaw(yaw>180) = yaw(yaw>180)-360;
        yaw(yaw<-180) = yaw(yaw<-180)+360;
        ndx = find(abs(yaw)>=90);

        %calculating displacement, thrust, slip
        dspmt = sqrt(diff(x).^2 + diff(y).^2);
        mvmtAngle = angVecCor(1:len-1) - myatan(diff(x),diff(y),'degrees',2);
        speed = dspmt.*f*cF;
        thrustV = dspmt.*cosd(mvmtAngle)*f*cF;
        
        
        for i = 1:length(ndx)
            
            angVecCor(1:ndx(i)) = angVecCor(1:ndx(i))+180;
            angVecCor(angVecCor>360) = angVecCor(angVecCor>360)-360;
            [angVecCor] = VelCheck(speed,thrustV,angVecCor,[0 len-1],1);
            
            %calculating displacement, thrust, slip
            dspmt = sqrt(diff(x).^2 + diff(y).^2);
            mvmtAngle = angVecCor(1:len-1) - myatan(diff(x),diff(y),'degrees',2);
            speed = dspmt.*f*cF;
            thrustV = dspmt.*cosd(mvmtAngle)*f*cF;
            
        end
        
        [angVecCor] = VelCheck(speed,thrustV,angVecCor,[0 len-1],1);
        
        %calculating displacement, thrust, slip
        mvmtAngleNew = angVecCor(1:len-1) - myatan(diff(x),diff(y),'degrees',2);
        thrustVNew = dspmt.*cosd(mvmtAngleNew)*f*cF;
        slipVNew = dspmt.*sind(mvmtAngleNew)*f*cF;
        yaw = diff(angVecCor);
        yaw(yaw>180) = yaw(yaw>180)-360;
        yaw(yaw<-180) = yaw(yaw<-180)+360;

        figure;subplot(2,1,1);plot(thrustVNew);
        title('Processed Thrust')
        subplot(2,1,2);plot(thrustV);
        title('PreProcessed Thrust')
        set(gcf,'Position',[1361 572 560 420])
        
        figure;plot(yaw);
        title('Processed yaw')
        set(gcf,'Position',[1361 46 560 420])
        
        %calculate major and minor axis endpoints
        xMajor1 = fly{k}.Center.x + mjrAL .* cosd(angVecCor);
        yMajor1 = fly{k}.Center.y + mjrAL .* sind(angVecCor);

        %calculate head
        fly{k}.Head.x = xMajor1;
        fly{k}.Head.y = yMajor1;
        fly{k}.ang = angVecCor;
        fly{k}.thrust = thrustVNew;
        fly{k}.slip = slipVNew;
        fly{k}.yaw = yaw;
        

        delX = (fly{k}.Center.x-sArena.arenaCent(1,k))*cF;
        delY = (fly{k}.Center.y-sArena.arenaCent(2,k))*cF;
        R = sqrt(delX.^2+delY.^2);

        fly{k}.xDist = delX./40;
        fly{k}.yDist = delY./40;
        fly{k}.DistanceR = R./40;
        
        fly{k}.Center.xUS = fly{k}.Center.x;
        fly{k}.Center.yUS = fly{k}.Center.y;
        
        xTemp = fly{k}.Center.x;
        yTemp = fly{k}.Center.y;
        fineX = smooth(xTemp,6,'loess')';
        roughX = smooth(xTemp,30,'loess')';
        fineY = smooth(yTemp,6,'loess')';
        roughY = smooth(yTemp,30,'loess')';
        temp3 = abs(roughX-xTemp);
        roughX(temp3>1.5)=fineX(temp3>1.5);
        fly{k}.Center.x = roughX;
        temp3 = abs(roughY-yTemp);
        roughY(temp3>1.5)=fineY(temp3>1.5);
        fly{k}.Center.y = roughY;

        fly{k}.Loss.x = sum((xTemp-fly{k}.Center.x).^2);
        fly{k}.Loss.y = sum((yTemp-fly{k}.Center.y).^2);
        
        fly{1,k}.LightOn = 5400;
        if fly{k}.LightOn>6000
            fly{k}.Status = 'bad';
        else
            fly{k}.Status = 'good';
        end
        
    end
end

fly2 = fly;
end


function [angVecCor] = correctForOpposites(angVecCor,durLim)

yaw = diff(angVecCor);
yaw(yaw>180) = yaw(yaw>180)-360;
yaw(yaw<-180) = yaw(yaw<-180)+360;
ndx = find(yaw>150);
ndx2 = find(yaw<-150);
if ~isempty(ndx) && ~isempty(ndx2)
    [PtBef,PtAft,~] = findBeforeAfter(ndx,ndx2,'both');
    [dur,ndx3] = min(abs([ndx-PtBef;ndx-PtAft]));
    
    comb2Cons = [PtBef(ndx3==1) ndx(ndx3==2);
        ndx(ndx3==1) PtAft(ndx3==2)];
    comb2Cons = sort(comb2Cons,2);
    comb2Cons(:,dur>durLim) = [];
    
    comb2Cons = sort(comb2Cons,2);
    
    for i = 1:size(comb2Cons,2)
        startNdx = comb2Cons(1,i);
        endNdx = comb2Cons(2,i);
        if i<size(comb2Cons,2)
            if all(endNdx<comb2Cons(1,i+1))
                angVecCor(startNdx+1:endNdx) = angVecCor(startNdx+1:endNdx)+180;
            end
        else
            angVecCor(startNdx+1:endNdx) = angVecCor(startNdx+1:endNdx)+180;
        end
    end
end
angVecCor(angVecCor>360) = angVecCor(angVecCor>360)-360;

end

function [angVecCor] = correctForSame(angVecCor,durLim)

yaw = diff(angVecCor);
yaw(yaw>180) = yaw(yaw>180)-360;
yaw(yaw<-180) = yaw(yaw<-180)+360;
ndx = find(abs(yaw)>100);
if ~isempty(ndx) && ~isempty(ndx)
    [PtBef,PtAft,~] = findBeforeAfter(ndx,ndx,'both');
    [dur,ndx3] = min(abs([ndx-PtBef;ndx-PtAft]));
    
    comb2Cons = [PtBef(ndx3==1) ndx(ndx3==2);
        ndx(ndx3==1) PtAft(ndx3==2)];
    comb2Cons = sort(comb2Cons,2);
    comb2Cons(:,dur>durLim) = [];
    
    comb2Cons = sort(comb2Cons,2);
    
    for i = 1:size(comb2Cons,2)
        startNdx = comb2Cons(1,i);
        endNdx = comb2Cons(2,i);
        if i<size(comb2Cons,2)
            if all(endNdx<comb2Cons(1,i+1))
                angVecCor(startNdx+1:endNdx) = angVecCor(startNdx+1:endNdx)+180;
            end
        elseif isnan(startNdx)
            angVecCor(endNdx) = angVecCor(endNdx)+180;
        elseif isnan(endNdx)
            angVecCor(startNdx+1) = angVecCor(startNdx+1)+180;
        else
            angVecCor(startNdx+1:endNdx) = angVecCor(startNdx+1:endNdx)+180;
        end
    end
end
angVecCor(angVecCor>360) = angVecCor(angVecCor>360)-360;

end


