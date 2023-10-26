%Plotting functions
function [] = PlottingFunctionsGUI(Data,Arena,fHandle)
fs = 30;
spd = sqrt(diff(Data.Center.x).^2+diff(Data.Center.y).^2).*fs*Arena.cF;
%spd = sqrt((Data.Center.x-Arena.arenaCent(1)).^2+(Data.Center.y-Arena.arenaCent(2)).^2);

cla(fHandle.axThrust)
plot(fHandle.axThrust,Data.Kinematics.thrust);
hold(fHandle.axThrust,'on');
plot(fHandle.axThrust,spd);
title(fHandle.axThrust,'Thrust');
hold(fHandle.axThrust,'off');

cla(fHandle.axYaw)
plot(fHandle.axYaw,Data.Kinematics.yaw);
title(fHandle.axYaw,'Yaw');

cla(fHandle.axXY)
plot3(fHandle.axXY,1:1:length(Data.Center.x),Data.Center.x,Data.Center.y);
view(fHandle.axXY,[90,0])
title(fHandle.axXY,'Centroid Positions');
 
end