function [] = plottingFunctions(s,sArena,len,fs,foldername)

if ~(exist(foldername,'dir')==7)
    mkdir(foldername)
end

time = (1:1:len)./fs;
x = s.Center.xUS;
y = s.Center.yUS;
speed = sqrt(s.Kinematics.thrust.^2+s.Kinematics.slip.^2);
mvmtAngle = s.AngVec(1:len-1) - myatan(diff(x),diff(y),'degrees',2);
mvmtAngle(mvmtAngle>180) = mvmtAngle(mvmtAngle>180)-360;
mvmtAngle(mvmtAngle<-180) = mvmtAngle(mvmtAngle<-180)+360;

figure(6)
subplot(221)
histogram(s.Kinematics.thrust(1:end))
title('Thrust')
xlabel('mm/s')
subplot(222)
histogram(s.Kinematics.slip(1:end))
title('Slip')
xlabel('mm/s')
subplot(223)
plot(time(1:end-1),s.Kinematics.thrust(1:end))
title('Thrust')
xlabel('Time(s)')
subplot(224)
histogram(mvmtAngle(1:end))
title('Movement Angle')
xlabel('Angle')
tempName = [foldername(10:end)]; underscoreIdx = find(tempName == '_'); tempName(underscoreIdx) = ' ';
sgtitle(tempName)
% print([foldername '/General_Histograms'],'-dpdf')
% savefig([foldername '/General_Histograms.fig'], 'fig')
print('-painters','-bestfit', '-dpdf',[foldername '/General_Histograms - ' foldername(10:end)]);
%     savefig([foldername '/General_Histograms' '.fig'])

 
figure(7)
plot(time(1:end-1),s.Kinematics.thrust(1:end));hold on
plot(time(1:end-1),speed(1:end),'k')
% plot(time,s.Flags.ratio(1:len),'*-r')
% plot(time,s.Flags.size(1:len),'*-g')
% s.Flags.thrust(s.Kinematics.thrust<-10) = 1;
% s.Flags.thrust(s.Kinematics.thrust>20) = 1;
% [~, idx] = find(s.Flags.thrust == 1);
% plot(idx/fs-1/fs, s.Kinematics.thrust(idx),'*g')
%legend('Crit Points', 'Thrust', 'Speed', 'Outliers')
legend('Thrust', 'Speed')
xlabel('Time (s)')
ylabel('Speed (mm/s)')
sgtitle(tempName)
hold off
% print([foldername '/Thrust_time_Course'],'-dpdf')
print('-painters','-bestfit', '-dpdf',[foldername '/Thrust_time_Course - ' foldername(10:end)]);
 
figure(8)
%plot(time,s.Flags.ratio(1:len),'*-r')
plot(time(1:end-1),s.Kinematics.slip(1:end))
hold on
plot(time(1:end-1),speed(1:end),'k')
s.Flags.slip(s.Kinematics.slip<-10) = 1;
s.Flags.slip(s.Kinematics.slip>10) = 1;
%[~, idx] = find(s.Flags.slip == 1);
%plot(idx/fs-1/fs, s.Kinematics.slip(idx),'*g')
%legend('Crit Points', 'Slip', 'Speed', 'Outliers')
legend('Slip', 'Speed')
xlabel('Time (s)')
ylabel('Speed (mm/s)')
sgtitle(tempName)
hold off
% print([foldername '/Slip_time_Course'],'-dpdf')
print('-painters','-bestfit', '-dpdf',[foldername '/Slip_time_Course - ' foldername(10:end)]);
 
%Plotting general features
figure(9)
subplot(221)
plot(time, s.MjrAxs(1,1:len).*sArena.cF, time, s.MinAxs(1,1:len).*sArena.cF)
axis([0 time(end) 0 20])
xlabel('Time (s)')
ylabel('Length (mm)')
legend('Major Axis', 'Minor Axis')
 
subplot(222)
plot(time(1:len-1), s.Kinematics.yaw(1:len-1))
axis([0 time(end) -180 180])
xlabel('Time (s)')
ylabel('Yaw (degrees)')
 
subplot(223)
plot(s.Center.x(1:s.LightOn),s.Center.y(1:s.LightOn),'b')
hold on
plot(s.Center.x(s.LightOn+1:end),s.Center.y(s.LightOn+1:end),'r')
h = ellipsev2(sArena.rad, sArena.rad, 0, sArena.arenaCent(1), sArena.arenaCent(2));
plot(h.y,h.x,'b')
set(gca,'YDir','Reverse')
xlabel('Pixels')
ylabel('Pixels')
sgtitle(tempName)
hold off
 
subplot(224)
avgMjr = mean(s.MjrAxs(1,1:len).*sArena.cF);
avgMin = mean(s.MinAxs(1,1:len).*sArena.cF);
stdMjr = std(s.MjrAxs(1,1:len).*sArena.cF);
stdMin = std(s.MinAxs(1,1:len).*sArena.cF);
c = categorical({'majorAxis','minorAxis'});
hold on
bar(c, [avgMjr, avgMin])
e = errorbar(c,[avgMjr, avgMin],[stdMjr, stdMin],'.');
e.Color = 'red';
e.CapSize = 30;
sgtitle(tempName)
hold off
ylim([0 15])
legend(['Mean: ' num2str(avgMjr) ' ' num2str(avgMin)], ...
    ['StdA: ' num2str(stdMjr) ' ' num2str(stdMin)], 'location', 'best')
% print([foldername '/Basic Stats'],'-dpdf')
print('-painters','-bestfit', '-dpdf',[foldername '/Basic Stats - ' foldername(10:end)]);
 
end