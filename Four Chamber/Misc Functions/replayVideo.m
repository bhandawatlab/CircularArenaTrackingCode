function [] = replayVideo(fly,Inside,filename,badFlies)
close all
[readframe,len,~,vid] = get_readframe_fcn(filename);

v = VideoWriter([filename(1:end-7) '_Processed.avi'],'Motion JPEG AVI');
v.Quality = 100;
v.FrameRate = 10;
open(v)
regions = 4;
h = figure(1);

%Use this code to replay back the fly video
set(gcf, 'Position', [400, 50, 1200, 915])

r=vid.max_height;
c=vid.max_width;

circleImage = false(r,c,regions);
for k = 1:regions
    I = Inside(:,:,k);
    stats = regionprops(I, 'EquivDiameter', 'Centroid');
    circles{k}.xy = stats.Centroid;
    circles{k}.r = stats.EquivDiameter./2;
    
    [x, y] = meshgrid(1:c, 1:r);
    circleImage(:,:,k) = (x - stats.Centroid(1)).^2 + (y - stats.Centroid(2)).^2 <= (stats.EquivDiameter./(2.*4/1.25)).^2;
end



for i = 1:len
    I =readframe(i);
    I2 = mean(double(I),3);
    I2(~any(Inside,3)) = 200;
    I5 = I2;
    I5(any(circleImage,3)) = 50;
    
    I3 = zeros(2048,2048,3);
    I3(:,:,1) = I5;
    I3(:,:,2) = I2;
    I3(:,:,3) = I2;
    imshow(uint8(I3),'InitialMagnification','fit');hold on
    for k = 1:regions
        if any(k == badFlies)
        
        else
            % major/minor axis
            MjrAxs = fly{k}.MjrAxs(i)./2;
            MinAxs = fly{k}.MinAxs(i)./2;

            %centroid
            centroid(1) = fly{k}.Center.x(i);
            centroid(2) = fly{k}.Center.y(i);

            %orientation (take orientation associated with largest major axis value)
            ang = fly{k}.ang(1,i);

            %calculate major and minor axis endpoints
            xMajor1 = centroid(1) + MjrAxs * cosd(ang);
            xMajor2 = centroid(1) - MjrAxs * cosd(ang);
            yMajor1 = centroid(2) + MjrAxs * sind(ang);
            yMajor2 = centroid(2) - MjrAxs * sind(ang);

            xMinor1 = centroid(1) + MinAxs * cosd(ang+90);
            xMinor2 = centroid(1) - MinAxs * cosd(ang+90);
            yMinor1 = centroid(2) + MinAxs * sind(ang+90);
            yMinor2 = centroid(2) - MinAxs * sind(ang+90);


            %calculate the pixels composing of the oval surrounding body
            body = ellipsev2(MjrAxs, MinAxs, (ang)*pi/180, centroid(1), centroid(2),'r',600);
            body.x(body.x<1) = 1;
            body.x(body.x>r) = r;
            body.y(body.y<1) = 1;
            body.x(body.x>c) = c;

            %Draw major and minor axis lines based on pixels
            [cx1,cy1,~] = improfile(I2,[xMajor1,xMajor2], [yMajor1,yMajor2],25);
            [cx2,cy2,~] = improfile(I2,[xMinor1,xMinor2], [yMinor1,yMinor2],25);

            % plot fly major axis, minor axis, body, and head
    %         plot(cx1,cy1,'color',[0.5 0.5 0.5],'linewidth',1);
    %         plot(cx2,cy2,'color',[0.5 0.5 0.5],'linewidth',1);
            plot(body.y,body.x,'color',[0 0.7 0],'linewidth',1);
            viscircles([xMajor1 yMajor1],1,'Color','r');
        end
    end
    
    % Add timestamp
    m = floor(i./(30.*60)); s = mod(floor(i./30),60);ms = round(mod(i,30).*1000./30);
    text(r./2-120,c./2,['Frame #: ' num2str(i)],'FontSize',10,'FontWeight','Bold')
    text(r./2-120,c./2+50,['Time (s): ' num2str(m) ':' num2str(s,'%02.f') ':' num2str(ms,'%03.f')],'FontSize',10,'FontWeight','Bold')
    
    F = getframe(h);
    writeVideo(v,F);
    hold off

    %used to show progress
    if mod(i,1000) == 0
        i
    end
end
close(v)
end