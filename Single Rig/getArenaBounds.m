function [Inside,vertices] = getArenaBounds(Inside,vertices,filename,regions)
[readframe,~,~,headerinfo] = get_readframe_fcn(filename);
img = readframe(1);

if isempty(Inside)
    [Inside,vertices] = collectArenaBounds(img,headerinfo,regions);
    save([filename(1:end-5) '_Inside.mat'],'Inside','vertices')
else
    if check == true
        %plot original video frame 1 and circles around arena
        imshow(uint8(img),'InitialMagnification','fit');hold on
        set(gcf,'Position',[2 42 958 954])
        for k = 1:regions
            I = Inside(:,:,k);
            stats = regionprops(I, 'EquivDiameter', 'Centroid');
            [xPts,yPts] = circle(stats.Centroid(1),stats.Centroid(2),stats.EquivDiameter./2);
            plot(xPts,yPts,'r','linestyle','--','LineWidth',1)
        end
        
        % ask the user if they accept these bounds. If not, they have to
        % manually find bounds
        x = input('Accept these arena bounds? [Y/N]');
        if strcmpi(x,'Y')
            
        elseif strcmpi(x,'N')
            [Inside,vertices] = collectArenaBounds(img,headerinfo,regions);
        else
            display('Not a valid option, assume Y')
        end
    end
end


figure;set(gcf,'position',[2 42 958 954])
imagesc(img)
hold on;
for i = 1:regions
    plot(vertices{i}(:,1),vertices{i}(:,2),'r','Linewidth',5);
    text(mean(vertices{i}(:,1)),mean(vertices{i}(:,2)),num2str(i),'fontsize',20)
end
print('-painters','-dpdf',['Analysis/' filename(8:end-20) '_Arena.pdf'])


end

function [Inside,vertices] = collectArenaBounds(img,v,regions)

figure;set(gcf,'Position',[2 42 958 954])
imagesc(img);
Inside = zeros(v.max_width,v.max_height,4);
for i = 1:regions
    h=imellipse; %%%%%%HPARENT IS NOT DEFINED imellipse(hparent,position)
    vertices{i}=wait(h);
    Inside(:,:,i) = h.createMask();
end
Inside = uint8(Inside);

end