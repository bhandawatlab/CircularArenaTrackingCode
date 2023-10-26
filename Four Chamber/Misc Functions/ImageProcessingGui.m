function [] = ImageProcessingGui()
% SpikeSortingGUI Comments/Documentations
close all
fileName = [];fileNameStim = [];DataExp = [];IdxInClusterUpdated = cell(1,4);

hs = addcomponents;
    function hs = addcomponents
        hs.fig = figure('Visible','off','Tag','fig','SizeChangedFcn',@resizeui);
        set(gcf,'position',[10 120 1000 800])
        hs.axThrust = axes('Parent',hs.fig,'Units','pixels', 'Tag','ax');
        hs.axXY = axes('Parent',hs.fig,'Units','pixels', 'Tag','ax');
        hs.axYaw = axes('Parent',hs.fig,'Units','pixels', 'Tag','ax');
        
        
        hs.btn = uicontrol(hs.fig,'String',...
            'Load Data',...
            'Callback',@loadData,...
            'Tag','button');
        
        hs.btn2 = uicontrol(hs.fig,'String',...
            'Update plots',...
            'Callback',@plotData,...
            'Tag','button');
        
        hs.btn3 = uicontrol(hs.fig,'String',...
            'Undo last action',...
            'Callback',@undo,...
            'Tag','button');
        
        hs.btn4 = uicontrol(hs.fig,'String',...
            'Save Data',...
            'Callback',@saveData,...
            'Tag','button');
        
        hs.txtbox1 = uicontrol(hs.fig,'Style','edit',...
            'String','Begin Position');
        
        hs.txtbox2 = uicontrol(hs.fig,'Style','edit',...
            'String','End Position');
        
        hs.pm = uicontrol(hs.fig,'Style','popupmenu',...
            'String',{'1','2','3','4'},...
            'Value',1);
        
    end

    function resizeui(hObject,event)
        % Get figure width and height
        figwidth = hs.fig.Position(3);
        figheight = hs.fig.Position(4);
        
        % Set button position
        bleftedge = 10;
        
        % Set axes position
        axheight1 = .2*figheight;
        axbottomedge = max(0,figheight - axheight1 - 30);
        axleftedge1 = bleftedge + 40;
        axwidth1 = max(0,(figwidth - 2*axleftedge1));
        hs.axThrust.Position = [axleftedge1 axbottomedge axwidth1 axheight1];
        
        % Set axes position S1
        axheight = .2*figheight;
        axbottomedge1 = max(0,figheight + axheight1-axbottomedge+30);
        axwidth = axwidth1./2-200;
        hs.axXY.Position = [axleftedge1 axbottomedge1 axwidth axheight];
        
        % Set axes position S2
        axheight = .2*figheight;
        axbottomedge1 = max(0,figheight + axheight1-axbottomedge+30);
        axwidth = axwidth1./2+170;
        axleftedge = bleftedge + (figwidth)/2-180;
        hs.axYaw.Position = [axleftedge axbottomedge1 axwidth axheight];
        
        axwidth = figwidth./3;
        
        y = axbottomedge1-60;
        x = axleftedge1;
        hs.btn.Position = [x y axwidth1./5 hs.btn.Position(4)];
        
        y = axbottomedge1-60;
        x = axleftedge1+axwidth/1.5;
        hs.btn2.Position = [x y axwidth1./5 hs.btn.Position(4)];
        
        y = axbottomedge1-120;
        x = axleftedge1;
        hs.btn3.Position = [x y axwidth1./5 hs.btn.Position(4)];
        
        y = axbottomedge1-120;
        x = axleftedge1+axwidth/1.5;
        hs.btn4.Position = [x y axwidth1./5 hs.btn.Position(4)];
        
        y = axbottomedge1-60;
        x = bleftedge + (figwidth)/2+20;
        hs.txtbox1.Position =  [x y axwidth hs.txtbox1.Position(4)];
        
        y = axbottomedge1-90;
        x = bleftedge + (figwidth)/2+20;
        hs.txtbox2.Position =  [x y axwidth hs.txtbox2.Position(4)];
        
    end

    function selectPoints(hObject,event)
        C{1} = get(hs.txtbox1,'String');
        C{2} = get(hs.txtbox2,'String');
        D = load(fileName{1});
        Data = D.s;
        Arena = D.sArena;
        fHandle = hs;
        PlottingFunctionsGUI(Data,Arena,fHandle);
        
        plotData(hObject,event);
        
        
        IdxInClusterUpdated = cell(1,4);
        for i = 1:4
            temp = strsplit(C{i},',');
            for j = 1:length(temp)
                temp2 = str2double(temp{j});
                if ~isnan(temp2)
                    IdxInClusterUpdated{i} = [IdxInClusterUpdated{i};IdxInCluster{temp2}];
                end
            end
        end
    end

    function loadData(hObject,event)
        [FileName,PathName] = uigetfile('*.mat','Select the Stim file',[pwd '\DataOld']);
        fileName = {[PathName FileName]};
    end

    function [] = plotData(hObject,event)
        D = load(fileName{1});
        if isfield(D, 's')
            Data = D.s;
            Arena = D.sArena;
            fHandle = hs;
            PlottingFunctionsGUI(Data,Arena,fHandle);
        end
    end

    function [] = updatePlot(hObject,event)
        D = load(fileName{1});
        if isfield(D, 's')
            Data = D.s;
            Arena = D.sArena;
            fHandle = hs;
            PlottingFunctionsGUI(Data,Arena,fHandle);
        end
    end

    function saveData(hObject,event)
        if ~isempty(DataExp)
            DataExp.IdxInClusterUpdated = IdxInClusterUpdated;
            fName = get(hs.txtbox5,'String');
            save('DataExp',fName)
        end
    end


hs.fig.Visible = 'on';
end




