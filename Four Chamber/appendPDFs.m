clear all;
close all;
clc;

outputFileName = 'Consolidated_Fly_Paths.pdf';
d ='PDFs';
files = dir(fullfile(d, '*.pdf'));
% Display the names
%disp(files.name)
%inputFiles = [];
inputFiles  = struct2cell(files); %vertcat(files.name);
%inputFiles.append(files.name)
input = [];
for i=1:size(files)
    input = [input, 'PDFs/'+string(inputFiles{1,i})];
end
mergePdfs(input, outputFileName);

