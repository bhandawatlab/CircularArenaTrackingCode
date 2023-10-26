function pdfmerge
%  inputs=uigetfile('*.pdf','Select the INPUT pdf FILE to merge','MultiSelect','on');
 [inputs, path2]=uigetfile('*.pdf','Select the INPUT pdf FILE to merge','MultiSelect','on')
 input_p =[]; 
 for i=1:length(inputs)
     input_p{i} =[path2 inputs{i}]
 end
[file,path] = uiputfile('*.pdf', 'Select output filename');
outputfile = [path file]; 
% following code is used from the reference 
% Benjamin Großmann (2021). Merge PDF-Documents (https://www.mathworks.com/matlabcentral/fileexchange/89127-merge-pdf-documents), MATLAB Central File Exchange. Retrieved August 24, 2021.
memSet = org.apache.pdfbox.io.MemoryUsageSetting.setupMainMemoryOnly();
merger = org.apache.pdfbox.multipdf.PDFMergerUtility;
cellfun(@(f) merger.addSource(f), input_p)
merger.setDestinationFileName(outputfile)
merger.mergeDocuments(memSet)
