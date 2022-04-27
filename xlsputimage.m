function xlsputimage(excelFileName, imageFileName, activeSheet, topLeftCorner, imageSize)
%XLSPUTIMAGE place an image inside an excel chart
%
%   XLSPUTIMAGE(EXCELFILENAME, IMAGEFILENAME) puts the image file IMAGEFILENAME  
%   inside the chart "Sheet1" at a default location. 
%
%   XLSPUTIMAGE(EXCELFILENAME, IMAGEFILENAME, ACTIVESHEET) puts the image file  
%   IMAGEFILENAME inside the chart "ACTIVESHEET" at a default location. 
%
%   XLSPUTIMAGE(EXCELFILENAME, IMAGEFILENAME, ACTIVESHEET, TOPLEFTCORNER) puts the  
%   image file IMAGEFILENAME inside the chart "ACTIVESHEET" where the top left 
%   corner is at the cell TOPLEFTCORNER.
%
%   XLSPUTIMAGE(EXCELFILENAME, IMAGEFILENAME, ACTIVESHEET, TOPLEFTCORNER) puts the 
%   image file IMAGEFILENAME inside the chart "ACTIVESHEET" where the top left 
%   corner is at the cell TOPLEFTCORNER and the image size is defined by the vector
%   IMAGESIZE (width and height respectively).
%
% references: 
%   https://stackoverflow.com/questions/7779880/insert-image-in-excel-through-matlab
%   https://docs.microsoft.com/en-us/office/vba/api/excel.shapes.addpicture
% 
% 
%  Written by Noam Greenboim
%
%% input assignment and assertion
if nargin < 5
    imInfo = imfinfo(imageFileName);
    imWidth = imInfo.Width;
    imHeight = imInfo.Height;
    if nargin < 4
        topLeftCorner = 'C9';
        if nargin < 3
            activeSheet = 'Sheet1';
            if nargin < 2
                error('excelFileName and imageFileName inputs must be specified')
            end
            if ~ischar(excelFileName) || isempty(excelFileName)
                error('excelFileName must be a string.')
            end
            if ~ischar(imageFileName) || isempty(imageFileName)
                error('imageFileName must be a string.')
            end
            if ~exist(excelFileName,'file')==2
                error(['The file ' imageFileName ' can''t be found'])
            end
            if ~exist(imageFileName,'file')==2
                error(['The file ' imageFileName ' can''t be found'])
            end
        else
            if ~ischar(activeSheet)
                error('Default sheet name must be a string.');
            end
        end
    else
        if ~ischar(topLeftCorner) || ~(regexpi(topLeftCorner,'^[a-zA-Z]{1,3}\d{1,8}$')==1)
            error('topLeftCorner must be a valid Excel cell. For example: A9, D15, etc.')
        end
    end
else
    if isnumeric(imageSize) && length(imageSize)==2
        imWidth = imageSize(1);
        imHeight = imageSize(2);
    else
        error('imageSize must be a 2-element vector, of width and height respectively')
    end
end
%%
try
    objExcel = actxserver('Excel.Application');
    objExcel.Workbooks.Open(excelFileName); % Open Excel file. Full path is necessary!
catch
    terminateSession(objExcel)
    error(['Unable to open the excel file ' excelFileName])
end
try
    objExcel.ActiveWorkbook.Worksheets.Item(activeSheet).Activate; % set active sheet
catch
    objExcel.ActiveWorkbook.Worksheets.Add;
    objExcel.ActiveWorkbook.ActiveSheet.Name = activeSheet;
end
LinkToFile = 0;
SaveWithDocument = 1;
left = objExcel.ActiveSheet.Range(topLeftCorner).Left;
top = objExcel.ActiveSheet.Range(topLeftCorner).Top;
try
    objExcel.ActiveSheet.Shapes.AddPicture(imageFileName,LinkToFile,SaveWithDocument,left,top,imWidth,imHeight);
catch
    terminateSession(objExcel)
    error(['Unable to place the image ' imageFileName ' in the Excel sheet'])
end
% Save, close and clean up.
objExcel.ActiveWorkbook.Save;
terminateSession(objExcel)
function terminateSession(objExcel)
objExcel.ActiveWorkbook.Close;
objExcel.Quit;
objExcel.delete;
