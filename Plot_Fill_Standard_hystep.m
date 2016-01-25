%% Written by T. Luetje on 12/22/2014
%Browse for File
%Convert File
%Get right Table
%Import Table
%Extract signals of converted File
%Plot Fill
%Calculate Results
%Show Results

function Plot_Fill_Standard_hystep(handles)

[filename, pathname] = uigetfile({'*.xlsx','xlsx-Files (*.xlsx)';'*.xls','xls-Files (*.xls)';'*.xlcs','xlcs-Files (*.xlcs)';'*.*','All Files'},....
    'Pick a Station Fill File',path,'MultiSelect', 'off');
if ischar(filename)==1
    [status, sheets] = xlsfinfo([pathname filename]);
    if ~isempty(status)
        [overviewData, overviewDataTxt, overviewDataRaw] = xlsread([pathname filename], sheets{1});
        numsheets = numel(sheets);
        numtests = numsheets-1;
        for i = 1:numtests
            testnames{i} = overviewDataTxt{4+i,1};
        end
        testnum = 1;
        while testnum ~= 0
            testnum = menu('Pick the test you want to process',testnames);
            if testnum ~= 0
                Plot_Fill_Standard_OneFill_hystep(handles, filename, pathname,testnum+1)
            end
        end
%         for i = 1:numsheets-1
%             reply = input(['Process test ' num2str(i) ', ' overviewDataTxt{4+i,16} '? Y/N [Y]:'],'s');
%             if isempty(reply)
%                 reply = 'Y';
%             end
%             if reply == 'Y'
%                 Plot_Fill_Standard_OneFill_hystep(handles, filename, pathname,i+1)
%             end
%         end
    end
end


