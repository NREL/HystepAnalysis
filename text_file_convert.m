%% Initialize variables.
[f, p] = uigetfile({'*.txt'},'Pick a Dispenser File',path,'MultiSelect', 'off');
savefile = strrep(f,'txt','mat');
savename = [p savefile];
delimiter = ',';
startRow = 2;
%% Format string for each line of text:
%   column1: text (%s)
%	column2: text (%s)
%   column3: text (%s)
%	column4: double (%f)
%   column5: text (%s)
%	column6: text (%s)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%s%s%s%f%s%s%[^\n\r]';

%% Open the text file.
fileID = fopen([p f],'r');

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false);

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Allocate imported array to column variable names
Date = dataArray{:, 1};
time = dataArray{:, 2};
Description = dataArray{:, 3};
Value = dataArray{:, 4};
Units = dataArray{:, 5};
Identifier = dataArray{:, 6};

%% Clear temporary variables
%clearvars filename delimiter startRow formatSpec fileID dataArray ans;
% %
% Remove NaN values
nanind = ~isnan(Value);
Date = Date(nanind);
Description = Description(nanind);
Identifier = Identifier(nanind);
time = time(nanind);
Units = Units(nanind);
Value = Value(nanind);

try
    datetime = datenum(cell2mat(Date)) + datenum(cell2mat(time)) - datenum('00:00', 'HH:MM');
catch
end
% x = [Category Date];

%%
% Sort into column matrices
%     category_names = unique(Category);
%     num_category = length(category_names);
%     M = {};
%     for k = 1:num_category
%         temp_ind = strcmpi(Category,category_names(k));
%         temp_cat = Category(temp_ind);
% %         temp_dat = Date(temp_ind);
%         temp_des = Description(temp_ind);
%         temp_ide = Identifier(temp_ind);
%         temp_tim = datetime(temp_ind);
%         temp_uni = Units(temp_ind);
%         temp_val = Value(temp_ind);

description_names = unique(Description);
num_identifier = length(description_names);
D = {};
for j = 1:num_identifier
    index = strcmpi(Description, description_names(j));
    name = strrep(cell2mat(description_names(j)),' ','');
    name = strrep(name,'#','');
    name = strrep(name,')','');
    name = strrep(name,'(','');
    name = strrep(name,'[','');
    name = strrep(name,']','');
    name = strrep(name,'-','');
    name = strrep(name,'\','');
    name = strrep(name,'/','');
    nametime = [name 'Time'];
    D.(name) = Value(index);
    D.(nametime) = datetime(index);
end
%     end
%save(savename, 'D')
