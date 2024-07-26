function Table = tableprint(Data,varargin)
%TABLEPRINT Print data in tabular form
%
% Syntax: 
%
%   tableprint(Data)
%   tableprint(Data,param,val,...)
%   Table = tableprint(...)
%
% Description:
%
%   TABLEPRINT prints a data table in a specified format.
%
%   Formatting options include the ability to:
%
%   o Use mixed numerical/character data
%   o Set numerical precision
%   o Convert large or small numbers to scientific format
%   o Align numbers to the left, right, or at the decimal point
%   o Customize the border and inner separator styles
%   o Adjust the column width
%   o Customize row and column indexing
%   o Leave empty entries/blocks in the table (look like merge and split)
%
% Input Arguments:
%
%   Data - m-by-n numeric matrix, cell array of heterogeneous data, 
%          string array, table or dataset.
%
% Optional Input Parameter Name/Value Pairs:
%
%   NAME                    VALUE
%
%   'title'                 String scalar or character vector to be used as
%                           the table title. 
%                           The default is '' (empty).
%
%   'NWStr'                 String scalar or character vector to be used in
%                           the northwest corner of the table.
%                           The default is '' (empty).
%
%   'colNames'              String vector or cell array of character
%                           vectors to be used as column names; or string
%                           scalar or character vector to be used as a
%                           common name to which indices (1),(2),... will
%                           be appended. The default is {} (empty).
%
%   'rowNames'              String vector or cell array of character
%                           vectors to be used as row names; or string
%                           scalar or character vector to be used as a
%                           common name to which indices (1),(2),... will
%                           be appended. The default is {} (empty).
%
%   'indexStyle'            Style for indexed header names. Values are
%                           '(1)', '[1]', '{1}', '.1', '1', 'none'. The
%                           default is '(1)'.
%
%   'indexStart'            Starting number for indexed header names. The
%                           default is 1.
%
%   'indexOrder'            Order of indexed header names. Values are
%                           'ascend' or 'descend'. The default is 'ascend'.
%
%   'numDigits'             A scalar that controls display precision
%                           (number of digits after the decimal point) for
%                           the entire table, or a vector that controls
%                           display precision for each column of the table.
%                           The default is 4.
%
%   'maxCutoff'             A scalar that controls the threshold that
%                           converts large numbers to the scientific format
%                           for the entire table, or a vector that controls
%                           the threshold for each column of the table. 
%                           The default is 1e5.                           
%
%   'minCutoff'             A scalar that controls the threshold that
%                           converts small numbers to the scientific format
%                           for the entire table, or a vector that controls
%                           the threshold for each column of the table.
%                           The default is -Inf.
%
%   'omitNaN'               Logical value indicating whether to skip NaN
%                           and display an empty entry. This is useful for
%                           making empty blocks of the table that looks 
%                           like merge and split. The default is false
%                           (print NaN as is).
%
%   'bracketData'           A vector of row numbers. In those rows, numeric
%                           numbers will be inside the (). This is useful
%                           in a statistic table that displays standard
%                           errors in parentheses. The default is empty.
%
%   'horzSep'               Symbol used repeatedly to separate column
%                           headers from the data. The default is '-'.
%
%   'vertSep'               Symbol used repeatedly to separate row headers
%                           from the data. The default is '|'.
%
%   'borderTop'             Symbol used repeatedly for the top border of
%                           the table. The default is '' (empty).
%
%   'borderBottom'          Symbol used repeatedly for the bottom border of
%                           the table. The default is '' (empty).
%
%   'borderLeft'            Symbol used repeatedly for the left border of
%                           the table. The default is '' (empty).
%
%   'borderRight'           Symbol used repeatedly for the right border of
%                           the table. The default is '' (empty).
%
%   'borderInnerHorz'       Symbol or m-by-1 cell array of Symbols for
%                           inner horizontal border of the table. The
%                           default is '' (empty).
%
%   'borderInnerVert'       Symbol or 1-by-n cell array of Symbols for
%                           inner vertical border of the table. The default
%                           is '' (empty).
%
%   'borderMaster'          Logical value indicating whether or not to add
%                           a table border. A value of false turns off all
%                           other border parameters. The default is true.
%
%   'alignNumber'           Align numbers in each column to the 'left',
%                           'right', or at 'decimal'. If columns need
%                           to be aligned differently, use a cell array of
%                           those values. The default is 'decimal'.
%
%   'alignRowName'          Align row names to the 'left', 'right', or
%                           'center'. The default is 'left'.
%
%   'alignColName'          Align column names to the 'left', 'right', or
%                           'center'. The default is 'center'.
%
%   'colPaddingLeft'        Number of leading white spaces in each column.
%                           The default is 1.
%
%   'colPaddingRight'       Number of trailing white spaces in each column.
%                           The default is 1.
%
%   'fileID'                Scalar file ID indicating where to write the
%                           table. The default is 1, which prints to the
%                           command window.
%
%   'fileName'              file name, of the form 'myTable.xlsx',
%                           to write to an EXCEL spreadsheet. The default
%                           is '' (empty).
%
% Output Arguments:
%
%   Table - Cell array containing both the data and the header information.

% Parse inputs and set defaults:

callerName = 'tableprint';
parseObj = inputParser;
parseObj.addParameter('title','',@(x)validateattributes(x,{'char','string'},{'scalartext'},callerName));
parseObj.addParameter('NWStr','',@(x)validateattributes(x,{'char','string'},{'scalartext'},callerName)); 
parseObj.addParameter('colNames',{});
parseObj.addParameter('rowNames',{});
parseObj.addParameter('indexStyle','(1)',@(x)validateattributes(x,{'char','string'},{'scalartext'},callerName));
parseObj.addParameter('indexStart',1,@isscalar);
parseObj.addParameter('indexOrder','ascend',@(x)validateattributes(x,{'char','string'},{'scalartext'},callerName));
parseObj.addParameter('numDigits',4,@(x)validateattributes(x,{'numeric','logical'},{'2d','nonnegative'},callerName));
parseObj.addParameter('maxCutoff',1e5,@isscalar);
parseObj.addParameter('minCutoff',-Inf,@isscalar);
parseObj.addParameter('omitNaN',false,@(x)validateattributes(x,{'numeric','logical'},{'scalar','binary'},callerName));
parseObj.addParameter('bracketData',[],@(x)validateattributes(x,{'numeric'},{'2d'},callerName));
parseObj.addParameter('horzSep','-');
parseObj.addParameter('vertSep','|');
parseObj.addParameter('borderTop','');
parseObj.addParameter('borderBottom','');
parseObj.addParameter('borderLeft','');
parseObj.addParameter('borderRight','');
parseObj.addParameter('borderInnerHorz',''); 
parseObj.addParameter('borderInnerVert','');
parseObj.addParameter('borderMaster',true,@(x)validateattributes(x,{'numeric','logical'},{'scalar','binary'},callerName));
parseObj.addParameter('alignNumber','decimal',@(x)validateattributes(x,{'char','string','cell'},{},callerName));
parseObj.addParameter('alignRowName','left',@(x)validateattributes(x,{'char','string'},{'scalartext'},callerName));
parseObj.addParameter('alignColName','center',@(x)validateattributes(x,{'char','string'},{'scalartext'},callerName));
parseObj.addParameter('colPaddingLeft',1,@isscalar);
parseObj.addParameter('colPaddingRight',1,@isscalar);
parseObj.addParameter('fileID',1);
parseObj.addParameter('fileName','',@(x)validateattributes(x,{'char','string'},{'scalartext'},callerName))

parseObj.parse(varargin{:});

title = parseObj.Results.title;
NWStr = parseObj.Results.NWStr;
colNames = parseObj.Results.colNames;
rowNames = parseObj.Results.rowNames;
indexStyle = parseObj.Results.indexStyle;
indexStart = parseObj.Results.indexStart;
indexOrder = parseObj.Results.indexOrder;
numDigits = parseObj.Results.numDigits;
maxCutoff = parseObj.Results.maxCutoff;
minCutoff = parseObj.Results.minCutoff;
omitNaN = parseObj.Results.omitNaN;
bracketData = parseObj.Results.bracketData;
horzSep = parseObj.Results.horzSep;
vertSep = parseObj.Results.vertSep;
borderTop = parseObj.Results.borderTop;
borderBottom = parseObj.Results.borderBottom;
borderLeft = parseObj.Results.borderLeft;
borderRight = parseObj.Results.borderRight;
borderInnerHorz = parseObj.Results.borderInnerHorz;
borderInnerVert = parseObj.Results.borderInnerVert;
borderMaster = parseObj.Results.borderMaster;
alignNumber = parseObj.Results.alignNumber;
alignRowName = parseObj.Results.alignRowName;
alignColName = parseObj.Results.alignColName;
colPaddingLeft = parseObj.Results.colPaddingLeft;
colPaddingRight = parseObj.Results.colPaddingRight;
fileID = parseObj.Results.fileID;
fileName = parseObj.Results.fileName;

title = convert2char(title);
NWStr = convert2char(NWStr);
colNames = convert2char(colNames);
rowNames = convert2char(rowNames);
indexStyle = convert2char(indexStyle);
indexOrder = convert2char(indexOrder);
horzSep = convert2char(horzSep);
vertSep = convert2char(vertSep);
borderTop = convert2char(borderTop);
borderBottom = convert2char(borderBottom);
borderLeft = convert2char(borderLeft);
borderRight = convert2char(borderRight);  
borderInnerHorz = convert2char(borderInnerHorz);
borderInnerVert = convert2char(borderInnerVert); 
alignNumber = convert2char(alignNumber);
alignRowName = convert2char(alignRowName); 
alignColName = convert2char(alignColName); 
fileName = convert2char(fileName); 
       
% If borderMaster is false, turn off all border settings:

if ~borderMaster
    horzSep = [];
    vertSep = [];
    borderTop = [];
    borderBottom = [];
    borderLeft = [];
    borderRight = [];
    borderInnerHorz = [];
    borderInnerVert = [];
end

% If there are no column/row headers, then there are no horizontal/vertical
% separators, except in the presence of NWStr:

if isempty(colNames) && isempty(NWStr)
    horzSep = [];    
end

if isempty(rowNames) && isempty(NWStr)
    vertSep = [];    
end

% All input types are converted to cell arrays for processing.

% Covert table objects:

if isa(Data,'table')
    
    % Extract data:
    
    TB = Data;
    Data = table2cell(TB);
    colNamesOverride = isempty(colNames);
    rowNamesOverride = isempty(rowNames);
    
    % Import column names from the table object:
    
    if colNamesOverride
        colNames = TB.Properties.VariableNames; % Always non-empty.
    end
    
    % Import row names from the table object:
    
    if rowNamesOverride && ~isempty(TB.Properties.RowNames)
        rowNames = TB.Properties.RowNames;
    end
    
    % Import variable descriptions from the table object: 
    
    if colNamesOverride && ~isempty(TB.Properties.VariableDescriptions)        
        temp = TB.Properties.VariableDescriptions;
        
        for m = 1:size(TB,2)
            
            nwords = length(temp{1,m});
            
            % Break long descriptions into two rows for better display:
            
            if (nwords > 20) && (nwords > 2*length(colNames{m}))
                
                cut = round(nwords/2);                
                loc = find(isspace(temp{1,m}));
                
                if ~isempty(loc)
                    
                    [~,ind] = min(abs(loc-cut));
                    cut = loc(ind);
                    temp{2,m} = temp{1,m}(cut+1:end);
                    temp{1,m} = temp{1,m}(1:cut);
                    
                else
                    
                    temp{2,m} = temp{1,m}(cut+1:end);
                    temp{1,m} = [temp{1,m}(1:cut),'-'];
                    
                end
                
            end
            
        end
        
        colNames = [colNames;temp];
        
    end
    
    % Import units from the table object:
    
    if colNamesOverride && ~isempty(TB.Properties.VariableUnits)
        colNames = [colNames;TB.Properties.VariableUnits];
    end
    
    % Import title from the table object:
    
    if isempty(title) && ~isempty(TB.Properties.Description)
        title = TB.Properties.Description;
    end
    
    % Import DimensionNames from the table object:
    
    if isempty(NWStr) && (~strcmp(TB.Properties.DimensionNames{1},'Row')...
            || ~strcmp(TB.Properties.DimensionNames{2},'Variable'))
        NWStr = [TB.Properties.DimensionNames{1},'\',TB.Properties.DimensionNames{2}];
    end
    
    % Check multiple column data:
    
    for m = size(Data,2):-1:1
        
        dim2 = size(Data{1,m},2);
        
        if dim2 > 1 && isnumeric(Data{1,m})
            
            Data = [Data(:,1:m-1),num2cell(table2array(TB(:,m))),Data(:,m+1:end)];
            colNames = [colNames(:,1:m-1),colNames(:,m*ones(1,dim2)),colNames(:,m+1:end)];
            
        end
        
    end
    
end

% Covert dataset arrays:

if isa(Data,'dataset')
    
    % Extract data:
    
    DS = Data;
    Data = dataset2cell(DS);
    
    if ~isempty(DS.Properties.VarNames)
        Data = Data(2:end,:);
    end
    
    if ~isempty(DS.Properties.ObsNames)
        Data = Data(:,2:end);
    end
    
    colNamesOverride = isempty(colNames);
    rowNamesOverride = isempty(rowNames);    
    
    % Import column names from the dataset object:
    
    if colNamesOverride && ~isempty(DS.Properties.VarNames)
        colNames = DS.Properties.VarNames;
    end
    
    % Import ObsNames from the dataset object:
    
    if rowNamesOverride && ~isempty(DS.Properties.ObsNames)
        rowNames = DS.Properties.ObsNames;
    end
    
    % Import variable descriptions from the dataset object:
    
    if colNamesOverride && ~isempty(DS.Properties.VarDescription)
        
        temp = DS.Properties.VarDescription;
        
        for m = 1:length(temp)
            
            nwords = length(temp{1,m});
            
            % Break long descriptions into two rows for better display:
            
            if (nwords > 20) && (nwords > 2*length(colNames{m}))
                
                cut = round(nwords/2);                
                loc = find(isspace(temp{1,m}));
                
                if ~isempty(loc)
                    
                    [~,ind] = min(abs(loc-cut));
                    cut = loc(ind);
                    temp{2,m} = temp{1,m}(cut+1:end);
                    temp{1,m} = temp{1,m}(1:cut);
                    
                else
                    
                    temp{2,m} = temp{1,m}(cut+1:end);
                    temp{1,m} = [temp{1,m}(1:cut),'-'];
                    
                end
                
            end
            
        end
        
        colNames = [colNames;temp];
        
    end
    
    % Import units from the dataset object:
    
    if colNamesOverride && ~isempty(DS.Properties.Units)
        colNames = [colNames;DS.Properties.Units];
    end
    
    % Import title from the dataset object:
    
    if isempty(title) && ~isempty(DS.Properties.Description)
        title = DS.Properties.Description;
    end
    
    % Import DimNames from the dataset object:
    
    if isempty(NWStr) && (~strcmp(DS.Properties.DimNames{1},'Observations')...
            || ~strcmp(DS.Properties.DimNames{2},'Variables'))
        NWStr = [DS.Properties.DimNames{1},'\',DS.Properties.DimNames{2}];
    end   
    
end

% Covert financial time series objects:

if isa(Data,'fints')
    
    Obj = Data;
    fields = ftsinfo(Obj);
    colNames = fields.seriesnames;
    nseries = fields.nseries;
    nobs = fields.ndata;
    Data = zeros(nobs,nseries);
    
    for m = 1:nseries
        Data(:,m) = getfield(Obj,colNames{m}); %#ok<GFLD>
    end
    
    rowNames = num2cell(datestr(Obj.dates),2);
    
    try %#ok<TRYNC>
        times = num2cell(datestr(Obj.times),2);
        rowNames = [rowNames,times];
    end
    
    NWStr = ' Dates';
    
    if ~isempty(Obj.desc)
        title = Obj.desc;
    end
        
end

% Covert matrix inputs to cell arrays:

if isnumeric(Data)
    Data = num2cell(Data);
end

if isstring(Data)
    Data = convert2char(Data);
end

if ~iscell(Data)
    error('Input data must be a numeric matrix, cell array, table or dataset array.')
end

% Limit display of very large data sets:

[numRows,numCols] = size(Data);

if numRows > 1e4
    
    warning('There are more than 10000 observations. Only the first and last 100 observations will be displayed.')
    
    Data = [Data(1:100,:);Data(end-99:end,:)];
    
    if iscell(rowNames) && size(rowNames,1)>200
        rowNames = [rowNames(1:100,:);rowNames(end-99:end,:)];
    end
    
    numRows = 200;
    
end

% Flags to print column/row headers:

if isempty(colNames) && isempty(NWStr)    
    printColumnHeader = false;
else
    printColumnHeader = true;
end

if isempty(rowNames) && isempty(NWStr)    
    printRowHeader = false;
else
    printRowHeader = true;
end

% Default (empty) column/row names:

if isempty(colNames)
    colNames = cell(1,numCols);
end

if isempty(rowNames)
    rowNames = cell(numRows,1);    
end

% Process multi-symbol column/row names:

if iscell(colNames) % Not a single symbol
    
    [colDim1,colDim2] = size(colNames);
    
    if colDim1 > 1 && colDim2 == 1
        colNames = colNames';
        colDim2 = colDim1;
    end
    
    if colDim2 < numCols
        error('Each column must have a name.')
    end
    
end

if iscell(rowNames) % Not a single symbol
    
    [rowDim1,rowDim2] = size(rowNames);
    
    if rowDim1 == 1 && rowDim2 > 1
        rowNames = rowNames';
        rowDim1 = rowDim2;
    end
    
    if rowDim1 < numRows
        error('Each row must have a name.')
    end
    
end
 
% Process single-symbol row names:

if ~iscell(rowNames) % Single symbol
    
    rowNamesStr = rowNames;
    rowNames = cell(numRows,1);
    rowNames(:) = {rowNamesStr};
    indexOrder = ~strcmpi(indexOrder,'descend')*2 - 1;
    
    switch indexStyle
        
        case '(1)'
            
            for m = 1:numRows
                rowNames{m} = [rowNames{m},'(',num2str(indexStart + indexOrder*(m-1)),')'];
            end
            
        case '[1]'
            
            for m = 1:numRows
                rowNames{m} = [rowNames{m},'[',num2str(indexStart + indexOrder*(m-1)),']'];
            end
            
        case '{1}'
            
            for m = 1:numRows
                rowNames{m} = [rowNames{m},'{',num2str(indexStart + indexOrder*(m-1)),'}'];
            end
            
        case '.1'
            
            for m = 1:numRows
                rowNames{m} = [rowNames{m},'.',num2str(indexStart + indexOrder*(m-1))];
            end
            
        case '1'
            
            for m = 1:numRows
                rowNames{m} = [rowNames{m},num2str(indexStart + indexOrder*(m-1))];
            end
            
    end    
end

% Process single-string column names:

if ~iscell(colNames) % Single string
    
    colNamesStr = colNames;
    colNames = cell(1,numCols);
    colNames(:) = {colNamesStr};
    indexOrder = ~strcmpi(indexOrder,'descend')*2 - 1;
    
    switch indexStyle
        
        case '(1)'
            
            for m = 1:numCols
                colNames{m} = [colNames{m},'(',num2str(indexStart + indexOrder*(m-1)),')'];
            end
            
        case '[1]'
            
            for m = 1:numCols
                colNames{m} = [colNames{m},'[',num2str(indexStart + indexOrder*(m-1)),']'];
            end
            
        case '{1}'
            
            for m = 1:numCols
                colNames{m} = [colNames{m},'{',num2str(indexStart + indexOrder*(m-1)),'}'];
            end
            
        case '.1'
            
            for m = 1:numCols
                colNames{m} = [colNames{m},'.',num2str(indexStart + indexOrder*(m-1))];
            end
            
        case '1'
            
            for m = 1:numCols
                colNames{m} = [colNames{m},num2str(indexStart + indexOrder*(m-1))];
            end
            
    end
    
end

% Number of digits, max and min cutoff points:

if isscalar(numDigits)
    numDigits = numDigits(1,ones(1,numCols));
end

if isscalar(maxCutoff)
    maxCutoff = maxCutoff(1,ones(1,numCols));
end

if isscalar(minCutoff)
    minCutoff = minCutoff(1,ones(1,numCols));
end

if length(numDigits) ~= numCols
    error('Length of numDigits must be equal to column numbers.')
end

if length(maxCutoff) ~= numCols
    error('Length of maxCutoff must be equal to column numbers.')
end

if length(minCutoff) ~= numCols
    error('Length of minCutoff must be equal to column numbers.')
end

% All numerical input is converted to character data for processing.

Contents = cell(numRows,numCols);

% Create display elements:

for n = 1:numCols
    
    for m = 1:numRows
        
        info = Data{m,n};
        
        if isempty(info)
            
            str = '';
            
        elseif isnumeric(info) || islogical(info)
            
            if ~isscalar(info)
                try
                    str = num2str(info);
                    Contents{m,n} = str;
                    continue
                catch
                    error('Data elements must be scalar.')
                end
            end
                        
            if info == round(info) && all(numDigits == numDigits(1)) && abs(info) < maxCutoff(n)
                str = num2str(info);
            elseif abs(info) > maxCutoff(n) || (abs(info) < minCutoff(n) && info ~= 0)
                str = num2str(info,['%10.',num2str(numDigits(n)),'e']);
            else
                str = num2str(info,['%10.',num2str(numDigits(n)),'f']);
            end
            
            if omitNaN && isnan(info)
                str = '';
            end
            
            if ~isempty(bracketData) && ismember(m,bracketData) && ~isempty(str)
                str = strcat('(', str, ')');
            end
            
        elseif ischar(info)
            
            str = info;
            
        elseif isstring(info)
            
            if isscalar(info)
                str = char(info);
            else
                str = strcat(info{:});
            end
            
        elseif isa(info,'nominal')
            
            str = char(info);
            
        elseif iscell(info)
            
            try
                str = num2str([info{:}]);
            catch
                error('Data element must be a scalar, string or character vector.')
            end
            
        else
            
            try
                str = char(info);
            catch
                error('Data element must be a scalar, string or character vector.')
            end
            
        end
        
        Contents{m,n} = str;
        
    end
    
end
        
% If numbers are aligned at the decimal, add spaces to align the padded
% values at the left:

if ~iscell(alignNumber)
    alignNumberCell = cell(1,numCols);
    alignNumberCell(:) = {alignNumber};
    alignNumber = alignNumberCell;
else
    if numel(alignNumber) ~= numCols
        error('Length of alignNumber must be equal to the number of columns.')
    end
end

for n = 1:numCols
    
    if strcmpi(alignNumber{n},'decimal')
        
        integerSize = zeros(numRows,1);
        
        for m = 1:numRows
            
            decPt = strfind(Contents{m,n},'.');
            
            if isempty(decPt)
                integerSize(m) = 0;
            else
                integerSize(m) = decPt(1)-1;
            end
            
        end
        
        integerSizeMax = max(integerSize);
        
        for m = 1:numRows
            
            integerPad = repmat(' ',1,integerSizeMax - integerSize(m));
            Contents{m,n} = [integerPad,Contents{m,n}];
            
        end
        
    end
      
end
    
% Compute column sizes:

ContentSize = zeros(numRows,numCols);

for n = 1:numCols
    for m = 1:numRows
        ContentSize(m,n) = length(Contents{m,n});
    end
end

[numRowsColHeader,numColsColHeader] = size(colNames);
ColHeaderSize = zeros(numRowsColHeader,numColsColHeader);

for n = 1:numColsColHeader
    for m = 1:numRowsColHeader
        ColHeaderSize(m,n) = length(colNames{m,n});
    end
end 

colSize = zeros(1,numCols);

for n = 1:numCols
    colSize(n) = max([ContentSize(:,n);ColHeaderSize(:,n)]);
end

% Align data:

for n = 1:numCols
    
    contentSizeMax = max(ContentSize(:,n));
    
    for m = 1:numRows
        
        contentPad = repmat(' ',1,contentSizeMax - ContentSize(m,n));
        colPad = repmat(' ',1,colSize(n) - contentSizeMax);
        colPadLength = round(length(colPad)/2);
        
        % Standardize lengths:
        
        switch lower(alignNumber{n})
            
            case {'decimal','left'}
                
                % Add white space to the right:
                Contents{m,n} = [colPad(1:colPadLength),Contents{m,n},contentPad,colPad(colPadLength+1:end)];
                
            case 'right'
                
                % Add white space to the left:
                Contents{m,n} = [colPad(1:colPadLength),contentPad,Contents{m,n},colPad(colPadLength+1:end)];
                
            case 'center'
                
                % Add white space to both the left and right:
                cut1 = round(length(contentPad)/2);
                Contents{m,n} = [colPad(1:colPadLength),contentPad(1:cut1),Contents{m,n},contentPad(cut1+1:end),colPad(colPadLength+1:end)];
                
            otherwise
                
                error('Align data to the ''left'',''center'', ''right'', or at the ''decimal''.')
                
        end
        
    end

end

% Align column headers:

for n = 1:numColsColHeader
    
    for m = 1:numRowsColHeader
        
        colHeaderPad = repmat(' ',1,colSize(n) - ColHeaderSize(m,n));
        
        % Standardize lengths:
        
        switch alignColName
            
            case 'left'
                
                % Add white space to the right:                
                colNames{m,n} = [colNames{m,n},colHeaderPad];
                
            case 'right'
                
                % Add white space to the left:
                colNames{m,n} = [colHeaderPad,colNames{m,n}];
                
            case 'center'
                
                % Add white space to both the left and right:
                colHeaderPadLength = round(length(colHeaderPad)/2);
                colNames{m,n} = [colHeaderPad(1:colHeaderPadLength),colNames{m,n},colHeaderPad(colHeaderPadLength+1:end)];
                
            otherwise
                
                error('Align column headers to the ''left'', ''right'', or ''center''.')
                
        end
        
    end
    
end

% Align row headers:

[numRowsRowHeader,numColsRowHeader] = size(rowNames);
RowHeaderSize = zeros(numRowsRowHeader,numColsRowHeader);

for n = 1:numColsRowHeader
    for m = 1:numRowsRowHeader
        RowHeaderSize(m,n) = length(rowNames{m,n});
    end    
end

rowHeaderSizeMax = max(RowHeaderSize,[],1);

if length(NWStr) > sum(rowHeaderSizeMax)
    rowHeaderSizeMax(:) = length(NWStr);
end

for n = 1:numColsRowHeader
    
    for m = 1:numRowsRowHeader
        
        rowHeaderPad = repmat(' ',1,rowHeaderSizeMax(n) - RowHeaderSize(m,n));
        
        % Standardize lengths:
                
        switch alignRowName
            
            case 'left'
                
                % Add white space to the right:                
                rowNames{m,n} = [rowNames{m,n},rowHeaderPad];
                
            case 'right'
                
                % Add white space to the left:
                rowNames{m,n} = [rowHeaderPad,rowNames{m,n}];
                
            case 'center'
                
                % Add white space to both the left and right:
                rowHeaderPadLength = round(length(rowHeaderPad)/2);
                rowNames{m,n} = [rowHeaderPad(1:rowHeaderPadLength),rowNames{m,n},rowHeaderPad(rowHeaderPadLength+1:end)];
                
            otherwise
                
                error('Align row headers to the ''left'', ''right'', or ''center''.')

        end
        
    end
    
end

% Pad each column:

if colPaddingLeft > 0 || colPaddingRight > 0
    
    columnPadLeft = repmat(' ',1,colPaddingLeft);
    columnPadRight = repmat(' ',1,colPaddingRight);
    
    for n = 1:numColsRowHeader
        
        for m = 1:numRowsRowHeader
            rowNames{m,n} = [columnPadLeft,rowNames{m,n},columnPadRight];
        end
        
    end
    
    for n = 1:numColsColHeader
        
        for m = 1:numRowsColHeader
            colNames{m,n} = [columnPadLeft,colNames{m,n},columnPadRight];
        end
        
    end    
    
    for m = 1:numRows
        
        for n = 1:numCols
            Contents{m,n} = [columnPadLeft,Contents{m,n},columnPadRight];
        end
        
    end
end

colSize = colSize+colPaddingLeft+colPaddingRight;

if printRowHeader
    rowHeaderSizeMax = rowHeaderSizeMax + colPaddingLeft + colPaddingRight;
end

% Compute the length of vertical inner border:

if ~iscell(borderInnerVert)
    
    borderInnerVerChar = borderInnerVert;
    borderInnerVert = cell(1,numCols-1);
    borderInnerVert(:) = {borderInnerVerChar};
    
end

nVerticalSpace = 0;

for n = 1:numCols-1
    nVerticalSpace = nVerticalSpace + length(borderInnerVert{n});
end

% Print the title:

if ~isempty(title)
    
    totalSpace = sum(rowHeaderSizeMax) + sum(colSize) + length(vertSep) + nVerticalSpace;
    titlePad = repmat(' ',1,max(0,totalSpace - length(title)));
    titlePadLength = round(length(titlePad)/2);
    fprintf(fileID,'<strong>%s%s%s</strong>\n',titlePad(1:titlePadLength),title,titlePad(titlePadLength+1:end));
    
end

% Print the top border:

if ~isempty(borderTop)
    
    borderTopBig = repmat(borderTop, 1, sum(rowHeaderSizeMax) + sum(colSize) + length(vertSep) + nVerticalSpace);
    
    if ~isempty(borderLeft)
        fprintf(fileID,'%s',' ');
    end
    
    fprintf(fileID,'%s',borderTopBig);
    
    if ~isempty(borderRight)
        fprintf(fileID,'%s',' ');
    end
    
    fprintf(fileID,'\n');
    
end

% Print the column headers and the northwest string:

if printColumnHeader
    
    for m = 1:numRowsColHeader
        
        fprintf(fileID,'%s',borderLeft);
        
        if m == numRowsColHeader
            
            NWPad = repmat(' ', 1, sum(rowHeaderSizeMax) - length(NWStr));
            fprintf(fileID,'%s%s%s',NWStr,NWPad,vertSep);
            
        else
            
            NWPad = repmat(' ', 1, sum(rowHeaderSizeMax));
            fprintf(fileID,'%s%s',NWPad,vertSep);
            
        end
        
        for n = 1:numCols-1
            fprintf(fileID,'%s%s',colNames{m,n},borderInnerVert{n});
        end
        
        fprintf(fileID,'%s%s\n',colNames{m,numCols},borderRight);
        
    end
end

% Print the horizontal separators:

if ~isempty(horzSep)
    
    if strcmp(horzSep,'_')
        
        fprintf(fileID,'%s',borderLeft);
        horizontalSep1 = repmat(horzSep,1,sum(rowHeaderSizeMax));
        fprintf(fileID,'%s%s',horizontalSep1,vertSep);
        
        for n = 1:numCols-1
            
            horizontalSepTemp = repmat(horzSep,1,colSize(n));
            fprintf(fileID,'%s%s',horizontalSepTemp,borderInnerVert{n});
            
        end
        
        horizontalSepTemp = repmat(horzSep,1,colSize(numCols));
        fprintf(fileID,'%s%s\n',horizontalSepTemp,borderRight);
        
    else
        
        horizontalSepBig = repmat(horzSep,1,sum(rowHeaderSizeMax) + sum(colSize) + length(vertSep) + nVerticalSpace);
        fprintf(fileID,'%s%s%s\n',borderLeft,horizontalSepBig,borderRight);
        
    end
    
end

% Print the contents:

for m = 1:numRows
    
    fprintf(fileID,'%s',borderLeft);
    
    if printRowHeader
        
        for n = 1:numColsRowHeader-1
            fprintf(fileID,'%s',rowNames{m,n});
        end
        
        fprintf(fileID,'%s%s',rowNames{m,numColsRowHeader},vertSep);
        
    end
        
    for n = 1:numCols-1
        fprintf(fileID,'%s%s',Contents{m,n},borderInnerVert{n});
    end
    
    fprintf(fileID,'%s%s\n',Contents{m,numCols},borderRight);
    
    % Inner horizontal border:
    
    if ~isempty(borderInnerHorz) && m < numRows
        
        if ischar(borderInnerHorz)
            borderInnerPrint = borderInnerHorz;
        elseif iscell(borderInnerHorz)
            borderInnerPrint = borderInnerHorz{m};
        else
            borderInnerPrint = ' ';
        end
        
        fprintf(fileID,'%s',borderLeft);
        innerSep1 = repmat(borderInnerPrint,1,sum(rowHeaderSizeMax));
        fprintf(fileID,'%s%s',innerSep1,vertSep);
        
        for n = 1:numCols-1
            
            innerSepTemp = repmat(borderInnerPrint,1,colSize(n));
            fprintf(fileID,'%s%s',innerSepTemp,borderInnerVert{n});
            
        end
        
        innerSepTemp = repmat(borderInnerPrint, 1, colSize(numCols));
        fprintf(fileID,'%s%s\n',innerSepTemp,borderRight);
        
    end
    
end
    
% Print the bottom border:

if ~isempty(borderBottom)
    
    if strcmp(borderBottom,'_')
        
        borderBottomBig = repmat(borderBottom, 1, sum(rowHeaderSizeMax) + sum(colSize) + length(vertSep) + nVerticalSpace);
        fprintf(fileID,'%s%s%s\n',borderLeft,borderBottomBig,borderRight);
        
    else
        
        borderBottomBig = repmat(borderBottom, 1, sum(rowHeaderSizeMax) + sum(colSize) + length(vertSep) + nVerticalSpace);
        
        if ~isempty(borderLeft)
            fprintf(fileID,'%s',' ');
        end
        
        fprintf(fileID,'%s',borderBottomBig);
        
        if ~isempty(borderRight)
            fprintf(fileID,'%s',' ');
        end
        
        fprintf(fileID,'\n');
        
    end
    
end

% Prepare tabular output:

if nargout > 0 || ~isempty(fileName)
    
    ind2D = 2*printColumnHeader+printRowHeader;
    
    switch ind2D
        
        case 0
            
            Table = Contents;
            
        case 1
            
            Table = [rowNames,Contents];
            
        case 2
            
            Table = [colNames;Contents];
            
        otherwise
            
            Table = cell(numRows + numRowsColHeader,numCols + numColsRowHeader);
            Table{1,1} = NWStr;
            Table(1:numRowsColHeader,numColsRowHeader+1:end) = colNames;
            Table(numRowsColHeader+1:end,1:numColsRowHeader) = rowNames;
            Table(numRowsColHeader+1:end,numColsRowHeader+1:end) = Contents;
            
    end
    
    Table = strtrim(Table);
    
    if ~isempty(fileName)
        xlswrite(fileName,Table);
    end
    
end

end

%------------------------------------------------------------------------
% Convert everything to character vector or cell array of character vectors
function B = convert2char(A)

if ischar(A)
    % Character vector -> character vector
    B = A;
elseif isstring(A)
    if isscalar(A)
        % String scalar -> character vector
        B = char(A);
    else
        % String vector/matrix -> cell array of character vectors
        % Temporary converter: B = cellstr(A)
        B = cell(size(A));
        for m = 1:numel(A)
            B{m} = char(A(m));
        end
    end
elseif iscell(A)
    % Cell array -> cell array of character vectors
    B = cell(size(A));
    for m = 1:numel(A)
        if isnumeric(A{m}) || islogical(A{m})
            B{m} = num2str(A{m});
        else
            try
                B{m} = char(A{m});
            catch
                error('Data elements in the cell array must be scalar.')
            end
            if size(B{m},1) > 1
                error('Data elements in the cell array must be scalar.')
            end
        end
    end
elseif isnumeric(A) || all(islogical(A(:)))
    if isempty(A) || isscalar(A)
        % Scalar number -> character vector
        B = num2str(A);
    else
        % Matrix -> cell array of character vectors
        B = cell(size(A));
        for m = 1:numel(A)
            B{m} = num2str(A(m));
        end
    end
end

end


 