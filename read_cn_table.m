function CN_table = read_cn_table(filename, firstLine)
%READ_CN_TABLE Read a single normal-force coefficient grid from a text file.
%   CN_table = read_cn_table('CN.txt') begins searching on line 46.
%   CN_table = read_cn_table('CN.txt', firstLine) changes the starting line.
%
%   A header such as "#ang/Mach  0.5  1.0  2.0" supplies the Mach values.
%   A separate "#" may precede it. Each following nonblank row contains
%   angle (deg) and one CN value per Mach. The next line starting with #
%   ends the table. Values may be separated by spaces, tabs, or commas.
%
%   CN_table.header       original header text
%   CN_table.Mach         increasing Mach breakpoints, 1-by-nMach
%   CN_table.data(:,1)    increasing angle breakpoints in degrees
%   CN_table.data(:,2:end)  CN, nAngles-by-nMach
if nargin < 2
    firstLine = 46;
end
validateattributes(firstLine, {'numeric'}, {'scalar','real','finite','integer','positive'});
lines = regexp(fileread(filename), '\r\n|\n|\r', 'split');

% Find the named header, so the opening "#" cannot be mistaken for the end.
pattern = '^\s*#?\s*(?:ang(?:le)?|alpha|aoa)\s*/\s*mach(?=[\s,;:]|$)[\s,;:]*';
headerLine = [];
for k = firstLine:numel(lines)
    labelEnd = regexpi(lines{k}, pattern, 'end', 'once');
    if ~isempty(labelEnd)
        headerLine = k;
        break
    end
end
if isempty(headerLine)
    error('read_cn_table:Header', ...
        'No angle/Mach header found in %s at or after line %d.', filename, firstLine);
end
header = strtrim(lines{headerLine});
mach = numericRow(lines{headerLine}(labelEnd+1:end), filename, headerLine, false);
if numel(mach) < 2 || any(mach < 0)
    error('read_cn_table:Mach', ...
        '%s line %d must list at least two nonnegative Mach values.', filename, headerLine);
end

endLine = headerLine + 1;
while endLine <= numel(lines) && isempty(regexp(lines{endLine}, '^\s*#', 'once'))
    endLine = endLine + 1;
end
if endLine > numel(lines)
    error('read_cn_table:EndMarker', ...
        'Missing closing # after the header on line %d of %s.', headerLine, filename);
end

width = numel(mach) + 1;
data = zeros(endLine - headerLine - 1, width);
nRows = 0;
for k = headerLine+1:endLine-1
    row = strtrim(lines{k});
    if isempty(row)
        continue
    end
    values = numericRow(row, filename, k, true);
    if numel(values) ~= width
        error('read_cn_table:RowWidth', ...
            '%s line %d has %d columns. Expected one angle and %d CN values.', ...
            filename, k, numel(values), numel(mach));
    end
    nRows = nRows + 1;
    data(nRows,:) = values;
end
if nRows < 2
    error('read_cn_table:Angles', '%s needs at least two angle rows.', filename);
end
data = data(1:nRows,:);

% Reorder the matching coefficients whenever a breakpoint axis is sorted.
[mach, columns] = sort(mach);
[angles, rows] = sort(data(:,1));
if any(diff(mach) <= 0) || any(diff(angles) <= 0)
    error('read_cn_table:Breakpoints', 'Duplicate angle or Mach breakpoint in %s.', filename);
end
CN_table = struct('header', header, 'Mach', mach, ...
                  'data', [angles, data(rows, 1 + columns)]);
end

function values = numericRow(line, filename, lineNumber, hasAngle)
tokens = regexp(strtrim(line), '[\s,;]+', 'split');
if hasAngle
    % Accept an optional degree symbol immediately after the angle.
    tokens{1} = regexprep(tokens{1}, [char(176) '$'], '');
end
% Some engineering tables write scientific notation with D in place of E.
values = str2double(regexprep(tokens, '[dD]', 'E'));
if isempty(values) || ~isreal(values) || any(~isfinite(values))
    error('read_cn_table:Number', ...
        'Expected finite real numeric entries in %s on line %d.', filename, lineNumber);
end
end
