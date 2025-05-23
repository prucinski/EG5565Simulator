%This script has been fully generated by ChatGPT :)
function builtString = extractDate(filename)    
    % Define the regular expression pattern to match the date and time part of the filename
    pattern = '\d{8}_\d{6}';
    % Use regular expression to extract the date and time from the filename
    dateTimeStr = regexp(filename, pattern, 'match');
    % Extract date and time components
    dateStr = dateTimeStr{1}(1:8); % Extract date
    timeStr = dateTimeStr{1}(10:end); % Extract time
    % Convert the date string into a MATLAB datetime object
    dateObj = datetime(dateStr, 'InputFormat', 'yyyyMMdd');
    % Extract hour, minute, and second components from the time string
    hour = str2double(timeStr(1:2));
    minute = str2double(timeStr(3:4));
    second = str2double(timeStr(5:6));
    % Format the time string as hh:mm:ss
    timeString = sprintf('%02d:%02d:%02d', hour, minute, second);
    % Convert the datetime object into a string format
    dateString = datestr(dateObj, 'yyyy-mm-dd');
    % Display the extracted date and time as strings
    builtString = strcat(dateString, ", ", timeString, newline);
end