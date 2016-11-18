function callMaxima(command_line, maxima_errlog_file)
  
maxima= getenv('MAXIMA');
if isempty(maxima)
    maxima= 'maxima';
end

command_str= [maxima ' ' command_line];

[status, output]= system(command_str);

fprintf('Maxima said:\n%s\n', output);
dd= dir(maxima_errlog_file);
if dd.bytes>0
    fprintf('There were errors during execution of maxima logged to file: %s\n', maxima_errlog_file);
end
