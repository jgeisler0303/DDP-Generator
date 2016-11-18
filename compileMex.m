function compileMex(files, build_dir, build_file, compile_opt, errlog_file, is_matlab)

if ~is_matlab
    outstr= ' -o ';
    compile_cmd= 'mkoctfile --mex ';
else
    outstr= ' -output ';
    compile_cmd= 'mex ';    
    clear mex
end

try
    delete(errlog_file);
end

obj_dir= fullfile(build_dir, 'obj');
if ~exist(obj_dir, 'dir'), mkdir(obj_dir); end

objs= cell(1, length(files));    
fprintf('\n');
for i= 1:length(files)
    f= files{i};
    if f(1)=='-', continue, end
    fprintf('Compiling %s\n', f);
    [~, fn]= fileparts(f);
    objs{i}= fullfile(obj_dir, [fn '.o']);
    compile_str= [compile_opt ' -c ' f outstr objs{i} ' 2>>' errlog_file];
    
    [status, output]= system([compile_cmd compile_str]);
    
    if ~isempty(strtrim(output))
      disp(output);
    end
end
compile_str= [sprintf('%s ', objs{:}) outstr build_file ' 2>>' errlog_file];
[status, output]= system([compile_cmd compile_str]);

dd= dir(errlog_file);
if dd.bytes>0
    fprintf('There were compile errors logged to file: %s\n', errlog_file);
else
    fprintf('Successful\n');
end


