function make_iLQG(problem, compile_switches, force_gen)

if ~exist(problem, 'file')
    problem= [problem '.mac'];
    if ~exist(problem, 'file')
        error(['could not find problem ' problem]);
    end
end

[problem_path, problem_name, problem_ext]= fileparts(problem);

if length(problem_name)>6 && strcmp(problem_name(1:6), 'optDef')
  problem_title= problem_name(7:end);
else
  problem_title= problem_name;
end

v= ver;
is_matlab= strcmp(v(1).Name, 'Octave');

base_dir= fileparts(mfilename ('fullpath'));
problem_file= which(fullfile(problem_path, problem));
target_dir= fullfile(fileparts(problem_file), [problem_title, '_gen_files']);
target_file= fullfile(fileparts(problem_file), ['iLQG' problem_title '.' mexext]);
maxima_errlog_file= fullfile(target_dir, 'maxima_errors.log');
mex_errlog_file= fullfile(target_dir, 'err.log');

if ~exist(target_dir, 'dir'), mkdir(target_dir); end

skip_gen= false;
if exist(fullfile(target_dir, 'iLQG_problem.h'), 'file') && (~exist('force_gen', 'var') || ~force_gen)
    dd_src= dir(problem_file);
    dd_tgt= dir(fullfile(target_dir, 'iLQG_problem.h'));
    if dd_src.datenum<dd_tgt.datenum
        skip_gen= true;
    end
end

old_dir= pwd;
cd(base_dir);

if ~skip_gen
    maxima= getenv('MAXIMA');
    if isempty(maxima)
        maxima= 'maxima';
    end
    if ispc
        command_str= [maxima ' -q --batch-string="problem_file:\"', strrep(problem_file, '\', '\\'), '\";target_dir:\"', strrep(target_dir, '\', '\\'), '\";batchload(\"make_iLQG.mac\");" 2>' maxima_errlog_file];
    else
        command_str= [maxima ' -q --batch-string="problem_file:\"', problem_file, '\"; target_dir:\"', target_dir, '\"; batchload(\"make_iLQG.mac\");" 2>' maxima_errlog_file];
    end
    [status, output]= system(command_str);
    
    fprintf('Maxima said:\n%s\n', output);
    dd= dir(maxima_errlog_file);
    if dd.bytes>0
        fprintf('There were errors during execution of maxima logged to file: %s\n', maxima_errlog_file);
    end
end

if ~exist('compile_switches', 'var')
    compile_switches= '-DDEBUG_BACKPASS=1 -DDEBUG_FORWARDPASS=1';
end

if is_matlab
    compile_cmd= 'mkoctfile --mex ';
    compile_opt= ['-DPRNT=mexPrintf ', compile_switches, ' -I. -I', target_dir, filesep, ' -o ', target_file];
else
    compile_cmd= 'mex ';    
    compile_opt= ['-DPRNT=mexPrintf ', compile_switches, ' -I. -I', target_dir, filesep, ' -output ', target_file];
    clear mex
end

compile_files{1}= 'iLQG.c';
compile_files{end+1}= 'iLQG_mex.c';
compile_files{end+1}= 'line_search.c';
compile_files{end+1}= 'back_pass.c';
compile_files{end+1}= 'matMult.c';
compile_files{end+1}= 'boxQP.c';
compile_files{end+1}= 'cholesky.c';
compile_files{end+1}= 'printMat.c';
compile_files{end+1}= fullfile(target_dir, 'iLQG_func.c');

compile_str= [compile_opt ' ' sprintf('%s ', compile_files{:}) ' 2>' mex_errlog_file];
fprintf('\nCompiling...\n');
[status, output]= system([compile_cmd compile_str]);
if status>0
    disp(output)
else
    dd= dir(mex_errlog_file);
    if dd.bytes>0
        fprintf('There were compile errors logged to file: %s\n', mex_errlog_file);
    else
        fprintf('Successful\n');
    end
end

cd(old_dir);
