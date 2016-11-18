function makeTarget(env, compile_switches, force_gen)

skip_gen= false;
if exist(fullfile(env.build_dir, 'iLQG_problem.h'), 'file') && (~exist('force_gen', 'var') || ~force_gen)
    dd_src= dir(env.problem_file);
    dd_tgt= dir(fullfile(env.build_dir, 'iLQG_problem.h'));
    if dd_src.datenum<dd_tgt.datenum
        skip_gen= true;
    end
end

if ~skip_gen
    if ispc
        command_line= ['-q --batch-string="problem_file:\"', strrep(env.problem_file, '\', '\\'), '\";target_dir:\"', strrep(env.build_dir, '\', '\\'), '\";batchload(\"make_iLQG.mac\");" 2>' env.maxima_errlog_file];
    else
        command_line= ['-q --batch-string="problem_file:\"', env.problem_file, '\"; target_dir:\"', env.build_dir, '\"; batchload(\"make_iLQG.mac\");" 2>' env.maxima_errlog_file];
    end

    callMaxima(command_line, env.maxima_errlog_file);
end

if ~exist('compile_switches', 'var')
    compile_switches= '-DDEBUG_BACKPASS=1 -DDEBUG_FORWARDPASS=1';
end

build_file= fullfile(fileparts(env.problem_file), ['ddp' env.problem_title '_mp.' mexext]);

compile_opt= ['-DPRNT=mexPrintf ', compile_switches, ' -v -I. -I../common -I', env.build_dir, filesep];

compile_files{1}= 'iLQG.c';
compile_files{end+1}= 'iLQG_mex.c';
compile_files{end+1}= 'line_search.c';
compile_files{end+1}= 'back_pass.c';
compile_files{end+1}= '../common/matMult.c';
compile_files{end+1}= '../common/boxQP.c';
compile_files{end+1}= '../common/cholesky.c';
compile_files{end+1}= '../common/printMat.c';
compile_files{end+1}= fullfile(env.build_dir, 'iLQG_func.c');
%compile_files{end+1}= '-lgomp';

old_cflags= getenv('CFLAGS');
old_xtracflags= getenv('XTRA_CFLAGS');
[status, cflags]= system('mkoctfile --print  CFLAGS');
cflags= strrep(cflags, '-fopenmp', '');
cflags= strrep(cflags, '-pthread', '');
cflags= [cflags ' -pthread -fopenmp -std=gnu99 '];
%cflags= [cflags '  -std=c99 '];
cflags= strrep(cflags, char(10), ' ');
setenv('CFLAGS', cflags);
[status, xtracflags]= system('mkoctfile --print  XTRA_CFLAGS');
xtracflags= strrep(xtracflags, '-fopenmp', '');
xtracflags= strrep(xtracflags, '-pthread', '');
xtracflags= strrep(xtracflags, char(10), ' ');
setenv('XTRA_CFLAGS', xtracflags);

compileMex(compile_files, env.build_dir, build_file, compile_opt, env.compile_errlog_file, env.is_matlab)

setenv('CFLAGS', old_cflags);
setenv('CFLAGS', old_xtracflags);