function mex_name= makeTarget(env, compile_switches, force_gen)

skip_gen= false;
if exist(env.maxima_errlog_file, 'file')
    dd= dir(env.maxima_errlog_file);
    last_err= dd.bytes>0;
end

if exist(fullfile(env.build_dir, 'iLQG_problem.h'), 'file') && (~exist('force_gen', 'var') || ~force_gen) && ~last_err
    dd_src= dir(env.problem_file);
    dd_tgt= dir(fullfile(env.build_dir, 'iLQG_problem.h'));
    if dd_src.datenum<dd_tgt.datenum
        skip_gen= true;
    end
end

if ~skip_gen
    if ispc
        command_line= ['-q --batch-string="batchload(\"make_iLQG.mac\")\$ gen_ddp(\"', strrep(env.problem_file, '\', '\\'), '\", \"', strrep(env.build_dir, '\', '\\'), '\")\$" 2>', env.maxima_errlog_file];
    else
        command_line= ['-q --batch-string="load(\"make_iLQG.mac\")\$ gen_ddp(\"', env.problem_file, '\", \"', env.build_dir, '\")\$" 2>', env.maxima_errlog_file];
    end

    maxima_res= callMaxima(command_line, env.maxima_errlog_file);
    
    if ~maxima_res
        mex_name= '';
        if exist(fullfile(env.target_dir, 'iLQG_problem.h'), 'file')
            movefile(fullfile(env.target_dir, 'iLQG_problem.h'), fullfile(env.target_dir, 'faulty_iLQG_problem.h')); 
        end
        if exist(fullfile(env.target_dir, 'iLQG_func.cpp'), 'file')
            movefile(fullfile(env.target_dir, 'iLQG_func.cpp'), fullfile(env.target_dir, 'faulty_iLQG_func.cpp'));
        end
            
        fprintf("Maxima quit with an error, stopping further processing\n");
        return
    end
    
end

if ~exist('compile_switches', 'var')
    compile_switches= '-DDEBUG_BACKPASS=1 -DDEBUG_FORWARDPASS=1';
    compile_switches= '';
end

mex_name= ['ddp' env.problem_title];
build_file= fullfile(fileparts(env.problem_file), [mex_name '.' mexext]);

compile_opt= ['-DPRNT=mexPrintf ', compile_switches, ' -v -I../common -I', env.build_dir, filesep];

compile_files{1}= '../common/iLQG.c';
compile_files{end+1}= '../common/iLQG_mex.c';
compile_files{end+1}= '../common/line_search.c';
compile_files{end+1}= '../common/back_pass.c';
compile_files{end+1}= '../common/matMult.c';
compile_files{end+1}= '../common/boxQP.c';
compile_files{end+1}= '../common/cholesky.c';
compile_files{end+1}= '../common/printMat.c';
compile_files{end+1}= fullfile(env.build_dir, 'iLQG_func.c');

old_cflags= getenv('CFLAGS');
old_xtracflags= getenv('XTRA_CFLAGS');
[status, cflags]= system('mkoctfile --print  CFLAGS');
cflags= strrep(cflags, '-fopenmp', '');
cflags= strrep(cflags, '-pthread', '');
cflags= [cflags ' -pthread -fopenmp -std=gnu99 '];
%cflags= [cflags ' -std=gnu99 '];
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
