function mex_name= makeTarget(env, compile_switches, force_gen)

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

mex_name= ['ddp' env.problem_title '_eigen3'];
build_file= fullfile(fileparts(env.problem_file), [mex_name '.' mexext]);

compile_opt= [compile_switches, ' -v -I', env.build_dir, filesep];

compile_files{1}= 'iterate.cpp';
compile_files{end+1}= 'iLQG_mex.cpp';
compile_files{end+1}= 'line_search.cpp';
compile_files{end+1}= 'back_pass.cpp';
compile_files{end+1}= 'boxQP.cpp';
compile_files{end+1}= 'printMat.cpp';
compile_files{end+1}= 'iLQG.cpp';
compile_files{end+1}= fullfile(env.build_dir, 'iLQG_func.cpp');

old_cflags= getenv('CXXFLAGS');
old_xtracflags= getenv('XTRA_CXXFLAGS');
[status, cflags]= system('mkoctfile --print  CXXFLAGS');
%cflags= strrep(cflags, '-fopenmp', '');
%cflags= strrep(cflags, '-pthread', '');
% enable SSE or NEON
% see http://eigen.tuxfamily.org/index.php?title=FAQ#How_can_I_enable_vectorization.3F
% use g++!
% -DEIGEN_DONT_PARALLELIZE
% maybe -O3 -ftree-vectorize
% EIGEN_NO_DEBUG // no range checking!

%cflags= [cflags ' '];
%cflags= [cflags ' -std=gnu99 '];
%cflags= strrep(cflags, char(10), ' ');
%setenv('CFLAGS', cflags);
%[status, xtracflags]= system('mkoctfile --print  XTRA_CFLAGS');
%xtracflags= strrep(xtracflags, '-fopenmp', '');
%xtracflags= strrep(xtracflags, '-pthread', '');
%xtracflags= strrep(xtracflags, char(10), ' ');
%setenv('XTRA_CFLAGS', xtracflags);



compileMex(compile_files, env.build_dir, build_file, compile_opt, env.compile_errlog_file, env.is_matlab)

setenv('CXXFLAGS', old_cflags);
setenv('XTRA_CXXFLAGS', old_xtracflags);