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

mex_name= ['ddp' env.problem_title '_eigen3'];
build_file= fullfile(fileparts(env.problem_file), [mex_name '.' mexext]);


compile_opt= [compile_switches, ' -v -DEIGEN_NO_MALLOC -DEIGEN_NO_DEBUG -DEIGEN_DONT_PARALLELIZE -I', env.build_dir, filesep];

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
cflags= strrep(cflags, '-fopenmp', '');
cflags= strrep(cflags, '-pthread', '');
cflags= [cflags ' -pthread -fopenmp '];
% enable SSE or NEON
% EIGEN_VECTORIZE is set when vectorization is on
% check for supported SIMD: http://stackoverflow.com/questions/28939652/how-to-detect-sse-avx-avx2-availability-at-compile-time
% gcc -march=native -dM -E - < /dev/null | egrep "SSE|AVX" | sort
% see http://eigen.tuxfamily.org/index.php?title=FAQ#How_can_I_enable_vectorization.3F
% maybe -O3 -ftree-vectorize
% EIGEN_NO_DEBUG // no range checking!

cflags= strrep(cflags, char(10), ' ');
setenv('CXXFLAGS', cflags);

[status, xtracflags]= system('mkoctfile --print  XTRA_CFLAGS');
xtracflags= strrep(xtracflags, '-fopenmp', '');
xtracflags= strrep(xtracflags, '-pthread', '');
xtracflags= strrep(xtracflags, char(10), ' ');
setenv('XTRA_CXXFLAGS', xtracflags);



compileMex(compile_files, env.build_dir, build_file, compile_opt, env.compile_errlog_file, env.is_matlab)

setenv('CXXFLAGS', old_cflags);
setenv('XTRA_CXXFLAGS', old_xtracflags);
