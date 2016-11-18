function make_target(problem, target, varargin)

base_dir= fileparts(mfilename ('fullpath'));
target_dir= fullfile(base_dir, 'target', target);

if ~exist(target_dir, 'dir')
    error('Target %s not found.', target);
end
if ~exist(fullfile(target_dir, 'makeTarget.m'), 'file')
    error('No function "makeTarget" found in target %s.', target);
end


if ~exist(problem, 'file')
    env.problem= [problem '.mac'];
    if ~exist(env.problem, 'file')
        error(['could not find problem ' env.problem]);
    end
else
    env.problem= problem;
end

[env.problem_path, env.problem_name, problem_ext]= fileparts(problem);

if length(env.problem_name)>6 && strcmp(env.problem_name(1:6), 'optDef')
  env.problem_title= env.problem_name(7:end);
else
  env.problem_title= env.problem_name;
end

v= ver;
env.is_matlab= ~strcmp(v(1).Name, 'Octave');

env.problem_file= which(fullfile(env.problem_path, env.problem));
env.build_dir= fullfile(fileparts(env.problem_file), [env.problem_title, '_' strrep(target, '/', '_')]);
env.maxima_errlog_file= fullfile(env.build_dir, 'maxima_errors.log');
env.compile_errlog_file= fullfile(env.build_dir, 'compile_err.log');

if ~exist(env.build_dir, 'dir'), mkdir(env.build_dir); end

old_dir= pwd;
cd(target_dir);

makeTarget(env, varargin{:});

cd(old_dir);
