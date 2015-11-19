function compileDDPProblem(fun, checkNaN)

if exist('checkNaN')~=1
    checkNaN= 0;
end

o= feval(fun);
f= o.f;
F= o.F;
L= o.L;
x_= o.x_(:)';
u_= o.u_(:)';
vars= [x_(:)' u_(:)'];
params= o.params;
param= [];
param_names= sortParams;

fn= 'ddpProblemDescription.h';
fid= fopen(fn, 'w');
fprintf(fid, '#ifdef DDP_PROBLEM\n');
if ~isempty(param_names)
    fprintf(fid, '    int n_params= %d;\n', length(param_names));
    for i= 1:length(param_names)
        fprintf(fid, '    tParamDesc p_name%d= {"%s", %d};\n', i, param_names{i}, getfield(param, param_names{i}));
    end    
    fprintf(fid, '    tParamDesc *paramdesc[]= { ');
    for i= 1:length(param_names)
        fprintf(fid, '&p_name%d', i);
        if i~=length(param_names), fprintf(fid, ', '); else fprintf(fid, ' };\n\n'); end
    end
else
    fprintf(fid, '    int n_params= %d;\n', 0);    
end
fprintf(fid, '    #define N_X %d\n', length(x_));
fprintf(fid, '    #define N_U %d\n', length(u_));
% writeHAS(L, x_, [], 'LX');
% writeHAS(L, u_, [], 'LU');
% writeHAS(L, x_, x_, 'LXX');
% writeHAS(L, u_, u_, 'LUU');
% writeHAS(L, x_, u_, 'LXU');
writeHAS(f, x_, x_, 'FXX');
writeHAS(f, u_, u_, 'FUU');
writeHAS(f, x_, u_, 'FXU');
if isTV
   fprintf(fid, '    #define IS_TV\n');
else
   fprintf(fid, '    #undef IS_TV\n');
end

fprintf(fid, '\n    double ddpJ(double x[], double u[], int k, double *p[], double N) {\n');
fprintf(fid, '        double t0;\n');
fprintf(fid, '        if(k>=0) {\n');
fprintf(fid, '            %s\n', replaceVars(L));
fprintf(fid, '        } else {\n');
fprintf(fid, '            %s\n', replaceVars(F));
fprintf(fid, '        }\n');
fprintf(fid, '        return t0;\n');
fprintf(fid, '    }\n');
    
fprintf(fid, '\n    int ddpf(double x_next[], double x[], double u[], int k, double *p[], double N) {\n');
for i= 1:length(x_)
    nx= ['x_next[' num2str(i-1) ']'];
    s= strrep(replaceVars(f(i)), 't0', nx);
    fprintf(fid, '        %s\n', s);
    writeCheckNaN(nx, '0', s);
end
fprintf(fid, '        return 1;\n');
fprintf(fid, '    }\n');
fprintf(fid, '#endif\n\n');

fprintf(fid, '#ifdef DDP_FXK\n');
writeFirst(fid, f, x_, 'fx', 1);
fprintf(fid, '\n');
writeFirst(fid, f, u_, 'fu', 1);
fprintf(fid, '\n');
writeSecVec(fid, f, x_, x_, 'fxx', 1);
fprintf(fid, '\n');
writeSecVec(fid, f, u_, u_, 'fuu', 1);
fprintf(fid, '\n');
writeSecVec(fid, f, x_, u_, 'fxu');
fprintf(fid, '#endif\n\n');

fprintf(fid, '#ifdef DDP_F\n');
writeFirst(fid, F, x_, 'Vx', 1);
fprintf(fid, '\n');
writeSecond(fid, F, x_, x_, 'Vxx', 1, 1);
fprintf(fid, '#endif\n\n');

fprintf(fid, '#ifdef DDP_L\n');
writeFirst(fid, L, x_, 'Qx', 1);
fprintf(fid, '\n');
writeSecond(fid, L, x_, x_, 'Qxx', 1, 1);
fprintf(fid, '\n');    
writeFirst(fid, L, u_, 'Qu', 1);
fprintf(fid, '\n');  
writeSecond(fid, L, u_, u_, 'Quu', 1, 1);
fprintf(fid, '\n');    
writeSecond(fid, L, x_, u_, 'Qxu', 0, 1);
fprintf(fid, '#endif\n');

fclose(fid);

fn= char(fun);
if length(fn)>6 & strcmp(fn(1:6), 'optDef')
    fn= fn(7:end);
end
mex('-output', ['ddp' fn], 'ddpMex.c')

    function s= replaceVars(f)
%         s= ccode(simplify(f));
        s= ccode(f);
        s= s(7:end);
        s= strrep(s, '~', '');
        for j= 1:length(x_)
            s= strrep(s, char(x_(j)), ['x[' num2str(j-1) ']']);
        end
        for j= 1:length(u_)
            s= strrep(s, char(u_(j)), ['u[' num2str(j-1) ']']);
        end
        for j= 1:length(params)
            pn= char(params(j));
            pos= strfind(pn, '_');
            if ~isempty(pos)
                idx= pn(pos+1:end);
                if isempty(idx)
                    idx= '0';
                else
                    if idx~='k'
                        idx= str2num(idx);
                        if isempty(idx)
                            idx= '0';
                        else
                            idx= num2str(idx-1);
                        end
                    end
                end
                pn_= pn(1: pos-1);
            else
                idx= '0';
                pn_= pn;
            end
            pi= find(strcmp(param_names, pn_));
            if length(pi)~=1
                fclose(fid);
                error('parameter nicht in liste gefunden');
            end
            pn_= ['(p[' num2str(pi-1) '])'];
            s= strrep(s, pn, ['(' pn_ '[' idx '])']);
        end            
    end

    function j= isTV
        j= 0;
        rem= findsym(f);
        while ~isempty(rem)
            [v, rem]= strtok(rem, ',');
            pos= strfind(v, '_');
            if length(pos)>1
                v= v(1:pos(1)-1);
            end
            if (isfield(param, v) & getfield(param, v)==-1) | v=='k' | ismember(sym(v), vars)
                j= 1;
                return;
            end
        end
    end

    function writeHAS(f, x1, x2, s)
        if isempty(x2)
            has= all(all(jacobian(f, x1)==0));
        else
            has= all(all(jacobian(jacobian(f, x1), x2)==0));
        end            
        if has
            fprintf(fid, '    #undef HAS_%s\n', s);
        else
            fprintf(fid, '    #define HAS_%s\n', s);
        end
    end

    function writeFirst(fid, f, dx, nm, writezeros)
        df= jacobian(f, dx);
        if all(all(df==0))
            if (exist('writezeros')~=1 | ~writezeros)
                return
            else
                fprintf(fid, '    for(i= 0; i<%d; i++) %s[i]= 0.0;\n', prod(size(df)), nm);
                return
            end
        end
        for j= 1:prod(size(df))
            nm_= [nm '[' num2str(j-1) ']'];
            s= strrep(replaceVars(df(j)), 't0', nm_);
            fprintf(fid, '    %s\n', s);
            writeCheckNaN(nm_, 'mxGetNaN()', s);
        end
    end

    function writeSecond(fid, f, dx1, dx2, nm, tri, writezeros)
        if exist('tri')~=1, tri= 0; end
        ddf= jacobian(jacobian(f, dx1), dx2);
        [n, m]= size(ddf);
        if tri & n~=m, error('nur symmetrische matrizen können in dreiecksform gespeichert werden'); end
        if all(all(ddf==0))
            if (exist('writezeros')~=1 | ~writezeros)
                return
            else
                if tri
                    nn= (n*(n+1))/2;
                else
                    nn= n*m;
                end
                fprintf(fid, '    for(i= 0; i<%d; i++) %s[i]= 0.0;\n', nn, nm);
                return
            end
        end
        l= 0;
        for k= 1:m
            for j= 1:n
                if ~tri | k>=j
                    nm_= [nm '[' num2str(l) ']'];
                    s= strrep(replaceVars(ddf(j, k)), 't0', nm_);
                    fprintf(fid, '    %s\n', s);
                    writeCheckNaN(nm_, 'mxGetNaN()', s);
                    l= l+1;
                end
            end
            if k~=m, fprintf(fid, '\n'); end
        end
    end

    function writeSecVec(fid, f, dx1, dx2, nm, tri, writezeros)
        if exist('tri')~=1, tri= 0; end
        n= length(dx1);
        m= length(dx2);
        if tri & n~=m, error('nur symmetrische matrizen können in dreiecksform gespeichert werden'); end
        if all(all(jacobian(jacobian(f, dx1), dx2)==0))
            if (exist('writezeros')~=1 | ~writezeros)
                return
            else
                if tri
                    nn= length(f)*(n*(n+1))/2;
                else
                    nn= length(f)*n*m;
                end
                fprintf(fid, '    for(i= 0; i<%d; i++) %s[i]= 0.0;\n', nn, nm);
                return
            end
        end
        l= 0;
        for kf= 1:length(f);
            ddf= jacobian(jacobian(f(kf), dx1), dx2);
            for k= 1:m
                for j= 1:n
                    if ~tri | k>=j
                        nm_= [nm '[' num2str(l) ']'];
                        s= strrep(replaceVars(ddf(j, k)), 't0', nm_);
                        fprintf(fid, '    %s\n', s);
                        writeCheckNaN(nm_, 'mxGetNaN()', s);
                        l= l+1;
                    end
                end
                if k~=m, fprintf(fid, '\n'); end
            end
            fprintf(fid, '//-\n');
        end
    end

    function pn= sortParams
        for j= 1:length(params)
            s= char(params(j));
            pos= strfind(s, '_');
            if length(pos)>1
                error('parameter dürfen nur einen unterstrich im namen tragen, der den namen vom index trennt');
            end
            if isempty(pos)
                nm= s;
                l= 1;
            else
                nm= s(1:pos-1);
                s= s(pos+1:end);
                if length(s)==1 & s=='k'
                    l= -1;
                else
                    l= str2num(s);
                    if isempty(l)
                        l= 1;
                    end
                end
            end
            if ~isfield(param, nm)
                param= setfield(param, nm, l);
            else
                param= setfield(param, nm, max(l, getfield(param, nm)));
            end
        end
        if ~isempty(param), pn= fieldnames(param); else pn= {}; end
    end

    function writeCheckNaN(vn, retval, line)
        if ~checkNaN, return; end
        fprintf(fid, '        if(mxIsNaN(%s) || mxIsInf(%s)) {\n', vn, vn);
        fprintf(fid, '            mexPrintf("    @k %%d: %s is nan or inf\\n", k);\n', line);
%         fprintf(fid, '            mexPrintf("    p[3]= %%g\\n", (p[3])[0]);\n');
%         fprintf(fid, '            mexPrintf("    p[1]= %%g\\n", (p[1])[0]);\n');
%         fprintf(fid, '            mexPrintf("    p[5]= %%g\\n", (p[5])[0]);\n');
%         fprintf(fid, '            mexPrintf("    N= %%g\\n", N);\n');
%         fprintf(fid, '            mexPrintf("    %s= %%g\\n", %s);\n', vn, vn);
%         fprintf(fid, '            mexPrintf("    %s= %%g\\n", %s);\n', 'x[5]', 'x[5]');
%         fprintf(fid, '            mexPrintf("    %s= %%g\\n", %s);\n', 'u[0]', 'u[0]');
%         fprintf(fid, '            mexPrintf("    %s= %%g\\n", %s);\n', '((p[3])[0])*(sqrt(u[0]*u[0]+((p[1])[0]))-((p[5])[0]))/N', '((p[3])[0])*(sqrt(u[0]*u[0]+((p[1])[0]))-((p[5])[0]))/N');
%         fprintf(fid, '            mexPrintf("    %s= %%g\\n", %s);\n', 'sqrt(u[0]*u[0]+((p[1])[0]))', 'sqrt(u[0]*u[0]+((p[1])[0]))');
%         fprintf(fid, '            mexPrintf("    %s= %%g\\n", %s);\n', 'u[0]*u[0]+((p[1])[0])', 'u[0]*u[0]+((p[1])[0])');
        fprintf(fid, '            return %s;\n', retval);
        fprintf(fid, '        }\n');
    end
end