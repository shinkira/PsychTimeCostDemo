cur_dir = pwd;
cd(fullfile(get_root,'DynamicProgramming'))

v1 = coder.typeof(1);
v2 = coder.typeof(1);
v3 = coder.typeof(1);
v4 = coder.typeof(1);
v5 = coder.typeof(0,[1,inf]);
v6 = coder.typeof(1);
v7 = coder.typeof(1);
v8 = coder.typeof(0,[1,4],[0 0]);
v9 = coder.typeof(1);
v10 = coder.typeof(0,[1,inf]);
v11 = coder.typeof(0,[1,2]);
v12 = coder.typeof(0,[inf,inf]);
v13 = coder.typeof(0,[inf,inf]);

codegen -config:mex d2b_dp.m -args {v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, v13} -report
cd(cur_dir);