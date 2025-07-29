p = mfilename('fullpath');
[p,~,~] = fileparts(p);
cd(p);
% addpath('/Volumes/Data/1 - LAB/01 - Proyectos/66 - UrgencyShin/models/generic');

if 1
    tl = [5,  0.1 ,0.001, 0, 0];
    th = [40, 0.4 ,0.08, 0, 0];
    tg = [15, 0.2 ,0.02, 0, 0];
else
    tl = [5,  0.1 ,0, 0, 0];
    th = [40, 0.4 ,0, 0, 0];
    tg = [15, 0.2 ,0, 0, 0];
end
fn_wrapper = @wrapper_dtb_empiricalbound_rt;
fit_diffusion(fn_wrapper,tl,th,tg)