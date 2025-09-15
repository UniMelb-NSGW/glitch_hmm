addpath(matlab_path)
addpath([matlab_path '/HMM'])
addpath([matlab_path '/Utilities'])

%addpath('/fred/oz022/ldunn_glitch/matlab_code_logs/')
%addpath('/fred/oz022/ldunn_glitch/matlab_code_logs/HMM')
%addpath('/fred/oz022/ldunn_glitch/matlab_code_logs/Utilities')


%global f_fiducial fd_fiducial fdd_fiducial sigma kappa_per_toa

freqs = load([working_prefix 'freqs.dat']);
fdots = load([working_prefix 'fdots.dat']);
zs = load([working_prefix 'zs.dat']);
kappa_per_toa = load([working_prefix 'kappas.dat']);

%f_fiducial = load('test_freq_fiducial.dat');
%fdd_fiducial = load('test_fddot_fiducial.dat');
%fd_fiducial = -load('test_fdot_fiducial.dat');

%sigma = max(1e-21, mean(diff(test_fdots))/mean(test_zs)^0.5)
%kappa_per_toa = test_kappas;

save([working_prefix 'init.mat']);

[J,max_BFs,E,E0] = seq_glitch_det(zs, freqs', fdots');
[path,full_posterior,alpha,beta,evidence,s, residuals] = HMM_Pulse_3dcol(zs, freqs', fdots',[J]);
f = fopen([out_prefix 'res.dat'], 'w');
for i = 1:max(size(J))
    fprintf(f, "%d %.16f\n", J(i), max_BFs(i));
end
fclose(f);

f = fopen([out_prefix 'fdot_path.dat'], 'w');
fprintf(f,"%d\n", fdots(path(:,1)));
fclose(f);

f = fopen([out_prefix 'f_path.dat'], 'w');
fprintf(f,"%d\n", freqs(path(:,2)));
fclose(f);

for n = 1:size(path,1) freq_posterior(:, n) = logsumexp(full_posterior(:, :, n),1); end;
writematrix(freq_posterior, [out_prefix 'f_posterior.dat'], 'Delimiter', ' ');

for n = 1:size(path,1) fdot_posterior(:, n) = logsumexp(full_posterior(:, :, n),2); end;
writematrix(fdot_posterior, [out_prefix 'fdot_posterior.dat'], 'Delimiter', ' ');

all_BFs = E;
for i = 1:(numel(J)+1)
    all_BFs(i,:) = all_BFs(i,:) - E0(i);
end

writematrix(all_BFs, [out_prefix 'bfs.dat'], 'Delimiter', ' ');

save([out_prefix 'analysis.mat'], '-v7')
