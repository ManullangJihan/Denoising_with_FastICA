%% Program Utama untuk FastICA Single Input
% Nama : Kharisma Putri Nabila
% NIM : 16101174

%% Environment
clear all; close all; close all hidden; clc

%% Menginput data
[fname,pname] = uigetfile('*.wav','Pilih sebuah Folder PCG');
pcgfile = fullfile(pname,fname);
[x, fs] = audioread(pcgfile);

% Mengatur Panjang Data
x = x(1:3000)';

%% Preprocessing
% Centering
x = x - mean(x);
figure;
plot((1:length(x))./fs,x)
title('Sinyal PCG Asli Hasil Centering');
xlabel('Waktu (detik)');
ylabel('Amplitudo');


%% Normalisasi
x = x./max(abs(x));
figure;
plot((1:length(x))./fs,x)
title('Sinyal PCG Asli Hasil Normalisasi');
xlabel('Waktu (detik)');
ylabel('Amplitudo');

% Mengatur Panjang Data
s1 = x;

% Menambahkan Noise
L = length(s1);
SNR_dB = 10;
SNR = 10^(SNR_dB/10);
Esym = sum(abs(s1).^2)/(L);
N0 = Esym/SNR;
if (isreal(s1))
    noiseSigma = sqrt(N0);
    s2 = noiseSigma*randn(1,L);
else
     noiseSigma = sqrt(N0/2);
    s2 = noiseSigma*(randn(1,L)+1i*randn(1,L));
end

% Centering Dan Normalisasi Sinyal AWGN
s2 = s2 - mean(s2); % Centering
s2 = s2./max(abs(s2)); % Normalisasi

%% Membuat Koefisien Vektor Sinyal Campuran
A = [0.8 0.2; 0.8 0.4];
% Mencampur kedua sumber sinyal
x = A*[s1;s2];

%% Menghitung Nilai MSE sebelum proses FastICA
n1 = numel(x(1,:));
n2 = numel(x(2,:));

% % Menghitung MSE di Campuran Pertama
% mse1_1 = ((norm(s1.*0.8 - x(1,:)).^2) / n1)*100;
% 
% % Menghitung MSE di Campuran Kedua
% mse1_2 = ((norm(s2.*0.8 - x(2,:)).^2) / n2)*100;
% 
% fprintf('>> MSE Campuran Pertama Sebelum Denoising %0.4f %%\n', mse1_1);
% fprintf('>> MSE Campuran Kedua Sebelum Denoising %0.4f %%\n', mse1_2);
% 
% % snr1_1 = 20*log10(rms(s1.*0.8)/rms(noise_est1_1));
% % snr1_2 = 20*log10(rms(s2.*0.8)/rms(noise_est1_2));
% Menghitung MSE Sinyal PCG di Campuran Pertama
            n1 = numel(x(1,:));
            s1m = 0.8*s1;
            mse1_1(ix) = ((norm(s1m - x(1,:)).^2) / n1)*100;
            fprintf('>> MSE Sinyal PCG Campuran Pertama Sebelum Denoising %0.4f\n', mse1_1(ix));

            % Menghitung MSE Sinyal MSE di campuran kedua
            n2 = numel(x(2,:));
            mse1_2(ix) = ((norm(s1m - x(2,:)).^2) / n2)*100;
            fprintf('>> MSE Sinyal PCG Campuran Kedua Sebelum Denoising %0.4f\n', mse1_2(ix));
%% Menghitung nilai SNR sebelum proses Denoising menggunakan FastICA
% SNR Sebelum Denoising pada Sinyal PCG Campuran Pertama
noise_est1_1 = s1m - x(1,:);
snr1_1= 20*log10(rms(s1m)/rms(noise_est1_1));
fprintf('>> SNR Campuran 1 Sebelum Denoising %0.4f dB\n', snr1_1);
% noise_est1_1 = s1.*0.8 - x(1,:);
% p2p11 = peak2peak(noise_est1_1);
% stpd11 = std(noise_est1_1);
% cpr11 = p2p11 ./ (4*stpd11);
% snr1_1 = abs(20* log10(cpr11));

% SNR sebelum Denoising pada Sinyal PCG Campuran Kedua
noise_est1_2 = s1m - x(2,:);
snr1_2 = 20*log10(rms(s1m)/rms(noise_est1_2));
fprintf('>> SNR Campuran 2 Sebelum Denoising %0.4f dB\n', snr1_2);

% noise_est1_2 = s1.*0.8 - x(2,:);
% p2p12 = peak2peak(noise_est1_2);
% stpd12 = std(noise_est1_2);
% cpr12 = p2p12 ./ (4*stpd12);
% snr1_2 = abs(20 * log10(cpr12));

% fprintf('>> The Signal Noise to ratio Before Denoising in Source 1 is %0.4f dB\n', snr1_1);
% fprintf('>> The Signal Noise to ratio Before Denoising in Source 2 is %0.4f dB\n', snr1_2); 


%% Melakukan Proses Whitening
cov_x = x*x'; % Mencari nilai kovariance matriks x
[eig_vec,eig_val] = eig(cov_x); % Mencari nilai eigen vector, dan eigen value dari kovariance matriks
W = (eig_val^-0.5)*eig_vec'; % Rumus mencari nilai w
P = W*x;
disp('Data Covariance sebelum proses whitening\n');
disp(x*x');
disp('Data Covariance setelah proses whitening\n');
disp(P*P');

%% Menghitung Dimensi sinyal
[num_sources,num_samples] = size(P);

%% Paramter Kontrol
max_iter = 100;
convergence_threshold = 0.00000000001;

%% FastICA main loop
% Inisiasi unmixing matrix
W = rand(num_sources,num_sources);

for source = 1:num_sources
    % Memilih satu baris dari Unmixing Matrix
    w = W(source,:)';
    nh = w'*P;
    [G, g, g_prime] = findHuber(nh);
    for iter = 1 : max_iter
        w_old = w;
        % Menyimpan nilai w terakhir untuk mengecek konvergensi
        w = P*g'; - mean(g_prime,2)*w;
        % Mengecek Ortogonalitas dengan sumber sinyal yang lain
        for i = 1:source-1
            w = w - (w'*W(i,:)')*W(i,:)';
        end
        % Projection on constraint
        w = w/norm(w);
        if norm(w_old - w) < convergence_threshold
            disp(['Source', num2str(iter),'iterations'])
            break
        end
    end
    W(source,:) = w;
end
%% Plotting Hasil Akhir
y = W*P;

%% Menghitung Nilai MSE setelah proses FastICA

% Menghitung MSE di campuran pertama
n3 = numel(y(1,:));
mse2_1 = ((norm(s1.*0.8 - y(1,:)).^2) / n3)*100;
fprintf('>> MSE Campuran Pertama Setelah Denoising %0.4f %%\n', mse2_1);

% Menghitung MSE di Campuran Kedua 
n4 = numel(y(2,:));
mse2_2 = ((norm(s2.*0.8 - y(2,:)).^2) / n4)*100;
fprintf('>> MSE Campuran Kedua Setelah Denoising %0.4f %%\n', mse2_2);

%% Menghitung nilai SNR setelah Proses FastICA
% % SNR Campuran Pertama
%   noise_est2_1 = s1.*0.8 - y(1,:);
% % snr2_1 = 20*log10(rms(s1.*0.8)/rms(noise_est2_1));
%   p2p1 = peak2peak(noise_est2_1);
%   stdp1 = std(noise_est2_1);
%   cpr1 = p2p1 ./ (4 * stdp1);
%   snr2_1 = abs(20 * log10(cpr1));
%   fprintf('>> SNR Sinyal PCG Setelah Denoising %0.4f dB\n', snr2_1);
% 
% % SNR Campuran Kedua
%   noise_est2_2 = s2.*0.8 - y(2,:);
% % snr2_2 = 20*log10(rms(s2.*0.8)/rms(noise_est2_2));
%   p2p2 = peak2peak(noise_est2_2);
%   stdp2 = std(noise_est2_2);
%   cpr2 = p2p2 ./ (4 * stdp2);
%   snr2_2 = abs(20 * log10(cpr2));
%   fprintf('>> SNR Sinyal AWGN Setelah Denoising %0.4f dB\n', snr2_2); 

 % SNR Sinyal PCG Campuran 1
            noise_est2_1 = s1m - y(1,:);
            snr2_1(ix) = 20*log10(rms(s1m)/rms(noise_est2_1));
            fprintf('>> SNR Sinyal PCG Campuran 1 Setelah Denoising  %0.4f\n', snr2_1(ix));
            
            % SNR sinyal AWGN pada Campuran 2
            noise_est2_2 = s2m - y(2,:);
            snr2_2(ix) = 20*log10(rms(s2m)/rms(noise_est2_2));
  fprintf('>> SNR Sinyal AWGN Setelah Denoising %0.4f\n', snr2_2(ix));            
% Plotting Gambar
%% Menampilkan Sinyal Asli
fprintf('1. Menampilkan Sinyal PCG Asli---\n');
%membuat folder penyimpanan hasil plotting
addpath('./plots');
bgdir = pwd;
out_folder = 'NABILA_PLOT';
if ~exist(out_folder, 'dir')
    mkdir(out_folder);
end
% membuat nama figure yang akan diplotting
out_full = fullfile(bgdir, out_folder);
ot2 = sprintf('1. Sinyal PCG Asli');
oname2 = fullfile(out_full, ot2);
% inisialisasi plotting
vlgn = 0; lgn = ''; closefig = 1; ord = 0; txt = 1;
x1 = [140/fs 1100/fs 1750/fs 2700/fs];
y1 = [0.49 0.4 0.52 0.49];
xlbl = 'Waktu (detik)'; ylbl = 'Amplitudo';
tlt = 'Sinyal PCG Asli';
str = {'\bf s1', '\bf s2', '\bf s1', '\bf s2'};
% inisialisasi nilai sumbu x dan y
x2 = (1:length(s1(550:3500)))/fs;
y2 = s1(550:3500);
% fungsi plotting figure
Plot_Signal(x2, y2, oname2, vlgn, lgn, xlbl, ylbl, tlt, ord, closefig, txt, x1, y1, str);


fprintf('2. Sinyal AWGN--\n');
% membuat nama figure yang akan diplotting
out_full = fullfile(bgdir, out_folder);
ot2 = sprintf('2. Sinyal AWGN');
oname2 = fullfile(out_full, ot2);
% inisialisasi plotting
vlgn = 0; lgn = ''; closefig = 1; ord = 0; txt = 0;
x1 = [140/fs 1100/fs 1750/fs 2700/fs];
y1 = [0.49 0.4 0.52 0.49];
xlbl = 'Waktu (detik)'; ylbl = 'Amplitudo';
tlt = 'Sinyal AWGN';
str = {'\bf s1', '\bf s2', '\bf s1', '\bf s2'};
% inisialisasi nilai sumbu x dan y
x2 = (1:length(s2(550:3500)))/fs;
y2 = s2(550:3500);
% fungsi plotting figure
Plot_Signal(x2, y2, oname2, vlgn, lgn, xlbl, ylbl, tlt, ord, closefig, txt, x1, y1, str);




fprintf('3. Mixed 1--\n');
% membuat nama figure yang akan diplotting
out_full = fullfile(bgdir, out_folder);
ot2 = sprintf('3. Mixed 1');
oname2 = fullfile(out_full, ot2);
% inisialisasi plotting
vlgn = 0; lgn = ''; closefig = 1; ord = 0; txt = 1;
x1 = [140/fs 1100/fs 1750/fs 2700/fs];
y1 = [0.49 0.38 0.49 0.49];
xlbl = 'Waktu (detik)'; ylbl = 'Amplitudo';
tlt = 'Mixed 1';
str = {'\bf s1', '\bf s2', '\bf s1', '\bf s2'};
% inisialisasi nilai sumbu x dan y
x_mixed = x(1, :);
x2 = (1:length(x_mixed(550:3500)))/fs;
y2 = x_mixed(550:3500);
% fungsi plotting figure
Plot_Signal(x2, y2, oname2, vlgn, lgn, xlbl, ylbl, tlt, ord, closefig, txt, x1, y1, str);

fprintf('4. Mixed 2--\n');
% membuat nama figure yang akan diplotting
out_full = fullfile(bgdir, out_folder);
ot2 = sprintf('4. Mixed 2');
oname2 = fullfile(out_full, ot2);
% inisialisasi plotting
vlgn = 0; lgn = ''; closefig = 1; ord = 0; txt = 1;
x1 = [140/fs 1100/fs 1750/fs 2700/fs];
y1 = [0.49 0.38 0.49 0.49];
xlbl = 'Waktu (detik)'; ylbl = 'Amplitudo';
tlt = 'Mixed 2';
str = {'\bf s1', '\bf s2', '\bf s1', '\bf s2'};
% inisialisasi nilai sumbu x dan y
x_mixed = x(2, :);
x2 = (1:length(x_mixed(550:3500)))/fs;
y2 = x_mixed(550:3500);
% fungsi plotting figure
Plot_Signal(x2, y2, oname2, vlgn, lgn, xlbl, ylbl, tlt, ord, closefig, txt, x1, y1, str);
 
fprintf('5. Mixed Sources Scatter Plot');
% membuat nama figure yang akan diplotting
out_full = fullfile(bgdir, out_folder);
ot2 = sprintf('5. Mixed Sources Scatter Plot');
oname2 = fullfile(out_full, ot2);
% plot figure
figure('Color', 'White');
scatter(x(1,:),x(2,:)); 
pbaspect([3.8 2 1]);
set(gca, 'XColor', 'black', 'YColor', 'black', 'LineWidth', 1.2, 'GridAlpha', 0.1);
set(gca, 'FontSize', 11, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
grid on;
title('Mixed Sources Scatter Plot', 'FontName', 'Rockwell', 'FontSize', 12, 'FontWeight', 'normal');
ylabel('Mixed 2 Amplitude Signal', 'FontName', 'Rockwell', 'FontSize', 12, 'FontWeight', 'normal');
xlabel('Mixed 1 Amplitude Signal', 'FontName', 'Rockwell', 'FontSize', 12, 'FontWeight', 'normal');
export_fig (oname2, '-jpg', '-r200', '-a4', '-painters', '-transparent');
close;

%Plotting sinyal whitening
fprintf('6. Whitened P1');
% membuat nama figure yang akan diplotting
out_full = fullfile(bgdir, out_folder);
ot2 = sprintf('6. Whitened P1');
oname2 = fullfile(out_full, ot2);
% inisialisasi plotting
vlgn = 0; lgn = ''; closefig = 1; ord = 0; txt = 0;
x1 = [140/fs 1100/fs 1750/fs 2700/fs];
y1 = [0.49 0.38 0.49 0.49];
xlbl = 'Waktu (detik)'; ylbl = 'Amplitudo';
tlt = 'Whitened P1';
str = {'\bf s1', '\bf s2', '\bf s1', '\bf s2'};
% inisialisasi nilai sumbu x dan y
P1 = P(1,:);
x2 = (1:length(P1(550:3500)))/fs;
y2 = P1(550:3500);
% fungsi plotting figure
Plot_Signal(x2, y2, oname2, vlgn, lgn, xlbl, ylbl, tlt, ord, closefig, txt, x1, y1, str);

fprintf('7. Whitened P2');
% membuat nama figure yang akan diplotting
out_full = fullfile(bgdir, out_folder);
ot2 = sprintf('7. Whitened P2');
oname2 = fullfile(out_full, ot2);
% inisialisasi plotting
vlgn = 0; lgn = ''; closefig = 1; ord = 0; txt = 1;
x1 = [140/fs 1100/fs 1750/fs 2700/fs];
y1 = [0.12 0.11 0.11 0.11];
xlbl = 'Waktu (detik)'; ylbl = 'Amplitudo';
tlt = 'Whitened P2';
str = {'\bf s1', '\bf s2', '\bf s1', '\bf s2'};
% inisialisasi nilai sumbu x dan y
P2 = P(2,:);
x2 = (1:length(P2(550:3500)))/fs;
y2 = P2(550:3500);
% fungsi plotting figure
Plot_Signal(x2, y2, oname2, vlgn, lgn, xlbl, ylbl, tlt, ord, closefig, txt, x1, y1, str);

fprintf('8. Whitened Scatter Plot');
% membuat nama figure yang akan diplotting
out_full = fullfile(bgdir, out_folder);
ot2 = sprintf('8. Whitened Scatter Plot');
oname2 = fullfile(out_full, ot2);
% plot figure
figure('Color', 'White');
scatter(P(1,:),P(2,:));
pbaspect([3.8 2 1]);
set(gca, 'XColor', 'black', 'YColor', 'black', 'LineWidth', 1.2, 'GridAlpha', 0.1);
set(gca, 'FontSize', 11, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
grid on;
title('Whitened Scatter Plot', 'FontName', 'Rockwell', 'FontSize', 12, 'FontWeight', 'normal');
ylabel('Whitenned Amplitude x2', 'FontName', 'Rockwell', 'FontSize', 12, 'FontWeight', 'normal');
xlabel('Whitenned Amplitude x1', 'FontName', 'Rockwell', 'FontSize', 12, 'FontWeight', 'normal');
export_fig (oname2, '-jpg', '-r200', '-a4', '-painters', '-transparent');
close;

fprintf('9. Sinyal PCG Hasil Denoising');
% membuat nama figure yang akan diplotting
out_full = fullfile(bgdir, out_folder);
ot2 = sprintf('9. Sinyal PCG Hasil Denoising');
oname2 = fullfile(out_full, ot2);
% inisialisasi plotting
vlgn = 0; lgn = ''; closefig = 1; ord = 0; txt = 1;
x1 = [140/fs 1100/fs 1750/fs 2700/fs];
y1 = [0.12 0.11 0.11 0.11];
xlbl = 'Waktu (detik)'; ylbl = 'Amplitudo';
tlt = 'Sinyal PCG Hasil Denoising';
str = {'\bf s1', '\bf s2', '\bf s1', '\bf s2'};
% inisialisasi nilai sumbu x dan y
denoised = y(1,:);
x2 = (1:length(denoised(550:3500)))/fs;
y2 = denoised(550:3500);
% fungsi plotting figure
Plot_Signal(x2, y2, oname2, vlgn, lgn, xlbl, ylbl, tlt, ord, closefig, txt, x1, y1, str);

fprintf('10. Sinyal AWGN Hasil DenoisinG');
% membuat nama figure yang akan diplotting
out_full = fullfile(bgdir, out_folder);
ot2 = sprintf('10. Sinyal AWGN Hasil Denoising');
oname2 = fullfile(out_full, ot2);
% inisialisasi plotting
vlgn = 0; lgn = ''; closefig = 1; ord = 0; txt = 1;
x1 = [140/fs 1100/fs 1750/fs 2700/fs];
y1 = [0.12 0.11 0.11 0.11];
xlbl = 'Waktu (detik)'; ylbl = 'Amplitudo';
tlt = 'Sinyal AWGN Hasil Denoising';
str = {'\bf s1', '\bf s2', '\bf s1', '\bf s2'};
% inisialisasi nilai sumbu x dan y
denoised = y(2,:);
x2 = (1:length(denoised(550:3500)))/fs;
y2 = denoised(550:3500);
% fungsi plotting figure
Plot_Signal(x2, y2, oname2, vlgn, lgn, xlbl, ylbl, tlt, ord, closefig, txt, x1, y1, str);

fprintf('11. Source Scatter Plot Denoised Signal');
% membuat nama figure yang akan diplotting
out_full = fullfile(bgdir, out_folder);
ot2 = sprintf('10. Source Scatter Plot Denoised Signal');
oname2 = fullfile(out_full, ot2);
% plot figure
figure('Color', 'White');
scatter(y(1,:), y(2,:));
pbaspect([3.8 2 1]);
set(gca, 'XColor', 'black', 'YColor', 'black', 'LineWidth', 1.2, 'GridAlpha', 0.1);
set(gca, 'FontSize', 11, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
grid on;
title('Denoised Scatter Plot', 'FontName', 'Rockwell', 'FontSize', 12, 'FontWeight', 'normal');
ylabel('Sinyal AWGN', 'FontName', 'Rockwell', 'FontSize', 12, 'FontWeight', 'normal');
xlabel('Sinyal Jantung', 'FontName', 'Rockwell', 'FontSize', 12, 'FontWeight', 'normal');
export_fig (oname2, '-jpg', '-r200', '-a4', '-painters', '-transparent');
close;