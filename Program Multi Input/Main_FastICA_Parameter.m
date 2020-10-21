%% Program Utama untuk FastICA Multi Input (Parameter)
% Nama : Kharisma Putri Nabila
% NIM : 16101174

%% Environment
clear all;
close all;
clc

%% Memilih folder untuk menyimpan
direk =  uigetdir('Choose a folder where you store the data');

if ~isequal(direk, 0)
    
    Nfiles = dir(fullfile(direk, '*.wav'));
    for ix = 1: numel(Nfiles)
            
            % Import data ke Matlab
            fname = Nfiles(ix).name;
            dname = fullfile(direk, fname);
            [x, fs] = audioread(dname);
            
            fprintf('%d) File: %s', ix, fname);
            
            %% Preprocessing
            % Centering
            xbar = (sum(x))/length(x);
            x = x - xbar;
            
            % Normalisasi
            x = x./max(abs(x));

            % Mengatur Panjang Data
            s1 = x(1:10000)';

            %% Menambahkan Noise
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

            %% Membuat Random Vector
            A = [0.8 0.2; 0.8 0.4];
           
            % Mencampur kedua sumber sinyal
            x = A*[s1;s2];

            %% Menghitung Nilai MSE Sebelum Proses Denoising dengan FastICA

            % Menghitung MSE Sinyal PCG di Campuran 1
            n1 = numel(x(1,:));
            s1m = 0.8*s1;
            mse1_1(ix) = ((norm(s1m - x(1,:)).^2) / n1)*100;
            fprintf('>> MSE Sinyal PCG Campuran 1 Sebelum Denoising %0.4f\n', mse1_1(ix));

            % Menghitung MSE Sinyal MSE di campuran kedua
            n2 = numel(x(2,:));
            mse1_2(ix) = ((norm(s1m - x(2,:)).^2) / n2)*100;
            fprintf('>> MSE Sinyal PCG Campuran 2 Sebelum Denoising %0.4f\n', mse1_2(ix));


            %% Menghitung nilai SNR Sebelum Proses FastICA
            % SNR Sinyal PCG Campuran 1
            noise_est1_1 = s1m - x(1,:);
            snr1_1(ix) = 20*log10(rms(s1m)/rms(noise_est1_1));
            fprintf('>> SNR Sinyal PCG Campuran 1 Sebelum Denoising  %0.4f\n', snr1_1(ix));
            
            % SNR Sinyal PCG Campuran 2
            noise_est1_2 = s1m - x(2,:);
            snr1_2(ix) = 20*log10(rms(s1m)/rms(noise_est1_2));
            fprintf('>> SNR Sinyal PCG Campuran 2 Sebelum Denoising %0.4f\n', snr1_2(ix));
            

            %% Melakukan Proses Whitening
            cov_x = x*x'; % Mencari nilai kovariance matriks x
            [eig_vec,eig_val] = eig(cov_x); % Mencari nilai eigen vector, dan eigen value dari kovariance matriks
            W = (eig_val^-0.5)*eig_vec'; % Rumus untuk mencari nilai w
            P = W*x;
            disp('Data Covariance sebelum proses whitening');
            disp(x*x');
            disp('Data Covariance setelah proses whitening');
            disp(P*P');

            %% Mengatur Parameter Kontrol
            convergence_threshold = 0.0001; %batas konvergensi sinyal
            max_iter = 10;  % Jumlah Maksimum perulangan

            % Ng_Entropy estimation functions.
            a = 1;

            % Melakukan Perhitungan Dimensi sinyal
            [num_sources,num_samples] = size(x); 

            %% FastICA 
            W = rand(num_sources,num_sources); % Initialize unmixing matrix

            for source = 1:num_sources
                w = W(source,:)'; % Selecting one row on of Unmixing matrix
                nh = w'*P;
                [G, g, g_prime] = findHuber(nh);
                for iter = 1:max_iter
                    w_old = w;  % Save last w to check convergence
                    w = P*g' - mean(g_prime,2)*w;
                    % Orthogonality check with other sources
                    for i = 1:source-1
                        w  = w - (w'* W(i,:)')*W(i,:)';
                    end
                    % Projection on constraint 
                    w = w/norm(w);
                     if norm(w_old - w) < convergence_threshold
                         disp(['Source ',num2str(source),' Found in ',num2str(iter),' iterations'])
                         break
                     end
                end
                W(source,:) = w; 
            end
            y = W*P;
    %% Menghitung Nilai MSE setelah proses FastICA

            % Menghitung MSE di Campuran Pertama
            n3 = numel(y(1,:));
            mse2_1(ix) = ((norm(s1m - y(1,:)).^2) / n3)*100;
            fprintf('>> MSE Campuran Pertama Setelah Denoising %0.4f\n', mse2_1(ix));

            % Menghitung MSE di Campuran Kedua
            n4 = numel(y(2,:));
            s2m = 0.4*s2;
            mse2_2(ix) = ((norm(s2m - y(2,:)).^2) / n4)*100;
            fprintf('>> MSE Campuran Kedua Setelah Denoising %0.4f\n', mse2_2(ix));


            %% Menghitung nilai SNR setelah Proses FastICA
            % SNR sinyal jantung pada Campuran 1
%             noise_est2_1 = s1m - y(1,:);
%             p2p1 = peak2peak(noise_est2_1);
%             stdp1 = std(noise_est2_1);
%             cpr1 = p2p1 ./ (4 * stdp1);
%             snr2_1(ix) = abs(20 * log10(cpr1));
%             fprintf('>> SNR Sinyal PCG Setelah Denoising %0.4f\n', snr2_1(ix));
            % SNR Sinyal PCG Campuran 1
            noise_est2_1 = s1m - y(1,:);
            snr2_1(ix) = 20*log10(rms(s1m)/rms(noise_est2_1));
            fprintf('>> SNR Sinyal PCG Campuran 1 Setelah Denoising  %0.4f\n', snr2_1(ix));
            
            % SNR sinyal AWGN pada Campuran 2
            noise_est2_2 = s2m - y(2,:);
            snr2_2(ix) = 20*log10(rms(s2m)/rms(noise_est2_2));
%           p2p2 = peak2peak(noise_est2_2);
%           stdp2 = std(noise_est2_2);
%           cpr2 = p2p2 ./ (4 * stdp2);
%           snr2_2(ix) = abs(20 * log10(cpr2));
            fprintf('>> SNR Sinyal AWGN Setelah Denoising %0.4f\n', snr2_2(ix));            
    end
end

 
%% Plotting Sinyal Data Parameter
    
    addpath('./plots');
    bgdir = pwd;
    out_folder = 'Data_Parameter_FastICA';
    if ~exist(out_folder, 'dir');
        mkdir(out_folder);
    end
    
    %% SNR sinyal PCG Campuran 1 Sebelum Denoising
    out_full = fullfile(bgdir, out_folder);
    ot2 = sprintf('1. SNR sinyal PCG Campuran 1 Sebelum Denoising');
    oname2 = fullfile(out_full, ot2);
    % inisialisasi nilai sumbu x dan y
    x_snr_1 = 1 : numel(Nfiles);
    y_snr_1 = snr1_1; 
    figure('Color', 'white');
    plot( x_snr_1,  y_snr_1, 'LineWidth', 2, 'LineStyle', ':', 'Marker', 'o', 'MarkerEdgeColor', 'red', 'MarkerFaceColor', 'red', 'MarkerSize', 6);
    pbaspect([3.6 2 1]);
    set(gca, 'XColor', 'black', 'YColor', 'black', 'LineWidth', 1.5, 'GridAlpha', 0.1);
    set(gca, 'FontSize', 12, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
    grid on;
    xlabel('Data PCG', 'FontName', 'Rockwell', 'FontSize', 12, 'FontWeight', 'normal');
    ylabel('SNR (dB)', 'FontName', 'Rockwell', 'FontSize', 12, 'FontWeight', 'normal');
    title('SNR Sinyal PCG pada Campuran 1 Sebelum Denoising ', 'FontName', 'Rockwell', 'FontSize', 12, 'FontWeight', 'normal');
    export_fig (oname2, '-jpg', '-r200', '-a4', '-painters', '-transparent');
    close;
    
    %% SNR Campuran 2 sebelum denoising
    out_full = fullfile(bgdir, out_folder);
    ot2 = sprintf('2. SNR Sinyal PCG Campuran 2 Sebelum Denoising');
    oname2 = fullfile(out_full, ot2);
    % inisialisasi nilai sumbu x dan y
    x_snr_2 = 1 : numel(Nfiles);
    y_snr_2 = snr1_2; 
    figure('Color', 'white');
    plot(x_snr_2, y_snr_2, 'LineWidth', 2, 'LineStyle', ':', 'Marker', 'o', 'MarkerEdgeColor', 'red', 'MarkerFaceColor', 'red', 'MarkerSize', 6);
    pbaspect([3.6 2 1]);
    set(gca, 'XColor', 'black', 'YColor', 'black', 'LineWidth', 1.5, 'GridAlpha', 0.1);
    set(gca, 'FontSize', 12, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
    grid on;
    xlabel('Data PCG', 'FontName', 'Rockwell', 'FontSize', 12, 'FontWeight', 'normal');
    ylabel('SNR (dB)', 'FontName', 'Rockwell', 'FontSize', 12, 'FontWeight', 'normal');
    title('SNR Sinyal PCG Campuran 2 Sebelum Denoising', 'FontName', 'Rockwell', 'FontSize', 12, 'FontWeight', 'normal');
    export_fig (oname2, '-jpg', '-r200', '-a4', '-painters', '-transparent');
    close;
    
    %% MSE Sinyal PCG Sebelum Denoising Campuran 1
    out_full = fullfile(bgdir, out_folder);
    ot2 = sprintf('3. MSE Campuran 1 Sebelum Denoising');
    oname2 = fullfile(out_full, ot2);
    % inisialisasi nilai sumbu x dan y
    x_mse_1 = 1 : numel(Nfiles);
    y_mse_1 = mse1_1; 
    figure('Color', 'white');
    plot(x_mse_1, y_mse_1, 'LineWidth', 2, 'LineStyle', ':', 'Marker', 'o', 'MarkerEdgeColor', 'red', 'MarkerFaceColor', 'red', 'MarkerSize', 6);
    pbaspect([3.6 2 1]);
    set(gca, 'XColor', 'black', 'YColor', 'black', 'LineWidth', 1.5, 'GridAlpha', 0.1);
    set(gca, 'FontSize', 12, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
    grid on;
    xlabel('Data PCG', 'FontName', 'Rockwell', 'FontSize', 12, 'FontWeight', 'normal');
    ylabel('MSE', 'FontName', 'Rockwell', 'FontSize', 12, 'FontWeight', 'normal');
    title('MSE Campuran 1 Sebelum Denoising', 'FontName', 'Rockwell', 'FontSize', 12, 'FontWeight', 'normal');
    export_fig (oname2, '-jpg', '-r200', '-a4', '-painters', '-transparent');
    close;
    
    %% MSE Sinyal PCG Sebelum Denoising Campuran 2
    out_full = fullfile(bgdir, out_folder);
    ot2 = sprintf('4. MSE Campuran 2 Sebelum Denoising');
    oname2 = fullfile(out_full, ot2);
    % inisialisasi nilai sumbu x dan y
    x_mse_2 = 1 : numel(Nfiles);
    y_mse_2 = mse1_2; 
    figure('Color', 'white');
    plot(x_mse_2, y_mse_2, 'LineWidth', 2, 'LineStyle', ':', 'Marker', 'o', 'MarkerEdgeColor', 'red', 'MarkerFaceColor', 'red', 'MarkerSize', 6);
    pbaspect([3.6 2 1]);
    set(gca, 'XColor', 'black', 'YColor', 'black', 'LineWidth', 1.5, 'GridAlpha', 0.1);
    set(gca, 'FontSize', 12, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
    grid on;
    xlabel('Data PCG', 'FontName', 'Rockwell', 'FontSize', 12, 'FontWeight', 'normal');
    ylabel('MSE', 'FontName', 'Rockwell', 'FontSize', 12, 'FontWeight', 'normal');
    title('MSE Campuran 2 Sebelum Denoising', 'FontName', 'Rockwell', 'FontSize', 12, 'FontWeight', 'normal');
    export_fig (oname2, '-jpg', '-r200', '-a4', '-painters', '-transparent');
    close;
    
    %% Hasil SNR Sinyal PCG Setelah Denoising
    out_full = fullfile(bgdir, out_folder);
    ot2 = sprintf('5. SNR Sinyal PCG Hasil Denoising');
    oname2 = fullfile(out_full, ot2);
    % inisialisasi nilai sumbu x dan y
    x_snr_after_1 = 1 : numel(Nfiles);
    y_snr_after_1 = snr2_1; 
    figure('Color', 'white');
    plot(x_snr_after_1,  y_snr_after_1, 'LineWidth', 2, 'LineStyle', ':', 'Marker', 'o', 'MarkerEdgeColor', 'red', 'MarkerFaceColor', 'red', 'MarkerSize', 6);
    pbaspect([3.6 2 1]);
    set(gca, 'XColor', 'black', 'YColor', 'black', 'LineWidth', 1.5, 'GridAlpha', 0.1);
    set(gca, 'FontSize', 12, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
    grid on;
    xlabel('Data PCG', 'FontName', 'Rockwell', 'FontSize', 12, 'FontWeight', 'normal');
    ylabel('SNR (dB)', 'FontName', 'Rockwell', 'FontSize', 12, 'FontWeight', 'normal');
    title('SNR Sinyal PCG Hasil Denoising', 'FontName', 'Rockwell', 'FontSize', 12, 'FontWeight', 'normal');
    export_fig (oname2, '-jpg', '-r200', '-a4', '-painters', '-transparent');
    close;
    
    %% Hasil SNR Sinyal AWGN Setelah Denoising
    out_full = fullfile(bgdir, out_folder);
    ot2 = sprintf('6. SNR Sinyal AWGN Hasil Denoising');
    oname2 = fullfile(out_full, ot2);
    % inisialisasi nilai sumbu x dan y
    x_snr_after_2 = 1 : numel(Nfiles);
    y_snr_after_2 = snr2_2; 
    figure('Color', 'white');
    plot(x_snr_after_2, y_snr_after_2, 'LineWidth', 2, 'LineStyle', ':', 'Marker', 'o', 'MarkerEdgeColor', 'red', 'MarkerFaceColor', 'red', 'MarkerSize', 6);
    pbaspect([3.6 2 1]);
    set(gca, 'XColor', 'black', 'YColor', 'black', 'LineWidth', 1.5, 'GridAlpha', 0.1);
    set(gca, 'FontSize', 12, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
    grid on;
    xlabel('Data PCG', 'FontName', 'Rockwell', 'FontSize', 12, 'FontWeight', 'normal');
    ylabel('SNR (dB)', 'FontName', 'Rockwell', 'FontSize', 12, 'FontWeight', 'normal');
    title('SNR Sinyal AWGN Hasil Denoising', 'FontName', 'Rockwell', 'FontSize', 12, 'FontWeight', 'normal');
    export_fig (oname2, '-jpg', '-r200', '-a4', '-painters', '-transparent');
    close;
    
    %% MSE Sinyal PCG Setelah Denoising Output 1
    out_full = fullfile(bgdir, out_folder);
    ot2 = sprintf('3. MSE Keluaran 1 Setelah Denoising');
    oname2 = fullfile(out_full, ot2);
    % inisialisasi nilai sumbu x dan y
    x_mse_after_1 = 1 : numel(Nfiles);
    y_mse_after_1 = mse2_1; 
    figure('Color', 'white');
    plot( x_mse_after_1, y_mse_after_1, 'LineWidth', 2, 'LineStyle', ':', 'Marker', 'o', 'MarkerEdgeColor', 'red', 'MarkerFaceColor', 'red', 'MarkerSize', 6);
    pbaspect([3.6 2 1]);
    set(gca, 'XColor', 'black', 'YColor', 'black', 'LineWidth', 1.5, 'GridAlpha', 0.1);
    set(gca, 'FontSize', 12, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
    grid on;
    xlabel('Data PCG', 'FontName', 'Rockwell', 'FontSize', 12, 'FontWeight', 'normal');
    ylabel('MSE', 'FontName', 'Rockwell', 'FontSize', 12, 'FontWeight', 'normal');
    title('MSE Keluaran 1 Setelah Denoising', 'FontName', 'Rockwell', 'FontSize', 12, 'FontWeight', 'normal');
    export_fig (oname2, '-jpg', '-r200', '-a4', '-painters', '-transparent');
    close;
    
    %% MSE Sinyal PCG Setelah Denoising Output 2
    out_full = fullfile(bgdir, out_folder);
    ot2 = sprintf('4. MSE Keluaran 2 Setelah Denoising');
    oname2 = fullfile(out_full, ot2);
    % inisialisasi nilai sumbu x dan y
    x_mse_after_2 = 1 : numel(Nfiles);
    y_mse_after_2 = mse2_2; 
    figure('Color', 'white');
    plot(x_mse_after_2, y_mse_after_2, 'LineWidth', 2, 'LineStyle', ':', 'Marker', 'o', 'MarkerEdgeColor', 'red', 'MarkerFaceColor', 'red', 'MarkerSize', 6);
    pbaspect([3.6 2 1]);
    set(gca, 'XColor', 'black', 'YColor', 'black', 'LineWidth', 1.5, 'GridAlpha', 0.1);
    set(gca, 'FontSize', 12, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
    grid on;
    xlabel('Data PCG', 'FontName', 'Rockwell', 'FontSize', 12, 'FontWeight', 'normal');
    ylabel('MSE', 'FontName', 'Rockwell', 'FontSize', 12, 'FontWeight', 'normal');
    title('MSE Keluaran 2 Setelah Denoising', 'FontName', 'Rockwell', 'FontSize', 12, 'FontWeight', 'normal');
    export_fig (oname2, '-jpg', '-r200', '-a4', '-painters', '-transparent');
    close;
    

