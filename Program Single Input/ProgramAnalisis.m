%% Program Untuk Analisis
% Nama : Kharisma Putri Nabila
% NIM : 16101174

% Environment
clear all;
close all;
clc;

%% Menginput data
[fname,pname] = uigetfile('*.wav','Pilih sebuah Folder PCG');
pcgfile = fullfile(pname,fname);
[s, fs] = audioread(pcgfile);
s = s(101:105)';

%% Preprocessing pada sinyal jantung
% Centering pada Sinyal Jantung
s1bar = sum(s)./length(s);
sc = s - s1bar;

% Proses Normalisasi pada sinyal Jantung
scmax = max(abs(sc));
s1 = sc./scmax;

%% Penambahan sinyal AWGN
s_2 = [0.1361 0.4642 -0.5717 0.2182 0.0807];

%% Preprocessing pada sinyal AWGN
s2bar = sum(s_2)./length(s_2);
sc2 = s_2 - s2bar;

scmax2 = max(abs(sc2));
s2 = sc2./scmax2;

% Penambahan Random Vektor Tercampur
A = [0.8 0.2; 0.8 0.4];
x = A*[s1;s2];


%% Melakukan Proses Whitening
cov_x = x*x'; % Mencari nilai kovariance matriks x
[eig_vec,eig_val] = eig(cov_x); % Mencari nilai eigen vector, dan eigen value dari kovariance matriks
W1 = (eig_val^-0.5)*eig_vec'; % Rumus yang di proposal nabila untuk mencari nilai w
P = W1*x;
disp('Data Covariance sebelum proses whitening');
disp(x*x');
disp('Data Covariance setelah proses whitening');
disp(P*P');

%% Mengatur Parameter Kontrol
convergence_threshold = 0.0001;
max_iter = 10;  % Jumlah Maksimum perulangan

% Ng_Entropy estimation functions.
a = 1;

% Melakukan Perhitungan Dimensi sinyal
[num_sources,num_samples] = size(x); 

%% ICA main loop
W = [0.8147 0.1270; 0.9058 0.9134]; % Initialize unmixing matrix

for source = 1:2
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