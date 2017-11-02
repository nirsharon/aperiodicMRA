% script name: "make_figure_counter_example"
%
% generate 4 figures for the counterexample of two different signals
% (here denoted by x and y) with the same second-order moment (and first-order too)
close all; clear;
to_save = 0;
% the signal
L   = 15;
ell = 5;    % to have integer L/ell

% creating periodic distribution (rho)
repeat_pat = rand(ell,1);
rho = repmat(repeat_pat,L/ell,1);
rho = rho/sum(rho);
C_Frho = gallery('circul',fft(rho))';

% generating the signals
x = rand(L,1).*sin((1:L)'/(2*pi)); x = x - ((L-2)/L)*mean(x);

C_x = gallery('circul',x)';
Fx = fft(x);

pm = rand(L/ell,1);
pm(pm<0.5) = -1;
pm(pm>=0.5) = 1;
pm = (-1)*ones(L/ell,1);
pm(1) =1;
sign_change = repmat(pm,ell,1);
Fy = Fx.*sign_change;

y = ifft(Fy);
C_y = gallery('circul',y)';
norm(C_x*diag(rho)*C_x' - C_y*diag(rho)*C_y','fro')

ln = 2; ms = 13; ax_fs = 22;

%norm(align_to_reference(x,y)-y)

if to_save
    folder_name = ['counter_example_figures_',date];
    mkdir(folder_name);
    cd(folder_name);
end

figure; plot(x,'LineWidth',ln); hold on; plot(y,'red','LineWidth',ln); xlim([1,L])
nameit = 'TwoSignals';
set(gca,'FontSize',ax_fs)
if to_save
    saveas(gcf,nameit,'fig');
    saveas(gcf,nameit,'jpg');
    print('-depsc2',nameit);
end

if 0
figure; imagesc(C_y*diag(rho)*C_y'); % title('second order matrix')
nameit = 'C_Fp_matrix';
set(gca,'FontSize',22)
if to_save
    saveas(gcf,nameit,'fig');
    saveas(gcf,nameit,'jpg');
    print('-depsc2',nameit);
end
end

if 0
figure; imagesc(abs(C_Frho));% title('C_{\mathcal{F}\rho}')
nameit = 'AbsCRho';
set(gca,'FontSize',22)
if to_save
    saveas(gcf,nameit,'fig');
    saveas(gcf,nameit,'jpg');
    print('-depsc2',nameit);
end
end


figure; plot(rho,'LineWidth',ln);
nameit = 'Rho';
set(gca,'FontSize',ax_fs);
xlim([1,L]);
if to_save
    saveas(gcf,nameit,'fig');
    saveas(gcf,nameit,'jpg');
    print('-depsc2',nameit);
end


figure; stem(real(fft(x)),'+','LineWidth',ln,'markersize',ms); hold on; 
stem(real(fft(y)),'s','LineWidth',ln,'markersize',ms);
nameit = 'real_fft_signals';
%title('Real part')
set(gca,'FontSize',22)
xlim([1,L]);
ymax = max(abs(real(fft(y))))*1.1;
ylim([-ymax,ymax]);
if to_save
    saveas(gcf,nameit,'fig');
    saveas(gcf,nameit,'jpg');
    print('-depsc2',nameit);
end

figure; stem(imag(fft(x)),'+','LineWidth',ln,'markersize',ms); hold on; 
stem(imag(fft(y)),'s','LineWidth',ln,'markersize',ms);
nameit = 'imag_fft_signals';
set(gca,'FontSize',ax_fs)
xlim([1,L]);
ymax = max(abs(imag(fft(y))))*1.1;
ylim([-ymax,ymax]);
if to_save
    saveas(gcf,nameit,'fig');
    saveas(gcf,nameit,'jpg');
    print('-depsc2',nameit);
end

if to_save
    save('data')
    cd '../'  
end


