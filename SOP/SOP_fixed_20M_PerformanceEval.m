
clc
clear all
% close all
format long
addpath('..','../NXPlib','../CW')
config = tx11b_config();
Gain_control = 0.2;
Preamble_20b = [1 0 -1 0.5 1 1 0.5 -0.5 0 1 1 1 1 0.5 -0.5 -1 -1 -1 -1 -0.5]/2;
plusplus = r_half([Preamble_20b,-Preamble_20b]*r_half(0.1)); %plusplus = [fliplr(Preamble_11b) fliplr(Preamble_11b)];    
load dnlos-ieee
SNRdB = 0; 
CFO = 0; %125e3;c
[tx_out,source_data,shortSYNC] = tx11b(config);
SopDetSymbol = 20;
Beta2E = r_half(0.55^2*plusplus*plusplus');
for loopSNR = 1:length(SNRdB)
    SNR = SNRdB(loopSNR);
    for loopParam = 1
        rng(1)
        OnTime_detect = 0;
        False_detect = 0;
        Missed_detect = 0;
        Len = 24;
        for pckt = 1:1000
            tau = round(1000*rand)+100;  
            SOP_ind = 1e6;
            if 1
                chan = squeeze(timeDomainChannel(1,1,:,pckt));
                tx160 = resample(tx_out,160,80);
                mp_out = filter(chan, 1, tx160);
                mp_out = resample(mp_out,80,160);   
            else
                mp_out = tx_out;
            end
            mp_out = [zeros(1,tau),mp_out];
            cfo_in = awgn(mp_out,SNR,'measured');
            Fs = config.fs_out;
            L = length(cfo_in);
            % ----------------- add CFO
            CFO = 150e3 - rand * 300e3;
            RxData = cfo_in.*exp(1i*2*pi*(1/Fs)*CFO*(0:L-1));            
             % ----------------- Resample data to 20MSPS             
            RxData = resample(RxData, 1, 4);
            RxData20 = r_half(0.5*RxData/max(abs(RxData)));
%             a(1:2:10*128) = real(RxData20(1:10*64));
%             a(2:2:10*128) = imag(RxData20(1:10*64));
%             [a.',ones(length(a),1)]

            RxData20  = r_half(RxData20*r_half(Gain_control));
            tau = round(tau/4);
            SOP_ind = 1e6;
            SOP_flag = 0;
            for k = 1:floor(length(RxData20)/64) - 1
                rssi = zeros(1,2*Len);
                frame64 = RxData20((k-1)*64 + 1:k*64);
                %frame32 = r_half((frame64(1:2:end) + frame64(2:2:end)));
                for kk = 1:40
                    rssi(1:2:Len*2) = rssi(1:2:Len*2) + real(frame64(kk:kk+(Len - 1))).^2;
                    rssi(2:2:Len*2) = rssi(2:2:Len*2) + imag(frame64(kk:kk+(Len - 1))).^2;
                end
                rssi = (rssi(1:2:end-1) + r_half(rssi(2:2:end)));      
                rssi = r_half(rssi);
                xcorr_pp = zeros(1,Len * 2);
                for kk = 1:40
                    xcorr_pp(1:2:Len * 2) = xcorr_pp(1:2:Len * 2) + real(frame64(kk:kk+(Len - 1)))*plusplus(kk);
                    xcorr_pp(2:2:Len * 2) = xcorr_pp(2:2:Len * 2) + imag(frame64(kk:kk+(Len - 1)))*plusplus(kk);
                end
                xcorr_pp2 = r_half(xcorr_pp).*r_half(xcorr_pp);
                xcorr_pp2 = (xcorr_pp2(1:2:end-1) + r_half(xcorr_pp2(2:2:end))).';
                Xcorr = r_half(r_half(xcorr_pp2).*r_rcp(rssi.'));
%                 dec2q(Xcorr)
                if max(Xcorr) > Beta2E
                   SOP_flag = 1;
                   SOP_ind = (k-1)*64;
                   break;
                end
            end  
            if SOP_ind >= (tau - 44) && SOP_ind <= (tau + SopDetSymbol*20)
               OnTime_detect = OnTime_detect + 1;
            elseif SOP_ind < (tau - 44)
                False_detect = False_detect + 1;
            elseif SOP_ind > (tau + SopDetSymbol*20)
                Missed_detect = Missed_detect + 1;
            end
        end
        Pd(loopSNR,loopParam) = OnTime_detect/pckt
        Pf(loopSNR,loopParam) = False_detect/pckt
        PM(loopSNR,loopParam) = Missed_detect/pckt
    end
    
end
subplot(1,2,1)
semilogy(SNRdB,Pd,'-o','LineWidth',2)
hold all
semilogy(SNRdB,PM,'-o','LineWidth',2)
xlabel('SNR(dB)')
ylabel('Probability')
legend('Detection','Miss Detection')
title('802.11b frame')
grid on