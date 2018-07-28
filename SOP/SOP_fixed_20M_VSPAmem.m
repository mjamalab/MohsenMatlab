
clc
clear all
close all
format long
addpath('../NXPlib','../CW')
init_Signal20;
Gain_control = 0.2;
Preamble_20b = [1 0 -1 0.5 1 1 0.5 -0.5 0 1 1 1 1 0.5 -0.5 -1 -1 -1 -1 -0.5]/2;
plusplus = r_half([Preamble_20b,Preamble_20b]*r_half(0.1)); %plusplus = [fliplr(Preamble_11b) fliplr(Preamble_11b)];    
plusminus = r_half([Preamble_20b,-Preamble_20b]*r_half(0.1));
Len = 32;
SopDetSymbol = 20;
Beta2E = r_half(0.55^2*plusminus*plusminus');
Rx = rx20(1:2:end) + 1i*rx20(2:2:end);
RxData20  = r_half(Rx.'*r_half(Gain_control));
SOP_ind = 1e6;
SOP_flag = 0;
% temp = dmemReadHexFile('sop.hex',0);
% RxData20 = (dmemReadComplex(temp,0,64,'half_fixed')).';
for k = 1:floor(length(RxData20)/64) - 1
    rssi = zeros(1,32);
    frame64 = RxData20((k-1)*64 + 1:k*64 + 8);
    %frame32 = r_half((frame64(1:2:end) + frame64(2:2:end)));
    for kk = 1:40
        rssi = rssi + real(frame64(kk:kk+(Len-1))).^2 + imag(frame64(kk:kk+(Len-1))).^2;
    end
    rssi = r_half(rssi);
    dec2q(rssi)
    xcorr_pp = zeros(1,32);
    for kk = 1:40
        xcorr_pp = xcorr_pp + frame64(kk:kk+(Len - 1))*plusplus(kk);
    end
    xcorr_pp.';
    xcorr_pp2 = abs(r_half(xcorr_pp)).^2
    Xcorr = r_half(r_half(xcorr_pp2).*r_rcp(rssi.'));
    dec2q(Xcorr)
    if max(Xcorr) > Beta2E
       SOP_flag = 1;
       SOP_ind = (k-1)*64;
       break;
    end
end  

