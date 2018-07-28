function out = sop_detect(RxData)
Preamble_11b = [1 0 -1 0.5 1 1 0.5 -0.5 0 1 1 1 1 0.5 -0.5 -1 -1 -1 -1 -0.5];
plusminus = [(Preamble_11b) (-Preamble_11b)]/10;
Beta = .55;
PP = dot(plusminus,plusminus);
out = 0;
for k = 1:24
    corr_11b_pp = dot(plusminus,RxData(k:k+39));
    rssi = abs(dot(RxData(k:k+39),RxData(k:k+39)));
    corr_11b_mag = abs(corr_11b_pp)^2/rssi;
    Thr = Beta^2*PP;
    if (corr_11b_mag) >= Thr
       out = 1;
    break;
    end
end