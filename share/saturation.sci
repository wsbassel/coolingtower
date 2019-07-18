function [ p ] = saturation(t_k )
// estimation of saturation pressure of a given saturation temperature in degrees K
tt = t_k;
t = t_k - 273.15;
f1 = -741.9242;
    f2 = -29.721;
    f3 = -11.55286;
    f4 = -0.8685635;
    f5 = 0.1094098;
    f6 = 0.439993;
    f7 = 0.2520658;
    f8 = 0.05218684;
    tc = 374.136;
    pc = 220.88;

       ftp = 1 / (100 * tt) * (tc - t) * (f1 + f2 * (0.65 - 0.01 * t) + f3 * (0.65 - 0.01 * t) ^ 2 ...
         + f4 * (0.65 - 0.01 * t) ^ 3 + f5 * (0.65 - 0.01 * t) ^ 4 + f6 * (0.65 - 0.01 * t) ^ 5 ...
         + f7 * (0.65 - 0.01 * t) ^ 6 + f8 * (0.65 - 0.01 * t) ^ 7);
       p = pc * exp(ftp);



endfunction
//comp(saturation )
//t_k=310
//[ p ] = saturation(t_k )
//disp (p)
