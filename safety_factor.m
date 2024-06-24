function MOS = safety_factor(sigma1_, sigma2_, sigma12_, Xt, Xc, Yt, Yc, S)
    F1 = 1/Xt + 1/Xc;
    F2 = 1/Yt + 1/Yc;
    F6 = 0;
    F11 = -1/(Xt*Xc);
    F22 = -1/(Yt*Yc);
    F66 = 1/(S^2);
    F12 = -1/(2*Xt*Xc);
   
    % polynomial: a*r^2 + B*r + c
    a = F11.*sigma1_.^2 + F22.*sigma2_.^2 + F66.*sigma12_.^2 + F12.*sigma1_.*sigma2_;
    b = F1.*sigma1_ + F2.*sigma2_ + F6.*sigma12_;
    c = -1;
    
    r = roots([a, b, c]);
    R = r(real(r)>0)
    MOS = R; % Margin of Safety
end
