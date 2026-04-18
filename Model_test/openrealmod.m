function realmod = openrealmod(asigma, dom)
    % 模型中的参数
    hA = 10; hf = 64; hp = 0.5; K = 90; 
    mA = 0.15; mf = 0.11; p = 7; rm = 0.3;
    
    Pri = 4;

    % a = V1 - k1*x; U0 = V3*(k4*x + V2)./(k2^2*k3^2 + (k4*x + V2).^2);
    d1 = @(T) Pri ./ (hp+Pri)*rm .* T .* (1 - T / K) - (mA * T * hA) ./ (T + hA) - (mf * T * hf^p)./(hf^p + T.^p);;  % 确定性部分
%     d2 = @(x) (asigma.*x);
    d2 = @(T) (asigma.*T).^2 ./ 2;   % 乘性噪声
    realmod = langevin_eq('D1', d1, 'D2', d2, 'domain', dom);
end