function [k1,k2] = makeKernel(K1,K2,Ke,dt)
%makeKernel discretizes the deconvolution kernel.
%  [k1,k2] = makeKernel(K1,K2,Ke,dt) is an utility function.
%
%   Copyright (c) 2024, Benjamin Zoller
%   All rights reserved.
%
%   This source code is licensed under the MIT license found in the
%   LICENSE file in the root directory of this source tree. 

kt = K1.X/Ke;
kk = K1.K;

Nk = round((kt(end))/dt);
k1 = interp1(kt,kk,dt*(0:Nk),'previous',1);

if ~isempty(K2)
    kt = (K2.X+K1.X(end)-K2.X(end))/Ke;
    kk = K2.K;
    
    k2 = interp1(kt,kk,dt*(0:Nk),'previous',1);
    k = find(k2<1,1,'first');
    k2(1:(k-1)) = 0;
else
    k2 = [];
end
end

