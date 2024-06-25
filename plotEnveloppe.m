function plotEnveloppe(h,X,Y,YE,c,nz,sty)
%plotEnveloppe plots confidence enveloppe around curves.
%   plotEnveloppe(h,X,Y,YE,c,nz,sty) is an utility function.
%
%   Copyright (c) 2024, Benjamin Zoller
%   All rights reserved.
%
%   This source code is licensed under the MIT license found in the
%   LICENSE file in the root directory of this source tree. 

if nargin < 6
    nz = 0;
    sty = '-';
end
if nargin < 7
    sty = '-';
end
I = ~isnan(Y) & ~isinf(Y);
X = X(I);
Y = Y(I);
YE = YE(:,I);

if size(YE,1)>1
    Ylow = YE(1,:);
    Yup = YE(2,:);
else
    Ylow = Y-YE;
    Yup = Y+YE;
end

if nz == 1
    Ylow(Ylow<1e-16) = 1e-16;
    Yup(Yup<1e-16) = 1e-16;
end

x = [X,X(end:-1:1)];
y = [Ylow,Yup(end:-1:1)];
patch(h,x,y,c,'FaceAlpha',0.3,'LineStyle','none');
plot(h,X,Y,sty,'color',c,'linewidth',2)
end

