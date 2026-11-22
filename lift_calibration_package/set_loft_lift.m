function [tNodes, CL] = set_loft_lift(out, loftDeg, t0, tf)
%SET_LOFT_LIFT Interpolate fitted effective lift CL between lofts.
% Uses normalized flight phase like the existing set_loft_drag.m.

[lofts,order] = sort([out.loftDeg]);
if any(diff(lofts) <= 0) || loftDeg < lofts(1) || loftDeg > lofts(end)
    error('Loft must be within the fitted range; lofts must be unique.');
end
if t0 == tf, error('Start and end times must differ.'); end

i = find(lofts <= loftDeg,1,'last');
j = find(lofts >= loftDeg,1,'first');
a = out(order(i)); b = out(order(j));

tauA = (a.knots(:)-min(a.time))/(max(a.time)-min(a.time));
tauB = (b.knots(:)-min(b.time))/(max(b.time)-min(b.time));
tau = unique([tauA;tauB]);

cA = interp1(tauA,a.knotCL(:),tau,'linear');
cB = interp1(tauB,b.knotCL(:),tau,'linear');
if i == j
    CL = cA;
else
    w = (loftDeg-lofts(i))/(lofts(j)-lofts(i));
    CL = (1-w)*cA + w*cB;
end

tNodes = t0 + (tf-t0)*tau;
lift_study_coefficient('set',tNodes,CL);
end
