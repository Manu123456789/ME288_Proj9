function lut = make_drag_loft_lut(samples)
% Build 3-D linear scattered interpolation from a sample table.
P = [samples.Loft,samples.Mach,samples.Re];
s = samples.Scale;
if any(~isfinite([P(:);s(:)])) || any(s<=0)
    error('Lookup samples must be finite and drag scales positive.');
end
lut.offset = min(P,[],1);
lut.span = max(P,[],1)-lut.offset;
if any(lut.span<=0), error('All three coordinates must vary.'); end
P = (P-lut.offset)./lut.span;
[P,~,g] = unique(P,'rows');
if size(P,1)<4 || rank(P-P(1,:))<3
    error('Samples do not span a 3-D volume; more independent data needed.');
end
sv = accumarray(g,s,[],@mean);
lut.duplicateSpread = max(accumarray(g,s,[],@max)-accumarray(g,s,[],@min));
lut.F = scatteredInterpolant(P(:,1),P(:,2),P(:,3),sv,'linear','none');
lut.samples = samples;
lut.lofts = unique(samples.Loft);
end
