function test_split_drag_runs
% Standalone test; no optimizer or propagator needed.
t = [30;43;59;61;75;92;100];
r = struct('label','synthetic','loftDeg',15,'time',t, ...
    'R',[t,2*t,3*t],'V',[4*t,5*t,6*t]);
c = split_drag_runs(r,30);
assert(numel(c)==3);
assert(isequal(vertcat(c.interval),[30 60;60 90;90 100]));
assert(isequal(c(2).R(1,:),[60 120 180]));
assert(isequal(c(2).V(1,:),[240 300 360]));
assert(isequal(c(1).R(end,:),c(2).R(1,:)));
assert(all([c.loftDeg]==15) && all([c.parentCase]==1));
for k = 1:numel(c)
    assert(all(diff(c(k).time)>0));
end
% Also check an exact final boundary and descending input timestamps.
r.time = flipud(t); r.R = flipud(r.R); r.V = flipud(r.V);
c = split_drag_runs(r,35);
assert(isequal(vertcat(c.interval),[100 65;65 30]));
assert(isequal(c(2).R(1,:),[65 130 195]));
assert(isequal(c(2).V(1,:),[260 325 390]));
assert(all(diff(c(2).time)<0));
fprintf('Chunk splitting tests passed.\n');
end
