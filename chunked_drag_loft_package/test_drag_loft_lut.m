function test_drag_loft_lut
% Standalone synthetic test; no Drag.m or propagator needed.
[L,M,R] = ndgrid([15 20],[1 2],[1e5 2e5]);
S = 1 + .002*(L-15) + .01*(M-1) + .02*(R-1e5)/1e5;
T = table(L(:),M(:),R(:),S(:),'VariableNames',{'Loft','Mach','Re','Scale'});
lut = make_drag_loft_lut(T);
cleanup = onCleanup(@() drag_loft_state_scale('reset'));
drag_loft_state_scale('set',lut,16.7);
s = drag_loft_state_scale(1.4,1.5e5);
expected = 1 + .002*1.7 + .01*.4 + .02*.5;
assert(abs(s-expected)<1e-10);
try
    drag_loft_state_scale(3,1.5e5);
    error('Test:ExpectedError','Out-of-domain query was accepted.');
catch err
    assert(strcmp(err.identifier,'DragLoft:OutsideCalibration'));
end
drag_loft_state_scale('reset');
assert(drag_loft_state_scale(1,1e5)==1);
fprintf('Synthetic lookup tests passed.\n');
end
