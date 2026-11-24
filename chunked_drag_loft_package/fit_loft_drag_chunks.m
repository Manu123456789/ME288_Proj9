function [out,chunkRuns] = fit_loft_drag_chunks(runs,dt,chunkSeconds,nKnots)
% New wrapper; original adaptive optimizer and LUT builder remain unchanged.
% nKnots is the number of initial optimization knots PER chunk.
if nargin<4, nKnots = []; end
chunkRuns = split_drag_runs(runs,chunkSeconds);
% Each chunk is an independent case for the existing optimizer.
drag_loft_state_scale('reset');
out = fit_loft_drag_time(chunkRuns,dt,nKnots);
for j = 1:numel(out)
    out(j).parentCase = chunkRuns(j).parentCase;
    out(j).chunkIndex = chunkRuns(j).chunkIndex;
    out(j).interval = chunkRuns(j).interval;
end
end
