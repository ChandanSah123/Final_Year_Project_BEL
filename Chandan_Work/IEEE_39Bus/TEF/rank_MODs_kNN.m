% rank_MODs_kNN.m
% Simple k-NN ranking of DB entries by KE_norm distance for a given fault & topology.
function candidates = rank_MODs_kNN(KE_norm_meas, DB, faultLoc, topologyHash, K)
% Return up to K candidate DB indices sorted by distance
dists = [];
idxs = [];
for i = 1:numel(DB)
    if ~strcmp(DB(i).fault, faultLoc)
        continue;
    end
    % optional topology filter
    if exist('topologyHash','var') && ~isempty(topologyHash)
        if ~strcmp(DB(i).topologyHash, topologyHash), continue; end
    end
    d = norm(DB(i).KE_norm(:) - KE_norm_meas(:));
    dists(end+1) = d;
    idxs(end+1) = i;
end
[sd, order] = sort(dists);
sel = idxs(order);
candidates = sel(1:min(K,numel(sel)));
end
