% Copyright (c) 2016, Paul Irofti <paul@irofti.net>
% 
% Permission to use, copy, modify, and/or distribute this software for any
% purpose with or without fee is hereby granted, provided that the above
% copyright notice and this permission notice appear in all copies.
% 
% THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
% WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
% MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
% ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
% WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
% ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
% OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.


function [p, ptrans] = all_junctions_get_element(d,n,K)
%% Walk each junction within the network and fetch the pressures
% INPUTS:
%   n -- total number of junctions
%
% OUTPUTS:
%   p -- pressure values
%   ptrans -- transitive pressure values
    tstep=1;
    itrans=1;   % transitive values counter

    % 10 = Don't save data and reinit links
    d.initializeHydraulicAnalysis(10);

    % 15 minutes is the observed length of time until the next hydraulic
	% event occurs that ENnextH() returns.
    % We stop when this value is no longer positive.
    while tstep>0
        d.runHydraulicAnalysis();

        for i=1:n
            % d.getNodesInfo can get all info
            ptrans(i,itrans)=d.getNodeHydaulicHead(i);
        end
        tstep=d.nextHydraulicAnalysisStep();
        itrans = itrans + 1;
    end
    p = sum(ptrans(:,K)')./numel(K);  % final pressure values
end
