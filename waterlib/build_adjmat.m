% Copyright (c) 2016,2018 Paul Irofti <paul@irofti.net>
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

function M = build_adjmat(network)
%% Build network weighted graph
% INPUTS:
%   network -- EPANET network input file (inp)
%
% OUTPUTS:
%   M -- adjacency matrix: each none zero entry constitutes a link whose
%        who's value represents the link's weight
%--------------------------------------------------------------------------
    d = epanet(network);

    % The number of junctions in a network equals the number of nodes minus
    % the number of tanks and reservoirs.
    junctions =  double(d.getNodeCount() - d.getNodeTankReservoirCount());
    M = zeros(junctions,junctions);

    links = d.getLinkCount();
    link2nodes = d.getLinkNodesIndex();

    for l=1:links
        %fprintf('link %d: ', l);
        % If link is not a pipe, skip it.
        if d.getLinkTypeIndex(l) ~= 1
            fprintf('not a pipe\n');
            continue;
        end

        i = link2nodes(l,1);
        j = link2nodes(l,2);
        % Found a pipe that connects junction to reservoir, skip it.
        if i > junctions || j > junctions
            %fprintf('does not connect two junctions\n');
            continue;
        end
        
        % Calculate link's weight.
        len = d.getLinkLength(l);
        diam = d.getLinkDiameter(l);
        %M(i,j) = len * pi * (diam/2)^2;
        M(i,j) = len;
        M(j,i) = M(i,j);
        
        %fprintf('%d <-> %d len %d diam %d\n', i, j, len, diam);
    end

    d.closeNetwork();
end
