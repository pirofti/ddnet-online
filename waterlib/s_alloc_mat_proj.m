% Copyright (c) 2020, Paul Irofti <paul@irofti.net>
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

function [I, err] = s_alloc_mat_proj(R,s,lambda,D)
%% Sensor Allocation using matrix approximation.
% INPUTS:
%   R -- residues matrix
%   s -- number of sensors
%   lambda -- distance penalty
%   D -- network node distances matrix
%
% OUTPUTS:
%   I -- sensor nodes
%   err -- frobenius norm of the approximation error
    [n,m] = size(R);
    m = m/n;
    %R = R(:,m*4:end);
    R = R./repmat(max(R),n,1);
    %%-------------------------------------------------------------------------
    I = [];     % selected nodes
    A = [];     % selected rows

    % Get node distances
    if nargin < 3
        lambda = 0;
        D = 0;
    end

    % Choose the first row corresponding to the largest 2-norm
    norms = sqrt(sum(R.^2,2));
    [~,i] = max(norms);
    A = [A R(i,:)'];
    A = normc(A);
    I = [I i];

    % Pick the next rows based on the smallest projection on the
    % existing selection.
    for k = 2:s   
        P = R*R(I,:)';  % projections on current vector set
        T = P*A';       % projections sum
        N = sqrt(sum(T.^2,2));  % projection norms
        if lambda
            N = N + lambda * sum(1./D(:,I),2);  % distance penalty
        end
        N(I) = Inf;     % do not pick from exiting set, is there a better way?
        [~,i] = min(N);
        a = R(i,:)' - T(i,:)';
        a = a/norm(a);
        A = [A a];
        I = [I i];
    end

    sort(I);
    R0 = R;
    R0(setdiff(1:size(R,1),I),:) = 0;
    err = norm(R0-R,'fro');
end
