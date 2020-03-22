% Copyright (c) 2020 Florin Stoican <florin.stoican@acse.pub.ro>
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

function [selected_sensors,failures]=sensor_selection_MI(M,s)

alpha  = binvar(size(M, 1), 1);
beta   = binvar(size(M, 2), 1);
cost   = sum(alpha) + 1000 * sum(beta);
constr = [sum(alpha) == s];
% constr=[];

for j = 1:size(M,2) % numberof_rows    
    constr = [constr, M(:,j)'*alpha >= 1 - beta(j)];
end
diagnostic = optimize(constr, cost);

selected_sensors = find(double(alpha));
failures=find(double(beta));
