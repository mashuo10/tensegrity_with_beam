function [Ia,Ib,a,b]=tenseg_boundary_2D(pinned_X,pinned_Y,nn)
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/.
%
% This function calculate the free and pinned number of nodal coordinate,
% as well as the transfer matrix of free and pinned coordinates
% Inputs:
%   pinned_X: number of nodes pinned in X direction
%	pinned_Y: number of nodes pinned in Y direction
%   pinned_Z: number of nodes pinned in Z direction
%
% Outputs:
%	Ia: transfer matrix of free coordinates
%	Ib: transfer matrix of pinned coordinates
%   a: number of free coordinates
%   b: number of pinned coordinates

%%
I=eye(2*nn);
b=sort([2*pinned_X-1;2*pinned_Y]);   %index of pinned nodes
a=setdiff((1:2*nn)',b);  %index of free node direction
Ia=I(:,a);  %free node index
Ib=I(:,b);  %pinned nod index
end
