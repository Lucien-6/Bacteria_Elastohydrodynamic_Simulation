function S = skew_matrix(v)
%SKEW_MATRIX  3x3 skew-symmetric (cross-product) matrix from a 3-vector
%
%  Syntax:
%    S = skew_matrix(v)
%
%  Inputs:
%    v - [3 x 1] Vector
%
%  Outputs:
%    S - [3 x 3] Skew-symmetric matrix such that S * u = v x u
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-02-06
%  Updated: 2026-02-11 — Standardised formatting

S = [  0,    -v(3),  v(2);
      v(3),   0,    -v(1);
     -v(2),  v(1),   0  ];

end
