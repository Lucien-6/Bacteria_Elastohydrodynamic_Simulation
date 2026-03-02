function S = skew_matrix(v)
%SKEW_MATRIX  3x3 skew-symmetric (cross-product) matrix from a 3-vector
%
%  Description:
%    Constructs the 3x3 skew-symmetric matrix S(v) such that S(v)*u = v x u
%    for any vector u.  Used to express the cross-product and angular-velocity
%    operator in matrix form throughout the Q and B matrix assemblies.
%
%    Explicit form:
%        S(v) = [  0   -v3   v2 ]
%               [  v3   0   -v1 ]
%               [ -v2   v1   0  ]
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
