function dZ = pack_state_rates(V, dqdt, pili_rates, N_segments)
%PACK_STATE_RATES  Assemble global state rate vector
%
%  Description:
%    Packs the individual velocity components (body translation, quaternion
%    rate, and per-pili Euler angle rates) into the single flat state rate
%    vector dZ required by ode15s.  Inverse operation of unpack_state.
%
%    Layout:   dZ = [V(3); dqdt(4); dtheta_1(N1); dphi_1(N1); dpsi_1(N1);
%                    ...; dtheta_M(NM); dphi_M(NM); dpsi_M(NM)]
%
%  Syntax:
%    dZ = pack_state_rates(V, dqdt, pili_rates, N_segments)
%
%  Inputs:
%    V           - [3 x 1] Body translational velocity
%    dqdt        - [4 x 1] Quaternion time derivative
%    pili_rates  - {M x 1} Cell of structs (.dtheta, .dphi, .dpsi)
%    N_segments  - [M x 1] Segments per pilus
%
%  Outputs:
%    dZ - [n x 1] Global state rate vector
%
%  Author:  Lucien (lucien-6@qq.com)
%  Created: 2026-02-06
%  Updated: 2026-02-11 — Standardised formatting

n_total = 7 + 3 * sum(N_segments);
dZ = zeros(n_total, 1);

dZ(1:3) = V;
dZ(4:7) = dqdt;

M   = length(N_segments);
idx = 8;

for k = 1:M
    Nk = N_segments(k);
    dZ(idx       : idx+Nk-1)   = pili_rates{k}.dtheta;
    dZ(idx+Nk    : idx+2*Nk-1) = pili_rates{k}.dphi;
    dZ(idx+2*Nk  : idx+3*Nk-1) = pili_rates{k}.dpsi;
    idx = idx + 3*Nk;
end

end
