omega := 2*Pi

(* Manufactured Timoshenko-beam solution for the 3d cross (cross2.geo, arms +-x, +-y, +-z).
   Global fields, sigma = signed axis coordinate of the arm:
       u_k = A_k Cos[omega sigma] Cos[omega t + a_k]
       r_k = B_k Cos[omega sigma] Cos[omega t + b_k]
   Convention (matches include/HyperHDG/local_solver/timowave.hxx):
     n = -du/ds - i x r,   m = -dr/ds,   f = dv/dt + dn/ds,   g = ds/dt + dm/ds + i x n.
   Nodal force balance at the center needs the arms in +-pairs: n(0) = -i x r(0) flips sign
   with the arm direction i, so opposite arms cancel and unpaired arms would not balance.
   Written with the signed global coordinate and the positive axis vector i = e_ax, the same
   f, g formulas hold on both arms of a pair (i and sin both flip sign; cos terms are even). *)

A = {1, 2, 3};  a = {a1, a2, a3};   (* u amplitudes and per-component time phases *)
B = {5, 7, 11}; b = {b1, b2, b3};   (* r amplitudes and per-component time phases *)

Uvec[t_] := Table[A[[k]] Cos[omega t + a[[k]]], {k, 3}]
Rvec[t_] := Table[B[[k]] Cos[omega t + b[[k]]], {k, 3}]

emit[lab_, vec_] := (
  Print["-- ", lab, " --"];
  Do[Print["res[", k - 1, "] = ", Simplify[vec[[k]]], ";"], {k, 1, 3}]);

Do[
  i = IdentityMatrix[3][[ax]];
  u[w_, t_] := Uvec[t] Cos[omega w];
  r[w_, t_] := Rvec[t] Cos[omega w];
  v[w_, t_] =  D[u[w, t], t];
  s[w_, t_] =  D[r[w, t], t];
  n[w_, t_] = -D[u[w, t], w] - Cross[i, r[w, t]];
  m[w_, t_] = -D[r[w, t], w];
  f[w_, t_] =  D[v[w, t], t] + D[n[w, t], w];
  g[w_, t_] =  D[s[w, t], t] + D[m[w, t], w] + Cross[i, n[w, t]];
  Print["==== arm axis e_", ax, " (w = point[", ax - 1, "]) ===="];
  emit["n", n[w, t]]; emit["m", m[w, t]];
  emit["f", f[w, t]]; emit["g", g[w, t]],
  {ax, 1, 3}]

Print["==== initial data (sig = point[0]+point[1]+point[2]) ===="]
emit["initial_u", Uvec[t] Cos[omega sig]]
emit["initial_v", D[Uvec[t] Cos[omega sig], t]]
emit["initial_r", Rvec[t] Cos[omega sig]]
emit["initial_s", D[Rvec[t] Cos[omega sig], t]]
