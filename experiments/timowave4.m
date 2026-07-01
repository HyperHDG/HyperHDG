omega := 2*Pi

(* Manufactured Timoshenko-beam solution for the planar cross (cross2.geo).
   Global fields (functions of sigma = x+y+z, restricting to the arm coordinate on
   each arm):
       u = {1,2,3} Cos[omega sigma] Cos[omega t]   (all components, nonzero at center)
       r = 0
   Convention (matches experiments/timowave.hxx):
     n = -du/ds - i x r,   m = -dr/ds,   f = dv/dt + dn/ds,   g = ds/dt + dm/ds + i x n.
   NB: a nonzero rotation r cannot be manufactured on this cross -- it makes n(0) = -i x r(0)
   arm-direction dependent, which cannot balance at the interior node (no nodal force), and
   the solver's rotational coupling is currently only 1st-order accurate in time (WIP). *)

emit[lab_, vec_] := (
  Print["-- CForm ", lab, " --"];
  Do[Print["res[", k - 1, "] = ", CForm[Simplify[vec[[k]]]], ";"], {k, 1, 3}]);

(* ============================= x-arm: i = e_x ============================= *)
i := {1, 0, 0}
u[x_, y_, z_, t_] := {1, 2, 3} Cos[omega x] Cos[omega t]
r[x_, y_, z_, t_] := {0, 0, 0}

v[x_, y_, z_, t_] =  D[u[x, y, z, t], t]
s[x_, y_, z_, t_] =  D[r[x, y, z, t], t]
n[x_, y_, z_, t_] = -D[u[x, y, z, t], x] - Cross[i, r[x, y, z, t]]
m[x_, y_, z_, t_] = -D[r[x, y, z, t], x]
f[x_, y_, z_, t_] =  D[v[x, y, z, t], t] + D[n[x, y, z, t], x]
g[x_, y_, z_, t_] =  D[s[x, y, z, t], t] + D[m[x, y, z, t], x] + Cross[i, n[x, y, z, t]]

Print["==== x-arm  u,r,n,m,v,s,f,g ===="]
Print[u[x, y, z, t]]; Print[r[x, y, z, t]]; Print[n[x, y, z, t]]; Print[m[x, y, z, t]]
Print[v[x, y, z, t]]; Print[s[x, y, z, t]]; Print[f[x, y, z, t]]; Print[g[x, y, z, t]]
emit["g x-arm (s=point[0])", g[x, y, z, t]]

(* ============================= y-arm: i = e_y ============================= *)
i := {0, 1, 0}
u[x_, y_, z_, t_] := {1, 2, 3} Cos[omega y] Cos[omega t]
r[x_, y_, z_, t_] := {0, 0, 0}

v[x_, y_, z_, t_] =  D[u[x, y, z, t], t]
s[x_, y_, z_, t_] =  D[r[x, y, z, t], t]
n[x_, y_, z_, t_] = -D[u[x, y, z, t], y] - Cross[i, r[x, y, z, t]]
m[x_, y_, z_, t_] = -D[r[x, y, z, t], y]
f[x_, y_, z_, t_] =  D[v[x, y, z, t], t] + D[n[x, y, z, t], y]
g[x_, y_, z_, t_] =  D[s[x, y, z, t], t] + D[m[x, y, z, t], y] + Cross[i, n[x, y, z, t]]

Print["==== y-arm  u,r,n,m,v,s,f,g ===="]
Print[u[x, y, z, t]]; Print[r[x, y, z, t]]; Print[n[x, y, z, t]]; Print[m[x, y, z, t]]
Print[v[x, y, z, t]]; Print[s[x, y, z, t]]; Print[f[x, y, z, t]]; Print[g[x, y, z, t]]
emit["g y-arm (s=point[1])", g[x, y, z, t]]

(* ============================= z-arm: i = e_z ============================= *)
i := {0, 0, 1}
u[x_, y_, z_, t_] := {1, 2, 3} Cos[omega z] Cos[omega t]
r[x_, y_, z_, t_] := {0, 0, 0}

v[x_, y_, z_, t_] =  D[u[x, y, z, t], t]
s[x_, y_, z_, t_] =  D[r[x, y, z, t], t]
n[x_, y_, z_, t_] = -D[u[x, y, z, t], z] - Cross[i, r[x, y, z, t]]
m[x_, y_, z_, t_] = -D[r[x, y, z, t], z]
f[x_, y_, z_, t_] =  D[v[x, y, z, t], t] + D[n[x, y, z, t], z]
g[x_, y_, z_, t_] =  D[s[x, y, z, t], t] + D[m[x, y, z, t], z] + Cross[i, n[x, y, z, t]]

Print["==== z-arm  u,r,n,m,v,s,f,g ===="]
Print[u[x, y, z, t]]; Print[r[x, y, z, t]]; Print[n[x, y, z, t]]; Print[m[x, y, z, t]]
Print[v[x, y, z, t]]; Print[s[x, y, z, t]]; Print[f[x, y, z, t]]; Print[g[x, y, z, t]]
emit["g z-arm (s=point[2])", g[x, y, z, t]]

(* ================= global primaries (sigma = x+y+z) ================= *)
uu[sig_, t_] := {1, 2, 3} Cos[omega sig] Cos[omega t]
emit["initial_u (sig=point[0]+point[1]+point[2])", uu[sig, t]]
emit["initial_v (dU/dt)", D[uu[sig, t], t]]
