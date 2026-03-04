i:= {1, 0, 0}

u[x_,t_] := {1t x^3, 2t x^3, 3t x^3}
r[x_,t_] := {5t x^3, 7t x^3, 11 t x^3}

v[x_,t_] =  D[u[x,t],t]
s[x_,t_] =  D[r[x,t], t]
n[x_,t_] = -D[u[x,t],x]-Cross[i, r[x,t]]
m[x_,t_] = -D[r[x,t], x]
f[x_,t_] =  D[v[x,t],t]+D[n[x,t],x]
g[x_,t_] =  D[s[x,t],t]+D[m[x,t],x]+Cross[i,n[x,t]]

Print["----u,r,n,m,f,g"]
Print[u[x,t]]
Print[r[x,t]]
Print[n[x,t]]
Print[m[x,t]]
Print[f[x,t]]
Print[g[x,t]]

P[n_,x_]  := LegendreP[n,2x-1]
Pn[n_,x_] := P[n,x]/Sqrt[Integrate[P[n,x]^2, {x,0,1}]]

Print["----orthonormal Legendre Polynomials on [0,1]"]
Print[Table[Pn[n,x],{n,0,3}]]

prod[f_,g_]:=Integrate[f*g, {x,0,1}]

Print["----inner products at t=0 for v, s"]
Print[Transpose[N[Table[prod[Pn[n,x],v[x,0]],{n,0,3}],4]]]
Print[Transpose[N[Table[prod[Pn[n,x],s[x,0]],{n,0,3}],4]]]
