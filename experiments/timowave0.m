i:= {1, 0, 0}

u[x_,t_] := {1, 1, 1}
r[x_,t_] := {0, 0, 0}

v[x_,t_] =  D[u[x,t],t]
s[x_,t_] =  D[r[x,t], t]
n[x_,t_] = -D[u[x,t],x]-Cross[i, r[x,t]]
m[x_,t_] = -D[r[x,t], x]
f[x_,t_] =  D[v[x,t],t]+D[n[x,t],x]
g[x_,t_] =  D[s[x,t],t]+D[m[x,t],x]+Cross[i,n[x,t]]

Print["----u,r,n,m,v,s,f,g"]
Print[u[x,t]]
Print[r[x,t]]
Print[n[x,t]]
Print[m[x,t]]
Print[v[x,t]]
Print[s[x,t]]
Print[f[x,t]]
Print[g[x,t]]

P[n_,x_]  := LegendreP[n,2x-1]
Pn[n_,x_] := P[n,x]/Sqrt[Integrate[P[n,x]^2, {x,0,1}]]

prod[f_,g_]:=Integrate[f*g, {x,0,1}]

Print["----inner products at t=0 for u"]
Print[Transpose[N[Table[prod[Pn[n,x],u[x,0]],{n,0,1}]]]]
