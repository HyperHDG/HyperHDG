i:= {1, 0, 0}

u[x_,t_] := {Cos[t], Cos[t], Cos[t]}
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
