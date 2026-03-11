i:= {1, 0, 0}

u[x_,y_,t_] := {0, 0, t (x+y)}
r[x_,y_,t_] := {0, 0, 0}

v[x_,y_,t_] =  D[u[x,y,t],t]
s[x_,y_,t_] =  D[r[x,y,t],t]
n[x_,y_,t_] = -D[u[x,y,t],x]-Cross[i, r[x,y,t]]
m[x_,y_,t_] = -D[r[x,y,t],x]
f[x_,y_,t_] =  D[v[x,y,t],t]+D[n[x,y,t],x]
g[x_,y_,t_] =  D[s[x,y,t],t]+D[m[x,y,t],x]+Cross[i,n[x,y,t]]

Print["----u,r,n,m,v,s,f,g"]
Print[u[x,y,t]]
Print[r[x,y,t]]
Print[n[x,y,t]]
Print[m[x,y,t]]
Print[v[x,y,t]]
Print[s[x,y,t]]
Print[f[x,y,t]]
Print[g[x,y,t]]

i:= {0, 1, 0}
v[x_,y_,t_] =  D[u[x,y,t],t]
s[x_,y_,t_] =  D[r[x,y,t],t]
n[x_,y_,t_] = -D[u[x,y,t],y]-Cross[i, r[x,y,t]]
m[x_,y_,t_] = -D[r[x,y,t],y]
f[x_,y_,t_] =  D[v[x,y,t],t]+D[n[x,y,t],y]
g[x_,y_,t_] =  D[s[x,y,t],t]+D[m[x,y,t],y]+Cross[i,n[x,y,t]]

Print["----u,r,n,m,v,s,f,g"]
Print[u[x,y,t]]
Print[r[x,y,t]]
Print[n[x,y,t]]
Print[m[x,y,t]]
Print[v[x,y,t]]
Print[s[x,y,t]]
Print[f[x,y,t]]
Print[g[x,y,t]]

