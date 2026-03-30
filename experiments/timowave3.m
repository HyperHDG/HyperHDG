i:= {0, 1, 0}

u[x_,y_,z_,t_] := {1t y^3, 2t y^3, 3t y^3}
r[x_,y_,z_,t_] := {5t y^3, 7t y^3, 11 t y^3}

v[x_,y_,z_,t_] =  D[u[x,y,z,t],t]
s[x_,y_,z_,t_] =  D[r[x,y,z,t],t]
n[x_,y_,z_,t_] = -D[u[x,y,z,t],y]-Cross[i,r[x,y,z,t]]
m[x_,y_,z_,t_] = -D[r[x,y,z,t],y]
f[x_,y_,z_,t_] =  D[v[x,y,z,t],t]+D[n[x,y,z,t],y]
g[x_,y_,z_,t_] =  D[s[x,y,z,t],t]+D[m[x,y,z,t],y]+Cross[i,n[x,y,z,t]]

Print["----u,r,n,m,v,s,f,g"]
Print[u[x,y,z,t]]
Print[r[x,y,z,t]]
Print[n[x,y,z,t]]
Print[m[x,y,z,t]]
Print[v[x,y,z,t]]
Print[s[x,y,z,t]]
Print[f[x,y,z,t]]
Print[g[x,y,z,t]]
