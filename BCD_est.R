source("Func.r")

n = 50							
p = 20		

## group 1
group = rep(4,5)
Omega_design = AR(0.8, p)	


## group 2
group = c(30, 60, 40, 70)	
Omega_design = GAR(0.8, p, p)


x = mvrnorm(n, rep(0,p), solve(Omega_design))                             
x = scale(x, scale = FALSE)	

Omega_hat = BCD(x, group)
Omega_hat





