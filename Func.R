library(MASS)
library(glmnet)
library(glasso)

BMCD.est = function(x, group, lam1, lam2, max_itr = 100)
{
	n = nrow(x)
	p = ncol(x)
	M = length(group)
	T = D.inverse = diag(p)
	for(j in c(M:2))
	{
		index = c(sum(group[1:(j-1)])+1):(sum(group[1:j]))
		y_temp1 = x[, index]							
		x_temp1 = x[, (1:sum(group[1:(j-1)]))]			

		p_j = ncol(y_temp1)							
		old_Omega = diag(p_j)
		old_A = matrix(0, nrow = p_j, ncol = ncol(x_temp1))
		diff_Omega = diff_A = 1
		maxiter = 1
		while((diff_Omega > 1e-3 | diff_A > 1e-3) & maxiter <= max_itr)
		{
			evalue = diag(eigen(old_Omega)$values)
			evector = eigen(old_Omega)$vectors
			halfOmega = evector %*% sqrt(evalue) %*% t(evector)
	
			y_temp2 = as.vector(y_temp1 %*% halfOmega)
			x_temp2 = kronecker(t(halfOmega), x_temp1)

			glmmod = glmnet(x_temp2, y_temp2, family = 'gaussian', alpha = 1)
			vec_A_transpose = as.vector(coef(glmmod, s = lam1))[-1]
			A_transpose = matrix(vec_A_transpose, ncol = p_j)
			new_A = t(A_transpose)

			S_temp = t(y_temp1 - x_temp1 %*% t(new_A)) %*% (y_temp1 - x_temp1 %*% t(new_A)) / n
			new_Omega = glasso(S_temp, lam2, penalize.diagonal=FALSE)$wi

			diff_Omega = sum((new_Omega - old_Omega)^2) / (p_j^2)
			diff_A = sum((new_A - old_A)^2) / (p_j * ncol(old_A))
			old_Omega = new_Omega
			old_A = new_A
			maxiter = maxiter + 1
		}
		T[index, 1:ncol(x_temp1)] = -old_A
		D.inverse[index,index] = old_Omega
	}
	index = 1:group[1]
	S_temp = t(x[, index]) %*% x[, index] / n
	D.inverse[index,index] = glasso(S_temp, lam2, penalize.diagonal=FALSE)$wi
	D.inverse[which(abs(D.inverse) < 1e-3)] = 0
	T[which(abs(T) < 1e-3)] = 0

	return(list(T = T, D.inverse = D.inverse))
}


BCD = function(x, group, lambda1 = NULL, lambda2 = NULL)
{
	if(is.null(lambda1))
	{
		lambda1 = c(0.05, 0.125, 0.2, 0.3, 0.5)
	}
	if(is.null(lambda2))
	{
		lambda2 = seq(0.3, 1.7, length.out = 6)
	}
	n = nrow(x)
	p = ncol(x)
	S = t(x) %*% x / n

	BIC = array(0, dim = c(length(lambda1), length(lambda2)))
	Omega_all = array(0, dim = c(p, p, length(lambda1), length(lambda2)))
	for(i in 1:length(lambda1))
	{
		for(k in 1:length(lambda2))
		{
			result = BMCD.est(x, group, lambda1[i], lambda2[k])
			Omega_all[,,i,k] = t(result$T) %*% result$D.inverse %*% result$T
			BIC[i,k] = -log(det(Omega_all[,, i, k])) + sum(diag(Omega_all[,, i, k] %*% S)) + log(n)*sum(abs(Omega_all[,, i, k]) > 1e-3)/n/2
		}
	}
	opt.row = which.min(BIC) %% nrow(BIC)
	opt.col = (which.min(BIC) %/% nrow(BIC)) + 1
	if(opt.row == 0)
	{
		opt.row = nrow(BIC)
		opt.col = which.min(BIC) %/% nrow(BIC)
	}
	Omega = Omega_all[,,opt.row, opt.col]
	Omega[which(abs(Omega) < 1e-3)] = 0
	return(Omega)	
}



AR = function(r, p)
{
	cov_design = array(0, dim = c(p,p))
	for(i in 1:p)
	{
		for(j in 1:i)
      	{
   			cov_design[i,j] = r^(abs(i-j))
			if(abs(cov_design[i,j]) < 1e-3)
			{
				cov_design[i,j] = 0
			}
			cov_design[j,i] = cov_design[i,j]
      	}
	}
	return(cov_design)
}

GAR = function(r, prow, pcol)
{
	A = matrix(0, nrow=prow, ncol=pcol)
	for(i in 1:min(prow, pcol))
	{
		for(j in 1:i)
      	{
   			A[i,j] = r^(abs(i-j))
			if(abs(A[i,j]) < 1e-3)
			{
				A[i,j] = 0
			}
			A[j,i] = A[i,j]
      	}
	}
	return(A)
}
