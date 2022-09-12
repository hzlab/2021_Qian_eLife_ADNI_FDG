################################################################
# myInput
# This function creates example used
#
#		  n: number of observations
#		  p: number of predictors
#		  T: number of time points sampled
#             sigma: error standard deviation
#      do.intercept: whether intercept should be included
#               rho: correlation of design matrix
#              zeta: similarity between coefficient fcns
################################################################
# Coefficient functions used; these are the examples in Daye et al. 2012 p11-12, not needed for any actual data analysis
# we need the sparse function from line 125 onwards only

bfcn_pc <- function(tVec) { return((-1*(tVec<.5) + (tVec>=.5))) }
bfcn_pn1 <- function(tVec) { return(2*tVec-1) }
bfcn_pn2 <- function(tVec) { return((1-2*tVec)^3) }
bfcn_lf1 <- function(tVec) { return(sin(2*pi*tVec)) }
bfcn_lf2 <- function(tVec) { return(cos(2*pi*tVec)) }
bfcn_hf1 <- function(tVec) { return(sin(8*pi*tVec)) }
bfcn_hf2 <- function(tVec) { return(cos(8*pi*tVec)) }
bfcn_cp <- function(tVec) { return( ((1-2*tVec)^3-cos(4*pi*tVec)) ) }	

myInput <- function(n, p, T, sigma, do.intercept, rho=0.5, zeta=NULL){
	tVec = seq(0,1,len=T)

	ygt = array(0,dim=c(n,T))
	X = array(0, dim=c(n,p,T));
	beta = array(0, dim=c(p,T));
	beta0 = rep(0,T)
	trueNNZ = rep(FALSE,p+do.intercept);

	# Example setup
	{
		if (do.intercept) {
	                beta0 = 2*bfcn_cp(tVec)
		}

		beta[1,] = (1-zeta)*bfcn_lf2(tVec) + zeta*bfcn_pn1(tVec)
		beta[2,] = (1-zeta)*bfcn_lf2(tVec) + zeta*bfcn_hf1(tVec)
		beta[3,] = (1-zeta)*bfcn_lf2(tVec) + zeta*bfcn_lf1(tVec)
		beta[4,] = (1-zeta)*bfcn_lf2(tVec) + zeta*bfcn_hf2(tVec)

		beta[5,] = (1-zeta)*bfcn_pn2(tVec) + zeta*bfcn_pn1(tVec)
		beta[6,] = (1-zeta)*bfcn_pn2(tVec) + zeta*bfcn_hf1(tVec)
		beta[7,] = (1-zeta)*bfcn_pn2(tVec) + zeta*bfcn_lf1(tVec)
		beta[8,] = (1-zeta)*bfcn_pn2(tVec) + zeta*bfcn_hf2(tVec)

		beta[9,] = 3*bfcn_hf1(tVec)
		beta[10,] = 1.5*bfcn_pn2(tVec)
		beta[11,] = 2*bfcn_lf1(tVec)
		beta[12,] = 2*bfcn_pc(tVec)

		Grp1 = 1:4
		Grp2 = 5:8

		X = array(0, dim=c(n,p,T))
		for (t in 1:T) {
			Mu = rep(0,p)
			Sigma = diag(p);
			for (i in 1:p) {
                		for (j in 1:p) {
                        		if (i!=j) {
                                		Sigma[i,j] = 0^abs(i-j);
                        		}
                		}
        		}

			for (i in Grp1) {
				for (j in Grp1) {
					if (i!=j) {Sigma[i,j]=rho}
				}
			}

			for (i in Grp2) {
				for (j in Grp2) {
					if (i!=j) {Sigma[i,j]=rho}
				}
			}

			X[,,t] = mvrnorm(n, Mu, Sigma);
		}

		# Sets weights and basis inner products
		W = matrix(0,p,p)

		for (i in Grp1) {
			for (j in Grp1) {
				if (i!=j) {W[i,j]=1}
			}
		}

		for (i in Grp2) {
			for (j in Grp2) {
				if (i!=j) {W[i,j]=1}
			}
		}

	}

	# Creates response
	for (t in 1:T) {
		ygt[,t] = ygt[,t] + beta0[t];
		if (any(beta0!=0)) {trueNNZ[1] = TRUE}
		for (j in 1:p) {
			ygt[,t] = ygt[,t] + beta[j,t]*X[,j,t] + sigma*rnorm(n)
			if (any(beta[j,]!=0)) {
				trueNNZ[j+do.intercept] = TRUE;
			}
		}
	}

	if (!do.intercept) {
		ygt = ygt - repmat(colMeans(ygt),nrow(ygt),1);
	}

	ret.obj = list("tVec"=tVec, "X" = X, "ygt" = ygt, "W"=W)

	return(ret.obj)

}


################################################################
# vc_sparse
# Performs sparse penalty varying-coefficient model
#
#            simLen: number of repetitions for the simulation study
#		  n: number of observations(subjects)
#		  T: number of time points sampled
#		  p: number of predictors
#             sigma: error standard deviation
#  	  do.alpha2: whether lambda_2 (structural penalty) should be used
#  	  do.alpha3: whether c (rescaling parameter) should be used
#               rho: correlation of design matrix (used for simulation)
#              zeta: similarity between coefficient fcns (used for simulation)
################################################################
vc_sparse <- function(simLen,n, T, p, 
	              X,ygt,W,tVec,
	              do.alpha2=FALSE, do.alpha3=TRUE, rootDir =''
            
                      )
{
# Loading packages
require(MASS)
require(boot)
require(splines)
require(statmod)
source('/mnt/isilon/CSC4/HelenZhouLab/HZLHD2/Data7/Members/Qianxing/ADNI_NgKP/svc/vc/misc.R')
source('/mnt/isilon/CSC4/HelenZhouLab/HZLHD2/Data7/Members/Qianxing/ADNI_NgKP/svc/vc/bsplines.R')
source('/mnt/isilon/CSC4/HelenZhouLab/HZLHD2/Data7/Members/Qianxing/ADNI_NgKP/svc/vc/vc.R') # calls vccd?
#system("rm -f *.so"); system("rm -f *.o"); system("rm -f *.mod"); 
#system("R CMD SHLIB vccd_module.f90"); system("R CMD SHLIB vccd.f90 vccd_module.f90 minpack/*.f")
dyn.load('/mnt/isilon/CSC4/HelenZhouLab/HZLHD2/Data7/Members/Qianxing/ADNI_NgKP/svc/vc/vccd.so')



# Default parameters
Tol = 1e-6
cvFold = 5
L = 4
do.intercept = FALSE
do.standardize = FALSE
degree=3
no.interior.knots = L-degree-1
nVarsLimit=50
maxIters=100
maxInnerIters=5

#cat("====================================================\n")

# Tuning parameters used
alpha_1 = 2^(seq(log2(1+Tol), log2(0+Tol), len=25))-Tol
if (do.alpha2) {
	alpha_2 = 2^(seq(log2(1.5+Tol), log2(0+Tol), len=25))-Tol #seq(1,0,len=100)
} else {
	alpha_2 = c(0)
}
if (do.alpha3) {
	alpha_3 = seq(0,1,len=25)
} else {
	alpha_3 = c(0)
}

VC_gamma = array(0,c(simLen,L*(p+do.intercept)))
VC_I = rep(0,simLen)
VC_C = rep(0,simLen)
VC_sens = rep(0,simLen)
VC_spec = rep(0,simLen)
VC_G = rep(0,simLen)
VC_NNZ = array(0,c(simLen,p+do.intercept))
VC_bias = array(0,c(simLen,p+do.intercept))
VC_var = array(0,c(simLen,p+do.intercept))

for (simI in 1:simLen) {
        #clock_tmp2 = tic()
	# Declares real data used 
	#ex.obj = myInput(n, p, T, sigma, do.intercept, rho, zeta)
	#tVec = ex.obj$tVec
	#sigma = ex.obj$sigma
	#trueNNZ = ex.obj$trueNNZ
	#Sigma = ex.obj$Sigma

	#X = ex.obj$X
	#ygt = ex.obj$ygt
	

	#ex.test.obj = myInput(500, p, T, sigma, do.intercept, rho, zeta)

	# Used in computing bias with T=500 time points.
	#ex.emp.obj = myInput(1, p, 500, sigma, do.intercept, rho, zeta)
	#betaEmp = ex.emp.obj$trueBeta
	#tVecEmp = ex.emp.obj$tVec

	# Performs cross-validation
	cv.obj = cv.vc(tVec, X, ygt, W, K = cvFold,
        	alpha_1=alpha_1, alpha_2=alpha_2, alpha_3=alpha_3,
        	degree=degree, no.interior.knots=no.interior.knots,
        	do.intercept=do.intercept, do.standardize=do.standardize,
        	nVarsLimit=nVarsLimit, maxIters=maxIters, maxInnerIters=maxInnerIters)

	if (length(alpha_2)==1) {
		if (length(alpha_3)==1) {
			cv.obj$mse_mean = array(cv.obj$mse_mean,dim(cv.obj$mse_mean)[c(1)])
		} else {
			cv.obj$mse_mean = array(cv.obj$mse_mean,dim(cv.obj$mse_mean)[c(1,3)])
		}
	} else {
		if (length(alpha_1)==1) {
			cv.obj$mse_mean = array(cv.obj$mse_mean,dim(cv.obj$mse_mean)[c(1,2)])
		}
	}
	gamma = cv.obj$gamma
	VC_gamma[simI,] = gamma;


        # obtain beta results
	beta_res = get_beta(tVec, gamma, p+do.intercept, 
		degree=degree, no.interior.knots=no.interior.knots)

	# Saves performance statistics
	#betaTmp = get_beta(tVecEmp, gamma, p+do.intercept, 
		#degree=degree, no.interior.knots=no.interior.knots)
	#VC_bias[simI,] = rowMeans((betaTmp - betaEmp)^2)
	#VC_var[simI,] = rowMeans((betaTmp - rowMeans(betaTmp))^2)

	#tmpGamma = gamma; dim(tmpGamma) = c(L,p+do.intercept); tmpGamma = t(tmpGamma);
	#tmpNNZ = rowSums(tmpGamma!=0)>0
	#VC_NNZ[simI,] = tmpNNZ
	#VC_C[simI] = sum(tmpNNZ[trueNNZ==TRUE])
	#VC_I[simI] = sum(tmpNNZ[trueNNZ!=TRUE])
	#VC_sens[simI] = VC_C[simI]/sum(trueNNZ)
	#VC_spec[simI] = (p+do.intercept-sum(trueNNZ) - VC_I[simI])/
	#			(p+do.intercept-sum(trueNNZ))
	#VC_G[simI] = sqrt(VC_sens[simI]*VC_spec[simI])

	# Display simI count
	#cat(simI,'  ')
	if (simI==simLen){ cat('\n') }
}

## Outputs performance results

# Sensitivity (bootsd)
#cat( sprintf("Sens \t\t%.3f (%.3f)\n", median(VC_sens), bootsd(VC_sens)) )
# Specificity (bootsd)
#cat( sprintf("Spec \t\t%.3f (%.3f)\n", median(VC_spec), bootsd(VC_spec)) )
# G-measure (bootsd)
#cat( sprintf("G-measure \t%.3f (%.3f)\n", median(VC_G), bootsd(VC_G)) )
# True positives (bootsd)
#cat( sprintf("TP \t\t%.1f (%.1f)\n", median(VC_C), bootsd(VC_C)) )
# False positives (bootsd)
#cat( sprintf("FP \t\t%.1f (%.1f)\n", median(VC_I), bootsd(VC_I)) )

cat("\n")

#VC_bias_sd = rep(0,sum(trueNNZ==TRUE))
#VC_var_sd = rep(0,sum(trueNNZ==TRUE))
#idx = 1
#for (j in which(trueNNZ==TRUE)) {
#	VC_bias_sd[idx] = bootsd(VC_bias[,j])
#	VC_var_sd[idx] = bootsd(VC_var[,j])
#	idx = idx + 1
#}

# Percent selected of true variables
#cat("\t\tbeta0 \tbeta1 \tbeta2 \tbeta3 \tbeta4 \tbeta5 \tbeta6 \tbeta7 \tbeta8 \tbeta9 \tbeta10 \tbeta11 \tbeta12 \n")
#cat("% Selected \t", sprintf("%.3f\t", colMeans(VC_NNZ[,trueNNZ==TRUE])), "\n", sep="")
# Bias squared of true variables (bootsd)
#cat("Bias^2 \t\t", sprintf("%.3f\t", colMeans(VC_bias[,trueNNZ==TRUE])), "\n", sep="")
#cat("(bootsd)\t", sprintf("%.3f\t", VC_bias_sd), "\n", sep="")
# Variance of true variables (bootsd)
#cat("Var \t\t", sprintf("%.3f\t", colMeans(VC_var[,trueNNZ==TRUE])), "\n", sep="")
#cat("(bootsd) \t", sprintf("%.3f\t", VC_var_sd), "\n", sep="")


ret.obj=list("beta"=beta_res)
return(ret.obj)

}

