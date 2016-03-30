pfix <- function(Sadj, Seq, beta, N){
	s = beta*Sadj-beta*Seq
	pfix = (1-exp(-2*s)) / (1-exp(-2*N*s))
	return(data.frame(s,pfix))
}

R0S = seq(1.2,3,length.out=1000)
	
smithConv = .07
distance = .6
birthRate = 1/(30*365)
deathRate = birthRate
gamma = 1/5
N = 15000000

meanStep = 0.6
sdStep = 0.4
alpha = (meanStep * meanStep) / (sdStep * sdStep)
beta = (sdStep * sdStep) / meanStep
x = seq(.0001,5,length=1000)
hx = dgamma(x,shape = alpha,rate = 1/beta)
hx = hx/sum(x*hx)

DISTANCES = c(0.1,0.3,0.6,1.0,2.0)

DISTANCES = c(0.6, 1.0, 2.0)
i=1
# for(distance in DISTANCES){
	# abline(v=distance,lty=i)
	# i=i+1
# }

pdf('R0_selection.pdf',3.425,4)

FIRST=TRUE
i=1
for(distance in DISTANCES){
	PFIXS = c()
	S = c()
	Scalc = c()
	for(R0 in R0S){
		beta = R0 * (gamma+birthRate)
		Seq = 1/R0
		Ieq = birthRate/beta*((R0)-1.0)
		Req = 1-Seq-Ieq
		
		Radj = Req*(1-smithConv*distance)
		Sadj = 1-Radj-1/N
		
		out = pfix(Sadj,Seq,beta,N)
	
		PFIXS = c(PFIXS, out$pfix)
		S = c(S, out$s)
		Scalc = c(Scalc,(R0-1)*(gamma*smithConv*distance+birthRate))
	}
	if(FIRST){
		plot(R0S,S,type='l',ylab=expression(paste('Selection coefficient ',italic('s'))),xlab=expression(italic('R'[0])),ylim=c(0,.06),lty=i)
		FIRST=FALSE
	}
	else{
		lines(R0S,S,type='l',lty=i)
	}
	i = i+1
}

dev.off()