beta = .36
nu = 0.2
N = 45000000
gamma = 9.1e-05

R0 = beta/(gamma+nu)
beta2 = beta*1.2
R0_2 = beta2/(gamma+nu)

S = 1/R0
I = gamma/beta*(R0-1)

S2 = 1/R0_2
I2 = gamma/beta2*(R0_2-1)

lambda = beta*I
lambda2 = beta2*I2
t=seq(0,20*365) #days

p1 = 1-exp(-lambda * t)
p2 = 1-exp(-lambda2* t)
pdf('SIR_seropositivity.pdf',width=3.425, height =4)
plot(t/365, p1,type='l', xlab = 'Age (years)',ylab='Fraction seropositive',ylim=c(0,1))
lines(t/365,p2,type='l',col='red')
lines(t/365,p2-p1,col='blue')
legend(x=0,y=1,legend = c(expression(paste(italic(R)[0],' = 1.8')),expression(paste(italic(R)[0],' = 2.16')),'difference'),col=c('black','red','blue'),lty=c(1,1,1))
dev.off()
alpha = .05
beta = .2

f.alpha.beta = (qnorm(1-alpha/2) + qnorm(1-beta))^2
age = 365*2
N = (p1[age]*(1-p1[age]) + p2[age]*(1-p2[age]))/(p1[age]-p2[age])^2 * f.alpha.beta
print(N)

f1 = .05
f2 = .03
N = (f1*(1-f1) + f2*(1-f2))/(f1-f2)^2 * f.alpha.beta

p1[365*2]
p2[365*2]

#t=seq(0,8*pi,length=100)
# par(mfrow=c(3,1))
# # plot(t,sin(t),type='l', xlab='',ylab='', xaxt='n', yaxt='n',bty='n',lwd=4,col=rgb(0.765, 0.728, 0.274))
# # plot(t,rep(1,length(t)),type='l', xlab='',ylab='', xaxt='n', yaxt='n',bty='n',lwd=4, col = rgb(0.324, 0.609, 0.708))
# # plot(t,sin(t+pi),type='l', xlab='',ylab='', xaxt='n', yaxt='n',bty='n',lwd=4, col=rgb(0.857, 0.131, 0.132))
