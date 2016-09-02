library(ggplot2)
library(RSQLite)
library(reshape2)
library(gridExtra)
library(grid)
library(scales)

figureWidth = 4.5 #(inches)
aspectRatio = .45

resultsDir = '../../analysis/univariate_analyses'
plotDir = './'

textSize = 9
pointSize = .5
lineSize = .3
plot_themes  = 	theme_classic() +
				theme(axis.line = element_line(size=.3)) +
				theme(axis.ticks = element_line(size=.3)) +
				theme(axis.ticks.length = unit(-.07,'cm')) +
				theme(axis.ticks.margin = unit(.11,'cm')) +
				theme(axis.title.x=element_text(size=textSize)) + 
				theme(axis.text.x=element_text(size= textSize)) + 
				theme(axis.title.y=element_text(size= textSize)) +
				theme(axis.text.y=element_text(size= textSize)) +
				theme(plot.title=element_text(size=textSize, hjust=-0.25)) +
				theme(plot.margin=unit(c(1,1,1,1),'mm')) +
				theme(legend.title=element_text(size=textSize)) +
				theme(legend.text=element_text(size=textSize)) +
				theme(legend.position ='bottom') +
				theme(legend.direction='horizontal') +
				theme(legend.margin = unit(-.5,'cm')) +
				theme(panel.border = element_rect(colour = "black", fill=NA, size=.5)) +
				theme(axis.line = element_blank()) +
				theme(text = element_text(family='serif')) 


resultsDb = paste(resultsDir,'results.sqlite',sep='/')
comboDb = dbConnect(SQLite(), dbname = resultsDb)
initExtension(comboDb)

parName = 'relativeR0'
PARNAMES = c('relativeN','tropicFractionI0','relativeR0','relativeTurnover','seasonalAmplitude')
antigenDF = dbGetQuery(comboDb, paste('SELECT',paste(PARNAMES,collapse=','),', sweepID,tropicsTrunkPro, tropicsAgLag FROM pooled_results',sep=' '))
plotDF = antigenDF[which(antigenDF$sweepID==parName),]

agLead.corr = cor.test(plotDF[,parName],plotDF$tropicsAgLag)
agLead.fit = lm(tropicsAgLag~plotDF[,parName], data=plotDF)
#Get 95% confidence interval for estimated increase in relative R0 for observed lead
predict.agLead.fit = lm(plotDF[,parName]~ tropicsAgLag, data=plotDF)
newdat = data.frame(tropicsAgLag = 0.25)
agLead.CI = predict(predict.agLead.fit, newdat, interval='confidence')

agLeadPlot = ggplot(plotDF, aes_string(x=parName,y="tropicsAgLag")) + 
	geom_point(size = pointSize) + 
	xlab(expression(paste("relative ",italic("R"[0]),sep=" "))) + 
	ylab("tropics antigenic lead") + 
	ylim(c(-0.4,0.4)) +
	stat_smooth(method="lm",se=FALSE,size=lineSize) +
	geom_hline(aes(yintercept=0.0),linetype=2, size = lineSize) + 
	plot_themes + 
	ggtitle(expression(paste("(",italic(b),")"))) 

trunkPro.corr = cor.test(plotDF[, parName],plotDF$tropicsTrunkPro)
trunkPro.fit = lm(tropicsTrunkPro~plotDF[,parName], data=plotDF)
#Get 95% confidence interval for estimated increase in relative R0 for observed trunk proportion
predict.trunkPro.fit = lm(plotDF[,parName]~ tropicsTrunkPro, data=plotDF)
newdat = data.frame(tropicsTrunkPro = 0.87)
trunkPro.CI = predict(predict.trunkPro.fit, newdat, interval='confidence')

trunkProPlot = ggplot(plotDF, aes_string(x=parName,y="tropicsTrunkPro")) + 
	geom_point(size = pointSize) + 
	xlab(expression(paste("relative ",italic("R"[0]),sep=" "))) + 
	ylab("tropics trunk proportion") + 
	ylim(c(0,1)) + 
	stat_smooth(method="lm",se=FALSE,size=lineSize) +
	geom_hline(aes(yintercept=1/3),linetype=2, size = lineSize) +
	plot_themes +
	ggtitle(expression(paste("(",italic(a),")")))

dbDisconnect(comboDb)

plots = list(trunkProPlot,agLeadPlot)
args_list = c(plots,1,2)
names(args_list)=c(letters[1:length(plots)],'nrow','ncol')

pdf(paste(plotDir,'relativeR0.pdf',sep='/'), width = figureWidth, height = figureWidth * aspectRatio)
do.call(grid.arrange,args_list)
dev.off()

