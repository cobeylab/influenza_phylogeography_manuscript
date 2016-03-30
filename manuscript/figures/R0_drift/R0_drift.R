#!/usr/bin/env Rscript


library(RSQLite)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(grid)
library(sensitivity)
library(scales)
library(xtable)
library(gtable)

resultsDir = '../../analysis/R0_drift'
resultsDb = paste(resultsDir,'results.sqlite',sep='/')
plotDir = './'

textSize = 6
pointSize = 1
lineSize = .8
plot_themes  = 	theme_classic() +
				theme(axis.line = element_line(size=.3)) +
				theme(axis.ticks = element_line(size=.3)) +
				theme(axis.ticks.length = unit(-.07,'cm')) +
				theme(axis.ticks.margin = unit(.11,'cm')) +
				theme(axis.title.x=element_text(size=textSize)) + 
				theme(axis.text.x=element_text(size= textSize)) + 
				theme(axis.title.y=element_text(size= textSize)) +
				theme(axis.text.y=element_text(size= textSize)) +
				theme(plot.title=element_text(size=textSize)) +
				theme(plot.margin=unit(c(3,3,3,3),'mm')) +
				theme(legend.title=element_text(size=textSize)) +
				theme(legend.text=element_text(size=textSize)) +
				theme(legend.position ='bottom') +
				theme(legend.direction='horizontal') +
				theme(legend.margin = unit(-.5,'cm')) +
				theme(panel.border = element_rect(colour = "black", fill=NA, size=.5)) +
				theme(axis.line = element_blank())


comboDb = dbConnect(SQLite(), dbname = resultsDb)
initExtension(comboDb)
parNames = c('R0')
plotDF = dbGetQuery(comboDb, paste('SELECT',paste(parNames,collapse=','),',meanFluxRate FROM pooled_results',sep=' '))

parName = 'R0'
fit = lm(meanFluxRate~R0,data=plotDF)

R2 = summary(fit)$r.squared
pearsonR = cor.test(plotDF$R0, plotDF$meanFluxRate)$estimate
CI = cor.test(plotDF$R0, plotDF$meanFluxRate)$conf.int

parName = 'R0'
betaLabel = paste('rho ==',round(rho,3))
pLabel = paste('p ==',round(p,3))
R2Label = paste('R^{2} ==',round(R2,3))
xLocation = (((max(plotDF[, parName])-min(plotDF[, parName]))/5+min(plotDF[, parName])))

R0Plot = ggplot(plotDF, aes_string(x=parName,y="meanFluxRate")) + 
	geom_point(size = pointSize*2) +
	xlab(expression(italic("R"[0]))) + 
	ylab("Mean antigenic drift rate\n(antigenic units per year)") + 
	ylim(c(0,2)) + 
	stat_smooth(method="lm",se=FALSE,size=lineSize) +
	geom_hline(aes(yintercept=1.01),linetype=2) +
	plot_themes +
				theme(axis.title.x=element_text(size= 10)) +
				theme(axis.text.x=element_text(size= 10)) + 
				theme(axis.title.y=element_text(size= 10)) +
				theme(axis.text.y=element_text(size= 10)) +
				theme(plot.title=element_text(size=10)) +
				theme(plot.margin=unit(c(2,2,2,2),'mm')) 
pdf('R0_drift.pdf',width=3.425, height = 3.425)
print(R0Plot)
dev.off()
