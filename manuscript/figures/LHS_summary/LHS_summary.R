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

resultsDir = '../../analysis/LHS_analysis'
plotDir = './'

figureWidth = 7 
textSize = 6
pointSize = 0.05
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

resultsDb = paste(resultsDir,'results.sqlite',sep='/')

parNames = c('relativeR0','seasonalAmplitude','relativeN','relativeTurnover','tropicFractionI0')
headerNames = c('\"Relative R\"[0]','\"Seasonal amplitude\"','\"Relative population size\"','\"Relative turnover\"','\"Fraction of initial infecteds in tropics\"')
headers = data.frame(parNames,headerNames)
comboDb = dbConnect(SQLite(), dbname = resultsDb)
initExtension(comboDb)

plotDF = dbGetQuery(comboDb, paste('SELECT',paste(parNames,collapse=','),', tropicsTrunkPro, tropicsAgLag FROM pooled_results',sep=' '))
tTPsens = pcc(plotDF[,1:5], plotDF$tropicsTrunkPro,nboot=100)
tALsens = pcc(plotDF[,1:5], plotDF$tropicsAgLag,nboot=100)
agAnova = aov(tropicsAgLag~relativeN+ tropicFractionI0 +relativeR0+relativeTurnover+seasonalAmplitude,data=plotDF)
trunkAnova = aov(tropicsTrunkPro~relativeN+ tropicFractionI0 +relativeR0+relativeTurnover+seasonalAmplitude,data=plotDF)

plots=c()

for(parName in parNames){
	print(parName)
	sens = tTPsens$PCC[parName,]$'original' 
	minCI = tTPsens$PCC[parName,]$'min. c.i.'
	maxCI = tTPsens$PCC[parName,]$'max. c.i.'

	trunkProPlot = ggplot(plotDF, aes_string(x=parName,y="tropicsTrunkPro")) + 
		geom_point(size = pointSize, alpha = 0.2) + 
		xlab(NULL) +
		ylab(NULL) +
		ylim(c(0,1)) + 
		theme_classic() + 
		ggtitle(paste('r = ',round(sens,2), '\n95% CI: (',round(minCI,2),', ',round(maxCI,2),')',sep='')) +
		plot_themes +
		theme(axis.text.x=element_blank()) 
	
	if(parName == parNames[1]){
		trunkProPlot = trunkProPlot + ylab("Tropics Trunk Proportion")
		trunkProPlotRow = ggplotGrob(trunkProPlot)
	}
	else{
		trunkProPlot = trunkProPlot + theme(axis.text.y=element_blank()) 
		trunkProPlotRow = cbind(trunkProPlotRow,ggplotGrob(trunkProPlot),size='first')
	}
}	

for(parName in parNames){
	sens = round(tALsens$PCC[parName,]$'original' ,2)
	minCI = round(tALsens$PCC[parName,]$'min. c.i.',2)
	maxCI = round(tALsens$PCC[parName,]$'max. c.i.',2)

	print(parName)		
	agLagPlot = ggplot(plotDF, aes_string(x=parName,y="tropicsAgLag")) + 
		geom_point(size = pointSize, alpha = 0.2) + 
		xlab(parse(text=paste(headers$headerNames[headers$parNames==parName]))) + 
		ylab(NULL) + 
		ylim(c(-0.4,0.4)) + 
		theme_classic() +
		ggtitle(paste('r = ',round(sens,2), '\n95% CI: (',round(minCI,2),', ',round(maxCI,2),')',sep='')) +
		plot_themes
		
	if(parName == parNames[1]){
		agLagPlot = agLagPlot + ylab("Tropics Antigenic Lead")
		agLagPlotRow = ggplotGrob(agLagPlot)
	}
	else{
		agLagPlot = agLagPlot + theme(axis.text.y=element_blank()) 
		agLagPlotRow = cbind(agLagPlotRow, ggplotGrob(agLagPlot), size='first')
	}
}

figureS2 = rbind(trunkProPlotRow,agLagPlotRow,size='first')
pdf(paste(plotDir,'LHS_summary.pdf',sep=''),figureWidth,figureWidth/2)
grid.draw(figureS2)
dev.off()
	
dbDisconnect(comboDb)

print(xtable(summary(agAnova)))
print(xtable(summary(trunkAnova)))