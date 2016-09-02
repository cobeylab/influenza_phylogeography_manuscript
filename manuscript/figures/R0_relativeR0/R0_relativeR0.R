#!/usr/bin/env Rscript
library(RSQLite)
library(ggplot2)
library(reshape2)
library(grid)
library(gridExtra)
library(scales)
library(plyr)

resultsDir = '../../analysis/R0_relativeR0/'
plotDir = './'

figureWidth = 4.5 #(inches)
aspectRatio = 0.6

textSize = 9
pointSize = 1.0
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
				theme(plot.title=element_text(size=textSize)) +
				theme(plot.title=element_text(size=textSize, hjust=-0.2)) +
				theme(plot.margin=unit(c(1,1,1,1),'mm')) +
				theme(legend.title=element_text(size=textSize)) +
				theme(legend.text=element_text(size=textSize)) +
				theme(legend.position ='bottom') +
				theme(legend.direction='horizontal') +
				theme(legend.margin = unit(0,'cm')) +
				theme(panel.border = element_rect(colour = "black", fill=NA, size=.5)) +
				theme(axis.line = element_blank()) +
				theme(text = element_text(family='serif')) 

resultsDb = paste(resultsDir,'results.sqlite',sep='')

parNames = c('relativeR0','R0')

comboDb = dbConnect(SQLite(), dbname = resultsDb)
initExtension(comboDb)


antigenDF = dbGetQuery(comboDb, paste('SELECT',paste(parNames,collapse=','),', runId, trunkProportions, antigenicLagAG1, extinct, excessDiversity, fluLike FROM pooled_results',sep=' '))

plotDF = ddply(antigenDF, .(relativeR0, R0), summarise, tropicsAgLag = mean(antigenicLagAg1,na.rm=TRUE), tropicsTrunkPro = mean(trunkProportions, na.rm=TRUE))
plotDF[is.na(plotDF)] = -1

parName1 = parNames[1]
parName2 = parNames[2]
trunkProPlot2D = ggplot(plotDF, aes_string(x=parName1,y=parName2)) + 
	geom_tile(aes(fill=tropicsTrunkPro)) + 
	scale_fill_gradientn("tropics trunk proportion",
					colours = c('#2A4A7F',"white","#800000"),
					values = rescale(c(0,1/3,1)),
					limits = c(0,1)) + 
	guides(fill=guide_colorbar(barwidth=8,barheight=0.4, title.position = 'top')) +
	xlab(expression(paste('relative ',italic('R'[0])))) + 
	ylab(expression(paste('baseline ',italic('R'[0])))) + 
	plot_themes + 
	ggtitle(expression(paste("(",italic(a),")"))) + 
	scale_x_continuous(expand = c(0,0)) +
	scale_y_continuous(expand = c(0,0)) +
	theme(panel.border = element_rect(colour = "black", fill=NA, size=.2)) +
	theme(axis.ticks = element_line(size=.1)) 


agLagPlot2D = ggplot(plotDF, aes_string(x=parName1,y=parName2)) + 
	geom_tile(aes(fill=tropicsAgLag)) + 
	scale_fill_gradientn("tropics antigenic lead",
						colours = c('#2A4A7F',"white","#800000"),
						values = rescale(c(-.4,0,.4)),
						limits = c(-.25,.25)) + 
	guides(fill=guide_colorbar(barwidth=8,barheight=0.4, title.position = 'top')) +
	xlab(expression(paste('relative ',italic('R'[0])))) + 
	ylab(expression(paste('baseline ',italic('R'[0])))) + 
	plot_themes + 
	ggtitle(expression(paste("(",italic(b),")"))) + 
	scale_x_continuous(expand = c(0,0)) +
	scale_y_continuous(expand = c(0,0)) +
	theme(panel.border = element_rect(colour = "black", fill=NA, size=.2)) +
	theme(axis.ticks = element_line(size=.15)) 

pdf(paste(plotDir,'R0_relativeR0.pdf',sep='/'),figureWidth, figureWidth*aspectRatio)
plots = list(trunkProPlot2D, agLagPlot2D)
args_list = c(plots,1,2)
names(args_list)=c(letters[1:length(plots)],'nrow','ncol')
do.call(grid.arrange,args_list)
dev.off()

dbDisconnect(comboDb)