#! /usr/bin/Rscript --vanilla

args <- commandArgs(TRUE)

if(!length(args) == 1) {
    message("Need one argument in the format: %.csv")
    q(status=1)
}

dat <- read.csv(args[1])

library(ggplot2)

totalArea <- function(dat) {
    return(  ggplot(dat, aes(AcquiredTime, TotalArea))
           + facet_grid(PrecursorMz+PeptideModifiedSequence~., scale="free")
           + stat_smooth(aes(group=1), method="lm", formula=y~1, level=0.99, colour="black", linetype=4, fill="yellow", alpha=.2)
           + stat_smooth(aes(group=1), method="lm", formula=y~1, level=0.95, colour="black", linetype=4, fill="blue", alpha=.3)
           + geom_point()
           + ggtitle(expression(atop(  "Total area for different runs"
                                     , "together with 95% and 99% confidence intervals")))
           + theme(axis.text.x = element_text(angle = 90, hjust = 1))
           )
}

isotopeDotProduct <- function(dat) {
    return(  ggplot(dat, aes(AcquiredTime, IsotopeDotProduct))
           + facet_grid(PrecursorMz+PeptideModifiedSequence~., scale="free")
           + stat_smooth(aes(group=1), method="lm", formula=y~1, level=0.99, colour="black", linetype=4, fill="yellow", alpha=.2)
           + stat_smooth(aes(group=1), method="lm", formula=y~1, level=0.95, colour="black", linetype=4, fill="blue", alpha=.3)
           + geom_point()
           + ggtitle(expression(atop(  "Isotope dot product for different runs"
                                     , "together with 95% and 99% confidence intervals")))
           + theme(axis.text.x = element_text(angle = 90, hjust = 1))
           )
}

bestRetentionTime <- function(dat) {
    return(  ggplot(dat, aes(AcquiredTime, BestRetentionTime))
           + facet_grid(PrecursorMz+PeptideModifiedSequence~., scale="free")
           + stat_smooth(aes(group=1), method="lm", formula=y~1, level=0.99, colour="black", linetype=4, fill="yellow", alpha=.2)
           + stat_smooth(aes(group=1), method="lm", formula=y~1, level=0.95, colour="black", linetype=4, fill="blue", alpha=.3)
           + geom_point()
           + ggtitle(expression(atop(  "Retention time for different runs"
                                     , "together with 95% and 99% confidence intervals")))
           + theme(axis.text.x = element_text(angle = 90, hjust = 1))
           )
}

maxFwhm <- function(dat) {
    return(  ggplot(dat, aes(AcquiredTime, MaxFwhm))
           + facet_grid(PrecursorMz+PeptideModifiedSequence~., scale="free")
           + stat_smooth(aes(group=1), method="lm", formula=y~1, level=0.99, colour="black", linetype=4, fill="yellow", alpha=.2)
           + stat_smooth(aes(group=1), method="lm", formula=y~1, level=0.95, colour="black", linetype=4, fill="blue", alpha=.3)
           + geom_point()
           + ggtitle(expression(atop(  "Width at half maximum intensity for different runs"
                                     , "together with 95% and 99% confidence intervals")))
           + theme(axis.text.x = element_text(angle = 90, hjust = 1))
           )
}

averageMassErrorPPM <- function(dat) {
    return(  ggplot(dat, aes(AcquiredTime, AverageMassErrorPPM))
           + facet_grid(PrecursorMz+PeptideModifiedSequence~., scale="free")
           + stat_smooth(aes(group=1), method="lm", formula=y~1, level=0.99, colour="black", linetype=4, fill="yellow", alpha=.2)
           + stat_smooth(aes(group=1), method="lm", formula=y~1, level=0.95, colour="black", linetype=4, fill="blue", alpha=.3)
           + geom_point()
           + ggtitle(expression(atop(  "Average mass error in ppm for different runs"
                                     , "together with 95% and 99% confidence intervals")))
           + theme(axis.text.x = element_text(angle = 90, hjust = 1))
           )
}

pdf(gsub('csv', 'pdf', args[1]), height=11.6, width=8.2)
print(totalArea(dat))
print(isotopeDotProduct(dat))
print(bestRetentionTime(dat))
print(maxFwhm(dat))
print(averageMassErrorPPM(dat))
dev.off()
