#! /usr/bin/Rscript --vanilla

# Copyright (C) 2014 Thomas W. D. Möbius (kontakt@thomasmoebius.de)
#
#     This program is free software: you can redistribute it and/or modify it
#     under the terms of the GNU General Public License as published by the
#     Free Software Foundation, either version 3 of the License, or (at your
#     option) any later version.
#
#     This program is distributed in the hope that it will be useful, but
#     WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
#     Public License for more details.
#
#     You should have received a copy of the GNU General Public License along
#     with this program. If not, see <http://www.gnu.org/licenses/>.

arg <- commandArgs(TRUE)

fileCSV <- grep('.csv$', arg, value=TRUE)
verbose <- is.element("verbose", arg)
pdfFileName <- gsub('csv', 'pdf', fileCSV)

if (verbose) {
    message("Command line arguments ", arg)
    message("  Report file: ", fileCSV)
    message("  Verbose: ", verbose)
    message("  Output file will be ", pdfFileName)
}

if(!length(fileCSV) == 1) {
    message("Need one single report file in the format: %.csv")
    q(status=1)
}

dat <- read.csv(fileCSV, na.strings="#N/A")
dat$ReplicateName <- NULL #factor(make.names(dat$ReplicateName))
dat$PeptideModifiedSequence <- factor(gsub("\\]|\\[|\\+", "", dat$PeptideModifiedSequence))
dat$AcquiredTime <- gsub("20([0-9][0-9])", "\\1", dat$AcquiredTime)
dat$AcquiredTime <- as.Date(dat$AcquiredTime, format="%m/%d/%y")
dat$AcquiredTime <- factor(dat$AcquiredTime, ordered=T)

if (verbose) {
    message("Successfully parsed report file having the following fields:")
    print(names(dat))
}

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(plyr))

if (verbose) message("Successfully loaded packages: ggplot2, plyr")

pseudoNames <- function (n, nump) {
    if (n < nump) {
        rep(1, n%%nump)
    } else {
        c(rep(1:(n%/%nump), each=nump), rep((n%/%nump)+1, n%%nump))
    }
}

createPseudoProteinNames <- function(oldnames, nump) {
    return(mapvalues(oldnames,
                     levels(oldnames),
                     pseudoNames(nlevels(oldnames), nump)))
}

nump = 10
if (is.null(dat$ProteinName))
    dat$ProteinName <- createPseudoProteinNames(dat$PeptideModifiedSequence, nump)

qcplot <- function(dat, qcstatistic) {
    p <-(  ggplot(dat, aes(x=AcquiredTime))
         + geom_point(aes_string(y=qcstatistic))
         + facet_grid(PrecursorMz+PeptideModifiedSequence~., scale="free")
         + theme(axis.text.x = element_text(angle = 90, hjust = 1))
         + geom_hline(aes(yintercept=qc.mean, colour="black"))
         + geom_hline(aes(yintercept=qc.lwr,  colour="darkgreen"))
         + geom_hline(aes(yintercept=qc.upr,  colour="darkgreen"))
         + geom_hline(aes(yintercept=qc.lwr2, colour="blue"))
         + geom_hline(aes(yintercept=qc.upr2, colour="blue")))
    return(p)
}

qcstat.TotalArea <- function (dat) {
    ds <- ddply(dat, .(PeptideModifiedSequence, PrecursorMz), summarise,
                qc.mean = mean(TotalArea, na.rm=T), qc.sd = sd(TotalArea, na.rm=T))
    ds$qc.upr <-  with(ds, qc.mean+qc.sd  )
    ds$qc.lwr <-  with(ds, qc.mean-qc.sd  )
    ds$qc.upr2 <- with(ds, qc.mean+2*qc.sd)
    ds$qc.lwr2 <- with(ds, qc.mean-2*qc.sd)
    dd <- merge(dat, ds)
    p  <- qcplot(dd, qcstatistic="TotalArea")
    p <- p + ggtitle(expression(atop(  "Total measured area for different runs together"
                                     , "with bands of one and two standard deviations")))
    plotlist <- dlply(dd, .(ProteinName), function(x) p %+% x)
    return(plotlist)
}

qcstat.BestRetentionTime <- function (dat) {
    ds <- ddply(dat, .(PeptideModifiedSequence, PrecursorMz), summarise,
                qc.mean = mean(BestRetentionTime, na.rm=T), qc.sd = sd(BestRetentionTime, na.rm=T))
    ds$qc.upr <-  with(ds, qc.mean+qc.sd  )
    ds$qc.lwr <-  with(ds, qc.mean-qc.sd  )
    ds$qc.upr2 <- with(ds, qc.mean+2*qc.sd)
    ds$qc.lwr2 <- with(ds, qc.mean-2*qc.sd)
    dd <- merge(dat, ds)
    p  <- qcplot(dd, qcstatistic="BestRetentionTime")
    p <- p + ggtitle(expression(atop(  "Best retention time for different runs together"
                                     , "with bands of one and two standard deviations")))
    plotlist <- dlply(dd, .(ProteinName), function(x) p %+% x)
    return(plotlist)
}

qcstat.MaxFwhm <- function (dat) {
    ds <- ddply(dat, .(PeptideModifiedSequence, PrecursorMz), summarise,
                qc.mean = mean(MaxFwhm, na.rm=T), qc.sd = sd(MaxFwhm, na.rm=T))
    ds$qc.upr <-  with(ds, qc.mean+qc.sd  )
    ds$qc.lwr <-  with(ds, qc.mean-qc.sd  )
    ds$qc.upr2 <- with(ds, qc.mean+2*qc.sd)
    ds$qc.lwr2 <- with(ds, qc.mean-2*qc.sd)
    dd <- merge(dat, ds)
    p  <- qcplot(dd, qcstatistic="MaxFwhm")
    p <- p + ggtitle(expression(atop(  "MaxFwhm for different runs together"
                                     , "with bands of one and two standard deviations")))
    plotlist <- dlply(dd, .(ProteinName), function(x) p %+% x)
    return(plotlist)
}

qcstat.MaxEndTime <- function (dat) {
    ds <- ddply(dat, .(PeptideModifiedSequence, PrecursorMz), summarise,
                qc.mean = mean(MaxEndTime, na.rm=T), qc.sd = sd(MaxEndTime, na.rm=T))
    ds$qc.upr <-  with(ds, qc.mean+qc.sd  )
    ds$qc.lwr <-  with(ds, qc.mean-qc.sd  )
    ds$qc.upr2 <- with(ds, qc.mean+2*qc.sd)
    ds$qc.lwr2 <- with(ds, qc.mean-2*qc.sd)
    dd <- merge(dat, ds)
    p  <- qcplot(dd, qcstatistic="MaxEndTime")
    p <- p + ggtitle(expression(atop(  "MaxEndTime for different runs together"
                                     , "with bands of one and two standard deviations")))
    plotlist <- dlply(dd, .(ProteinName), function(x) p %+% x)
    return(plotlist)
}

qcstat.AverageMassErrorPPM <- function (dat) {
    ds <- ddply(dat, .(PeptideModifiedSequence, PrecursorMz), summarise,
                qc.mean = mean(AverageMassErrorPPM, na.rm=T), qc.sd = sd(AverageMassErrorPPM, na.rm=T))
    ds$qc.upr <-  with(ds, qc.mean+qc.sd  )
    ds$qc.lwr <-  with(ds, qc.mean-qc.sd  )
    ds$qc.upr2 <- with(ds, qc.mean+2*qc.sd)
    ds$qc.lwr2 <- with(ds, qc.mean-2*qc.sd)
    dd <- merge(dat, ds)
    p  <- qcplot(dd, qcstatistic="AverageMassErrorPPM")
    p <- p + ggtitle(expression(atop(  "AverageMassErrorPPM for different runs together"
                                     , "with bands of one and two standard deviations")))
    plotlist <- dlply(dd, .(ProteinName), function(x) p %+% x)
    return(plotlist)
}

qcstat.IsotopeDotProduct <- function (dat) {
    ds <- ddply(dat, .(PeptideModifiedSequence, PrecursorMz), summarise,
                qc.mean = mean(IsotopeDotProduct, na.rm=T), qc.sd = sd(IsotopeDotProduct, na.rm=T))
    ds$qc.upr <-  with(ds, qc.mean+qc.sd  )
    ds$qc.lwr <-  with(ds, qc.mean-qc.sd  )
    ds$qc.upr2 <- with(ds, qc.mean+2*qc.sd)
    ds$qc.lwr2 <- with(ds, qc.mean-2*qc.sd)
    dd <- merge(dat, ds)
    p  <- qcplot(dd, qcstatistic="IsotopeDotProduct")
    p <- p + ggtitle(expression(atop(  "IsotopeDotProduct for different runs together"
                                     , "with bands of one and two standard deviations")))
    plotlist <- dlply(dd, .(ProteinName), function(x) p %+% x)
    return(plotlist)
}

ps1 <- try(qcstat.TotalArea(dat), silent=T)
ps2 <- try(qcstat.BestRetentionTime(dat), silent=T)
ps3 <- try(qcstat.MaxFwhm(dat), silent=T)
ps4 <- try(qcstat.MaxEndTime(dat), silent=T)
ps5 <- try(qcstat.AverageMassErrorPPM(dat), silent=T)

if (verbose) message("Opening PDF file for writing ", pdfFileName)

pdf(pdfFileName, height=11.6, width=8.2)
if(!inherits(ps1, "try-error")) print(ps1)
if(!inherits(ps2, "try-error")) print(ps2)
if(!inherits(ps3, "try-error")) print(ps3)
if(!inherits(ps4, "try-error")) print(ps4)
if(!inherits(ps5, "try-error")) print(ps5)
msg <- dev.off()

if (verbose) {
    message("Successfully closed PDF file ", pdfFileName)
}

if (verbose) {
    Sys.sleep(2)
    browseURL(pdfFileName)
}
