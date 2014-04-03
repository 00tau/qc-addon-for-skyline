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

####################################
## Reading command line arguments ##
####################################

arg <- commandArgs(TRUE)

fileCSV <- grep('.csv$', arg, value=TRUE)
verbose <- is.element("verbose", arg)
pdfFileName <- gsub('csv', 'pdf', fileCSV)

nump = 10
identifiers <- c("ProteinName", "PeptideModifiedSequence", "PrecursorMz", "AcquiredTime")

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

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(chron))

if (verbose) message("Successfully loaded packages: ggplot2, plyr")

##########################################
## Adding Optional Pseudo Protein Names ##
##########################################

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

##################
## Data parsing ##
##################

parseTime <- function (chr) {
    a <- as.character(chr)
    a <- matrix(unlist(strsplit(a, " ")), nrow=3)
    dts <- chron(a[1,], a[2,])
    dts <- dts + as.numeric(mapvalues(a[3,], c("AM", "PM"), c(0,1/2)))
    return(dts)
}

cleanup <- function (dat) {
    dat$ReplicateName <- NULL #factor(make.names(dat$ReplicateName))
    dat$PeptideModifiedSequence <- factor(gsub("\\]|\\[|\\+", "", dat$PeptideModifiedSequence))
    dat$AcquiredTime <- factor(parseTime(dat$AcquiredTime), ordered=TRUE)

    if (is.null(dat$ProteinName)) {
        dat$ProteinName <- createPseudoProteinNames(dat$PeptideModifiedSequence, nump)
    }
    return(dat)
}

dat <- cleanup(read.csv(fileCSV, na.strings="#N/A"))

if (verbose) {
    message("Successfully parsed report file having the following fields:")
    print(names(dat))
}

#######################
## Outlier detection ##
#######################

#' @export
grupp.test <- function(xvec, sgnf) {
    n <- length(na.omit(xvec))
    candidate <- abs(xvec - mean(xvec, na.rm=TRUE))
    gstat <- max(candidate, na.rm=TRUE) / sd(xvec, na.rm=TRUE)
    tstat <- qt(p=1-(sgnf / (2*n)), df=n-2)
    crval <- ((n-1)/sqrt(n)) * sqrt(tstat  / (n-2+tstat))
    if (gstat > crval) {
        return(which.max(candidate))
    } else {
        return(integer(0))
    }
}

#' @export
grupp.iterate <- function(xvec, sgnf) {
    pvec <- integer(0)
    while(sum(!is.na(xvec)) > 3 & !identical(sd(xvec, na.rm=T), 0)) {
        nvec <- grupp.test(xvec, sgnf)
        if (length(nvec) == 0) break()
        xvec[nvec] <- NA
        pvec <- c(pvec, nvec)
    }
    return(pvec)
}

#' @export
grupp <- function(xvec, sgnf) {
    outlier <- rep(FALSE, length(xvec))
    pvec <- grupp.iterate(xvec, sgnf)
    outlier[pvec] <- TRUE
    return(outlier)
}

###################################################
## Control Statistics and plotting functionality ##
###################################################

goodTheme <- theme(  axis.text.x = element_text(angle = 45, hjust = 1)
                   , legend.position="bottom"
                   , axis.text = element_text(colour = "black")
                   , panel.background = element_rect(fill='grey80', colour='white')
                   , panel.grid.major.x = element_line(size = 1.42)
                   , strip.text.y = element_text(size=6, angle=45)
                   #, plot.background = element_rect(fill='grey80', colour='black')
                   )

qcFlowChart <- function(dat, stat) {
    p <-(  ggplot(dat, aes(x=AcquiredTime, fill=quality))
         + geom_point(aes_string(y=stat), shape=21, size=3)
         + scale_shape_discrete(solid=F)
         + facet_grid(PeptideModifiedSequence+PrecursorMz~., scale="free")
         + geom_hline(aes(yintercept=mean), colour="black")
         + geom_hline(aes(yintercept=mean - 1*sd),  colour="green3")
         + geom_hline(aes(yintercept=mean + 1*sd),  colour="green3")
         + geom_hline(aes(yintercept=mean - 2*sd),  colour="yellow3")
         + geom_hline(aes(yintercept=mean + 2*sd),  colour="yellow3", show_guide=T)
         + scale_fill_manual(expression("Observation within")
         , values=c("one"="green", "two"="yellow", "more"="blue", "outlier"="red")
         , labels=c("one"="one std dev",
                    "two"="two std dev",
                    "more"="more that two std dev",
                    "outlier"="possible outlier")
         , breaks=c("one", "two", "more", "outlier"))
         )
    return(p)
}

qcstat <- function (dat, stat) {
    dat <- dat[c("ProteinName", "PeptideModifiedSequence", "PrecursorMz", "AcquiredTime", stat)]
    dat <- ddply(  dat
                  , .(PeptideModifiedSequence, PrecursorMz)
                  , function(d) {
                      tmp <- data.frame(  ProteinName = d$ProteinName
                                        , AcquiredTime = d$AcquiredTime)
                      tmp[[stat]] <- d[[stat]]
                      tmp$ugly <- grupp(d[[stat]], 0.0002)
                      return(tmp)
                  }
                  ) # do not particularly like this hack
    datm <- dat[dat$ugly == FALSE,]
    frm <- as.formula(paste0(stat, " ~ PeptideModifiedSequence + PrecursorMz"),
                             env=new.env())
    dat.mean <- aggregate(frm, data=datm, mean, na.rm=T)
    dat.sd <- aggregate(frm, data=datm, sd, na.rm=T)
    names(dat.mean)  [3] <- "mean"
    names(dat.sd)    [3] <- "sd"
    dat <- Reduce(merge, list(dat, dat.mean, dat.sd))
    good <- abs(dat[[stat]] - dat$mean) > dat$sd
    bad  <- abs(dat[[stat]] - dat$mean) > 2*dat$sd
    dat$quality <- factor(  pmax(good + bad, 3*dat$ugly)
                          , levels=0:3
                          , labels=c("one", "two", "more", "outlier"))
    return(dat)
}

customTitle <- function (stat) {
    return(ggtitle(paste(strwrap(paste(stat, "for different runs together with
                                       bands of one and two standard
                                       deviations"), width=72),
                                       collapse="\n")))
}

qcplot <- function (dat, stat) {
    dat <- qcstat(dat, stat)
    flowChart <- qcFlowChart(dat, stat) + customTitle(stat) + goodTheme
    ps <- dlply(dat, .(ProteinName), function(x) flowChart %+% x)
}

qcmake <- function (dat, ident) {
    plotMe <- function (l) {
        if (is.numeric(dat[[l]])) return(qcplot(dat, l))
    }
    lapply(setdiff(names(dat), ident), plotMe)
}

#########################################
## The actual plotting into a PDF file ##
#########################################

ps <- qcmake(dat, identifiers)

if (verbose) message("Opening PDF file for writing ", pdfFileName)

pdf(pdfFileName, height=11.6, width=8.2)
ps
dev.off()

if (verbose) {
    message("Successfully closed PDF file ", pdfFileName)
}

if (verbose) {
    Sys.sleep(2)
    browseURL(pdfFileName)
}
