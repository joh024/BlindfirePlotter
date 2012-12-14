BFPenv = new.env()

BlindfirePlot =
  ## -respform-
  ## response variable (y), formula
  ## can also be specified as a character vector (of length > 1),
  ##   which will be collapsed into a formula
  ## e.g. c("A", "B") becomes A ~ B
  ##
  ## -orfac-
  ## ordered explanatory variables
  ## factors, but has ordered meaning
  ## maximum 1 for any given analysis
  ##  (any more can be considered a unorfac)
  ##
  ## -unorfac-
  ## unordered explanatory variables
  ## factors, with no ordered meaning
  ## maximum 2 for any given analysis
  ##  (including any flow-over from
  ##   multiple respnum, orfac)
  ##
  ## -scaleresp-
  ## scaling variable for each resp element, formula
  ## specified as ~ expr
  ##
  ## -datdf-
  ##
  ##
  ## -returnind & returnfac-
  ##
  ##
  ## -pdfname-
  ## The name of the pdf file to generate.
  ## Special cases:
  ##   NULL (default): Auto-generate a pdfname
  ##   NA: Do not generate a pdf file
  function(respform, orfac = NULL, unorfac = NULL, datdf,
           scaleresp = NULL, returnind = NULL, returnfac = NULL,
           pdfname = NULL, width = 16 * 0.9, height = 9 * 0.9,
           runtime = TRUE, plotfunc = NULL){
    if(runtime == TRUE) starttime = proc.time()
        
    ## if respform is a character vector, convert to a formula
    if(is.character(respform) && length(respform) > 1)
      respform = parse(text =
        paste(respform, collapse = "~"))[[1]]
    
    ## Make a nice character vector to represent respform
    ##   for use in graph titles
    if(is.null(scaleresp))
      respchar = paste(deparse(substitute(respform)), "\n")
    else
      respchar = paste(deparse(substitute(respform)),
        "scaled by", deparse(substitute(scaleresp)), "\n")
    ## Parse the formula
    respnum = BFPparseFormula(respform, scaleresp)

    ## Compute the number of factors
    numfacs = (if(length(respnum) > 1) 1 else 0) +
      length(orfac) + length(unorfac)

    ## Construct code fragments
    codefrags = list()
    for(i in 1:length(respnum))
      codefrags[[i]] = substitute(eval(respnum[[i]]$expr,
                 lapply(respnum[[i]]$els, function(x)
                        datagrab(datdf, eval(x, datdf), facs))
                 ), list(i = i))


    ## Open pdf
    if(is.null(pdfname))
      pdfname = paste("BlindfirePlotResult-",
        format(Sys.time(), "%H.%M.%S-%d-%b-%Y"), ".pdf", sep = "")
    if(!is.na(pdfname))
      pdf(pdfname, width = width, height = height)

    ## Clear BFPenv
    rm(list = ls(envir = BFPenv), envir = BFPenv)
    ## Keep track of curind in BFPenv
    assign("curind", 1, envir = BFPenv)
    assign("usedind", 0, envir = BFPenv)
    
    ## Single variable plots
    if(length(respnum) == 1 && numfacs >= 1){
      BFPuI(respnum, respchar, orfac, unorfac, datdf,
        codefrags, returnind, returnfac, plotfunc[["uI"]])
    }

    ## Two variable plots
    if(!is.null(orfac) && numfacs >= 2){
      BFPoIuI(respnum, respchar, orfac, unorfac, datdf,
        codefrags, returnind, returnfac, plotfunc[["oIuI"]])
    }
    # if(numfacs >= 2){
      # BFPuII(respnum, respchar, orfac, unorfac, datdf,
        # codefrags, returnind, returnfac, plotfunc[["uII"]])
    # }
    
    ## Three variable plots
    if(!is.null(orfac) && numfacs >= 3){
      BFPoIuII(respnum, respchar, orfac, unorfac, datdf,
        codefrags, returnind, returnfac, plotfunc[["oIuII"]])
    }
    
    ## Close pdf
    if(!is.na(pdfname))
      dev.off()

    ## Print runtime
    if(runtime == TRUE)
      print(structure(proc.time() - starttime, class = "proc_time"))
    
    ## Return datasets
    rm("curind", "usedind", envir = BFPenv)
    invisible(as.list(BFPenv))
  }

BFPparseFormula =
  ## Blindfire Plot supporting function for
  ##   parsing of respform formula
  local({
    ##########################
    ## Supporting functions ##
    SplitByTilde =
      ## Split given expr by `~`
      ##  running recursively as necessary
      ## Returns list, each element containing
      ##  the parts of expr as split by `~`
      function(expr){
        if(is.call(expr)){
          if(expr[[1]] == as.name("~")){
            if(length(expr) == 3)
              ## of form a ~ b
              c(SplitByTilde(expr[[2]]), list(expr[[3]]))
            else
              ## of form ~ b
              list(expr[[2]])
          } else list(expr)
        } else list(expr)
      }

    GrabElements =
      ## Grab the individual elements
      ##  of the given expr
      function(expr){
        if(is.call(expr)){
          if(length(expr) == 3)
            c(GrabElements(expr[[2]]), GrabElements(expr[[3]]))
          else
            GrabElements(expr[[2]])
        } else if(is.name(expr))
          list(expr)
        else NULL
      }

    ###############################
    ## Main function starts here ##
    function(respform, scaleresp){
      ## Split by ~
      exprlist = SplitByTilde(respform)

      ## Add scaling
      if(!is.null(scaleresp)){
        if(is.call(scaleresp))
          if(scaleresp[[1]] == as.name("~"))
            scaleresp = scaleresp[[2]]
        exprlist = lapply(exprlist, function(x)
          call("/", x, scaleresp))
      }
      names(exprlist) = as.character(exprlist)

      ## Add els list
      comblist = lapply(exprlist, function(x){
        els = GrabElements(x)
        els = els[!duplicated(as.character(els))]
        names(els) = as.character(els)
        list(expr = x, els = els)
      })

      ## Return combined list
      comblist
    }})

BFPinnerfunc =
  ## The `inner' function for the BFP plotting
  ##   sub-functions, handling routine procedures
  ##   common across all sub-functions
  function(respnum, respchar, facs, datdf,
           codefrags, returnind, returnfac, plotfunc){
    ## Grab current curind from BFPenv and increment it by 1
    curind = get("curind", envir = BFPenv)
    assign("curind", curind + 1, envir = BFPenv)

    ## If check passes, eval and draw
    if(BFPretcheck(returnind, returnfac, curind, facs)){
      if(length(codefrags) == 1){
        curdat = eval(codefrags[[1]], list(facs = facs))
        if(length(facs) == 3){
          rawdat = curdat
          curdat = list()
          for(cursec in dimnames(rawdat)[[3]])
            curdat = c(curdat, list(rawdat[,,cursec]))
          names(curdat) = dimnames(rawdat)[[3]]
        }
      } else{
        datraw = lapply(codefrags, function(x)
          eval(x, list(facs = facs)))
        names(datraw) = names(respnum)
        if(length(facs) == 1)
          curdat = do.call("cbind", datraw)
        else if(length(facs) == 2)
          curdat = lapply(datraw, function(x) x[,])
      }
      plotfunc(curdat)
      title(main = paste(respchar, "by",
              paste(facs, collapse = " and ")))
      title(main = paste("index", curind),
            adj = 0, cex.main = 0.8, font.main = 1)

      ## Save restructured data to BFPenv
      assign(paste("index", curind, sep = "."), curdat, envir = BFPenv)
    }
  }

BFPuI =
  ## Blindfire Plot supporting function for
  ##   one variable plots, with 1 unorfac
  ## Canonical form:
  ## structure(rbinom(5, 10, 0.5), .Names = sample(LETTERS, 5))
  ## only works for length(respnum) == 1
  function(respnum, respchar, orfac, unorfac, datdf,
           codefrags, returnind, returnfac, plotfunc){
    if(is.null(plotfunc))
      plotfunc = barplot
    
    ## loop first by orfac
    if(!is.null(orfac))
      for(curorf in 1:length(orfac)){
        facs = orfac[curorf]
        BFPinnerfunc(respnum, respchar, facs, datdf,
                     codefrags, returnind, returnfac, plotfunc)
      }
    ## then, unorfac
    if(!is.null(unorfac))
      for(curuno in 1:length(unorfac)){
        facs = unorfac[curuno]
        BFPinnerfunc(respnum, respchar, facs, datdf,
                     codefrags, returnind, returnfac, plotfunc)
      }

    ## Save indices to BFPenv
    BFPvecstore("uI")
  }

BFPoIuI =
  ## Blindfire Plot supporting function for
  ##   two variable plots, with 1 orfac and 1 unorfac
  ## Canonical form:
  ## matrix(rbinom(9, 10, 0.5), nrow = 3,
  ##        dimnames = list(orfac = LETTERS[1:3],
  ##        unorfac = c("Baseball", "Football", "Golf")))
  ## where length(respnum) > 1, unorfac will always be
  ##   names(respnum)
  function(respnum, respchar, orfac, unorfac, datdf,
           codefrags, returnind, returnfac, plotfunc){
    if(is.null(plotfunc))
      plotfunc = lineplot
    
    for(curorf in 1:length(orfac)){
      if(length(respnum) > 1){
        ## our only unorfac will be names(respnum)
        facs = orfac[curorf]
        BFPinnerfunc(respnum, respchar, facs, datdf,
                     codefrags, returnind, returnfac, plotfunc)
      } else{
        ## loop through other factors
        ## first, remaining orfacs
        ##   orfac[-curorf]
        combfacs = NULL
        if(length(orfac) > 1)
          combfacs = c(combfacs, orfac[-curorf])
        if(!is.null(unorfac))
          combfacs = c(combfacs, unorfac)
        for(posscomb in combfacs){
          facs = c(orfac[curorf], posscomb)
          BFPinnerfunc(respnum, respchar, facs, datdf,
                       codefrags, returnind, returnfac, plotfunc)
        }
      }
    }
    
    ## Save indices to BFPenv
    BFPvecstore("oIuI")
  }

BFPuII =
  ## Blindfire Plot supporting function for
  ##   two variable plots, with 2 unorfac
  ## Canonical form:
  ## matrix(rbinom(9, 10, 0.5), nrow = 3,
  ##        dimnames = list(unorfacA = c("NZ", "Aus", "USA"),
  ##        unorfacB = c("Baseball", "Football", "Golf")))
  ## where length(respnum) > 1, 1 of the 2 unorfac
  ##   will always be names(respnum)
  function(respnum, respchar, orfac, unorfac, datdf,
           codefrags, returnind, returnfac, plotfunc){
    if(is.null(plotfunc))
      plotfunc = testplot
    
    if(length(respnum) > 1){
      ## SKIPPED FOR NOW
      print("SKIPPED, BFPuII with respnum > 1")
    } else{
      combfacs = c(orfac, unorfac)
      possfacs = outer(combfacs, combfacs, function(x, y)
        paste("c('", x, "', '", y, "')",
              sep = ""))[!diag(length(combfacs))]
      for(posscomb in possfacs){
        facs = eval(parse(text = posscomb))
        BFPinnerfunc(respnum, respchar, facs, datdf,
                     codefrags, returnind, returnfac, plotfunc)
      }
    }

    ## Save indices to BFPenv
    BFPvecstore("uII")
  }

BFPoIuII =
  ## Blindfire Plot supporting function for
  ##   three variable plots, with 1 orfac and 2 unorfac
  ## Canonical form:
  ##  list(Australia =
  ##    matrix(rbinom(9, 10, 0.5), nrow = 3,
  ##           dimnames = list(orfac = LETTERS[1:3],
  ##           unorfac = c("Baseball", "Football", "Golf"))),
  ##       New.Zealand =
  ##    matrix(rbinom(9, 10, 0.5), nrow = 3,
  ##           dimnames = list(orfac = LETTERS[1:3],
  ##           unorfac = c("Baseball", "Football", "Golf"))))
  function(respnum, respchar, orfac, unorfac, datdf,
           codefrags, returnind, returnfac, plotfunc){
    if(is.null(plotfunc))
      plotfunc = function(datlist) lineplot(tlist(datlist))
    
    for(curorf in 1:length(orfac)){
      if(length(respnum) > 1){
        ## one of our unorfac is names(respnum)
        ## loop through other factors
        ## first, remaining orfacs
        ##   orfac[-curorf]
        combfacs = NULL
        if(length(orfac) > 1)
          combfacs = c(combfacs, orfac[-curorf])
        if(!is.null(unorfac))
          combfacs = c(combfacs, unorfac)
        for(posscomb in combfacs){
          facs = c(orfac[curorf], posscomb)
          BFPinnerfunc(respnum, respchar, facs, datdf,
                       codefrags, returnind, returnfac, plotfunc)
        }
      } else{
        combfacs = c(orfac[-curorf], unorfac)
        possfacs = outer(combfacs, combfacs, function(x, y)
          paste("c('", orfac[curorf], "', '", x, "', '", y, "')",
                sep = ""))[!diag(length(combfacs))]
        for(posscomb in possfacs){
          facs = eval(parse(text = posscomb))
          BFPinnerfunc(respnum, respchar, facs, datdf,
                       codefrags, returnind, returnfac, plotfunc)
        }
      }
    }

    ## Save indices to BFPenv
    BFPvecstore("oIuII")
  }

BFPretcheck =
  ## Logical check on return requirements being met
  ## Return TRUE IF:
  ##   No return requirements (both NULL)
  ## OR
  ##   Either curind can be found in returnind
  ##     or any facs can be found in returnfac
  ##   (hence the length of the intersections > 0)
  ## Otherwise, return FALSE.
  function(returnind, returnfac, curind, facs){
    if(is.null(returnind) && is.null(returnfac))
      TRUE
    else
      length(c(intersect(curind, returnind),
               intersect(facs, returnfac))) > 0
}

BFPvecstore =
  ## Checks the last used index (usedind)
  ##   against current index (curind)
  ## Saves sequence vector of used indices on most recent BFP
  ##   plot run to variable of name "storename" in BFPenv
  function(storename){
    curind = get("curind", envir = BFPenv)
    usedind = get("usedind", envir = BFPenv)
    if(curind - 1 > usedind){
      assign(storename, (usedind + 1):(curind - 1), envir = BFPenv)
      assign("usedind", curind - 1, envir = BFPenv)
    }
  }

datagrab =
  ## A wrapper function of sorts for `by`
  ## Biggest thing it does is to automatically convert the resulting
  ##   array into a vector, matrix or list of matrices, based on
  ##   the length of x (and hence the dimensions of the resulting array)
  ## y is actual vector of data
  ## x is vector of column names
  function(data, y, x, FUN = sum){
    res = by(y, data[,x], FUN)
    if(length(x) == 1)
      res = structure(as.numeric(res), .Names = dimnames(res)[[1]])
    else if(length(x) == 2)
      res = res[,]
    res
  }

tlist = function(datlist){
    datlist = lapply(datlist, as.matrix)
    newlist = list()
    for(i in 1:ncol(datlist[[1]]))
      newlist[[i]] = sapply(datlist, function(x) x[,i])
    names(newlist) = dimnames(datlist[[1]])[[2]]
    newlist
  }

lineplot = 
  ## A barplot with beside = TRUE
  ## but as a line plot
  function(heightlist, legendpos = "topright"){
    if(!is.list(heightlist)) heightlist = list(heightlist)

    ## orfac legend width calc
    orfaclegs = rownames(heightlist[[1]])
    widthaslines = max(strwidth(c("Orfac", orfaclegs),
      units = "inches")/par("csi") + 1)
    newmar = par("mar")
    newmar[4] = widthaslines
    opar = par(mar = newmar)
    
    ## Treat Inf as NA
    heightlist = lapply(heightlist, function(x) {x[x == Inf] = NA; x})
    
    maxval = sapply(heightlist, function(x) max(x, na.rm = TRUE))
    maxind = which(maxval == max(maxval, na.rm = TRUE))[1]
    if(length(maxind) == 0) maxind = 1
    barplot(heightlist[[maxind]] * 1.03,
            beside = TRUE, col = NA, border = NA)

    npts = length(heightlist)
    if(npts > 2){
      lty = rep(1, length = npts)
      bg = hsv(1:npts/npts, s = 0.5, v = 0.8)
      col = hsv(1:npts/npts, s = 0.5, v = 0.8)
    } else{
      lty = 1:2
      bg = c("black", "white")
      col = c("black", "black")
    }
    
    for(i in 1:length(heightlist)){
      height = heightlist[[i]]
      nbars = nrow(height)
      abline(v = 0.5, lty = 3)
      for(j in 1:ncol(height)){
        abline(v = (nbars + 1) * j + 0.5, lty = 3)
        lines(1:nbars + 0.5 + (nbars + 1) * (j - 1), height[,j],
              lty = lty[i], col = col[i])
        points(1:nbars + 0.5 + (nbars + 1) * (j - 1), height[,j],
               pch = 21, bg = bg[i])
      }
    }

    ## orfac legend plot
    mtext(paste("Orfac\n\n",
                paste(orfaclegs, collapse = "\n"), sep = ""),
          side = 4, las = 1, line = 0.5)

    ## superpos legend
    if(length(heightlist) > 1)
      legend(legendpos, legend = names(heightlist),
             lty = lty, pt.bg = bg, pch = 21, bg = "white")

    ## return par to normal
    par(opar)
  }

testplot =
  function(datamat){
    xvec = seq(from = 1, to = ncol(datamat))
    yvec = seq(from = 1, to = nrow(datamat))

    plot.new()
    plot.window(xlim = c(0, ncol(datamat)), xaxs = "i",
                ylim = c(0, nrow(datamat)), yaxs = "i")
    box()

    sum.narm = function(...) sum(..., na.rm = TRUE)
    ## Sort by order of marginals
    datamat = datamat[order(apply(datamat, 1, sum.narm)),
      order(apply(datamat, 2, sum.narm))]
    
    axis(1, at = xvec - 0.5, labels = colnames(datamat))
    axis(2, at = yvec - 0.5, labels = rownames(datamat))

    minval = min(datamat, na.rm = TRUE)
    maxval = max(datamat, na.rm = TRUE) - minval

    gray2 = function(level)
      if(is.na(level)) "#CC8F8F" else gray(level)
    
    for(j in xvec)
      for(i in yvec)
        rect(j - 1, i - 1, j, i,
             col = gray2(1 - (datamat[i, j] - minval)/maxval))
  }

testplot2 =
  function(datamat, ...){
    sum.narm = function(...) sum(..., na.rm = TRUE)
    ## Sort by order of marginals
    datamat = datamat[order(apply(datamat, 1, sum.narm)),
      order(apply(datamat, 2, sum.narm))]
    persp(z = datamat, ...)
  }
