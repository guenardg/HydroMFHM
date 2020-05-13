##
### Auxiliary functions
##
### Connecting to the PostgreSQL database, opening a tcp tunnel if necessary.
connectToDB <-  function(cfg, remote=FALSE) {
  if(remote) {
    ## Connect on a remotely hosted database
    ## Opening a tcp tunnel for PostgreSQL
    "ssh -L %d:127.0.0.1:%d -N %s@%s &" %>%
      sprintf(
        cfg$psql$tcp_socket$client_port,
        cfg$psql$tcp_socket$host_port,
        cfg$ssh$user,
        cfg$ssh$host) %>%
      system
    Sys.sleep(1L)  ## Waiting for the tunnel to become available
    ##
    ## Opening a PostgreSQL connection through the tcp tunnel
    DBI::dbConnect(RPostgreSQL::PostgreSQL(),
                   host="localhost",
                   port=cfg$psql$tcp_socket$client_port,
                   user=cfg$psql$login$user,
                   password=cfg$psql$login$password,
                   dbname=cfg$psql$dbname$hydro)
    ##
  } else {
    ## Connect on a locally hosted database
    DBI::dbConnect(RPostgreSQL::PostgreSQL(),
                   host="localhost",
                   port=cfg$psql$tcp_socket$host_port,
                   user=cfg$psql$login$user,
                   password=cfg$psql$login$password,
                   dbname=cfg$psql$dbname$hydro)
  }
}
##
### Disconnecting from the PostgreSQL database, closing the tcp tunnel if
### necessary.
disconnectFromDB <-  function(psqlCon, remote=FALSE) {
  DBI::dbDisconnect(psqlCon)
  if(remote)
    "fuser -k %d/tcp" %>%
    sprintf(cfg$psql$tcp_socket$client_port) %>%
    system(intern = TRUE) %>%
    as.integer
}
##
gcd.hf2 <- function (from, to, radius = 6371) 
{
  if(is.null(ncol(from))||is.null(ncol(to)))
    stop("Both parameters 'from' and 'to' must have columns.")
  if((ncol(from)!=2L)||(ncol(to)!=2L))
    warning("Parameters 'from' and/or 'to' has more than 2 columns, only the first 2 are used.")
  rr <- matrix(NA, nrow(to), nrow(from), dimnames = list(rownames(to), rownames(from)))
  for (j in 1L:nrow(from)) {
    delta.lon <- (to[, 2L] - from[j, 2L]) * pi/180
    delta.lat <- (to[, 1L] - from[j, 1L]) * pi/180
    a <- sin(delta.lat/2)^2 + cos(from[j, 1L] * pi/180) * cos(to[, 1L] * pi/180) * sin(delta.lon/2)^2
    d <- numeric(length(a))
    for (i in 1L:length(a)) d[i] <- 2 * asin(min(1, sqrt(a[i]))) * radius
    rr[, j] <- d
  }
  return(rr)
}
##
PoissonPDF_log <- function(x,lambda) {
  ans <- rep(NA,length(x))
  whcalc <- (lambda!=0)&(x!=0)
  whzero <- (lambda==0)&(x==0)
  ans[lambda!=0] <- x[lambda!=0]*log(lambda[lambda!=0])-lambda[lambda!=0]-lgamma(x[lambda!=0]+1)
  ans[whzero] <- 0
  ans
}
##
realpha <- function(grp, cl, alpha, penalty.coef) {
  cat("Re-calculating glmnet with alpha:",alpha,"and par:",
      paste(penalty.coef,collapse=","))
  glm <-
    grp %>%
    parLapply(
    cl=cl,
    X=.,
    fun=function(X,alpha,penalty.coef) {
      penalty.factor <- c(0.5,(1+exp(-X$mm%*%penalty.coef))^-1)
      glmnet(
        y=X$y,
        x=X$descriptors[,-1L],
        alpha=alpha,
        family="poisson",
        penalty.factor=penalty.factor
      )
    },
    alpha=alpha,
    penalty.coef=penalty.coef
  )
  cat("\n")
  return(glm)
}
##
objf_lambda <- function(par,cl,grp,glm,obs,L_perfect,L_null_CV) {
  cat("Trying lambda =",par,"\n")
  prd <- c()
  for(i in 1L:length(grp)) {
    cat(i,"| ")
    prd %<>%
      c(exp(grp[[i]][["predictors"]]%*%as.matrix(coef(glm[[i]],s=exp(par[1L])))))
  }
  L_model_CV <- sum(PoissonPDF_log(x=obs,lambda=prd))
  Rsqr_pred_CV <- 1-(L_perfect-L_model_CV)/(L_perfect-L_null_CV)
  cat("Got:",Rsqr_pred_CV,"\n")
  return(-Rsqr_pred_CV)
}
##
objf_alpha <- function(par,cl,grp,obs,L_perfect,L_null_CV) {
  glm <- realpha(grp,cl,alpha=(1+exp(-par[1L]))^-1,penalty.coef=par[-1L])
  res <- optim(par=-5,f=objf_lambda,method="Brent",cl=cl,grp=grp,glm=glm,
               obs=obs,L_perfect=L_perfect,L_null_CV=L_null_CV,
               control=list(maxit=1000),lower=-10,upper=0)
  return(res$value)
}
##
Apply_nonzero_B <- function(x,Xterm,Zterm,fun)
  fun(x$reg$B[dat$X$term==Xterm,dat$Z$term==Zterm,drop=FALSE]!=0)
##
list_which_nonzero <- function(x, Xterm, Zterm, byrow=FALSE) {
  res <- data.frame(X=character(),Y=character(),Value=double())
  if(byrow) {
    rs <- rowSums(!!x$reg$B[dat$X$term==Xterm,dat$Z$term==Zterm,drop=FALSE])
    wnz <- names(which(rs!=0))
    for(i in wnz) {
      ## i=wnz[1L]
      sel <- x$reg$B[i,dat$Z$term==Zterm]
      wnz <- which(sel!=0)
      res %<>%
        rbind(
          data.frame(
            X=i,
            Z=names(wnz),
            Values=sel[wnz]
          )
        )
    }
  } else {
    cs <- colSums(!!x$reg$B[dat$X$term==Xterm,dat$Z$term==Zterm,drop=FALSE])
    wnz <- names(which(cs!=0))
    for(i in wnz) {
      ## i=wnz[1L]
      sel <- x$reg$B[dat$X$term==Xterm,i]
      wnz <- which(sel!=0)
      res %<>%
        rbind(
          data.frame(
            X=names(wnz),
            Z=i,
            Values=sel[wnz]
          )
        )
    }
  }
  rownames(res) <- NULL
  res
}
##
interaction_display <- function(dat, mod, X, Z, xlab, labels=X,
                                prob=0.95, plot=TRUE, rng_stretch=0.015,
                                yoffset=0.15, angle=30, length=0.1,
                                signif=rep(2L,length(X))) {
  qtn <-
    dat$X$data[,X] %>%
    apply(2L,quantile,prob=c(0.5*(1-prob),0.5,1-0.5*(1-prob)))
  tmp <- qtn %>% {matrix(NA,2*ncol(.),ncol(.))}
  for(i in 1L:ncol(qtn))
    for(j in 1L:ncol(qtn))
      tmp[2L*(i-1L)+(1L:2L),j] <- qtn[if(i==j) c(1L,3L) else rep(2L,2L),j]
  tmp <-
    tmp %*%
    dat[[mod]]$B[X,dat$Z$term!="Constant"&dat$Z$term!="Phylogeny"] +
    dat[[mod]]$B[dat$X$term=="Constant",Z]
  tmp %<>%
    matrix(
      length(.)/2L,2L,byrow=TRUE,
      dimnames=list(colnames(qtn),c("2.5%","97.5%")))
  ##
  if(plot) {
    plot(NA, xlim=range(tmp)+rng_stretch*c(-1,1), ylim=c(0.5,nrow(tmp)+0.5),
         yaxt="n", xlab=xlab, ylab="")
    axis(2L,at=1L:nrow(tmp),labels=labels,las=1L)
    for(i in 1L:nrow(tmp)) {
      arrows(x0=tmp[i,1L], x1=tmp[i,2L], y0=i, y1=i, angle=angle,
             length=length)
      points(x=tmp[i,1L],y=i,pch=21L,bg="black")
      text(x=tmp[i,1L],y=i+yoffset,
           labels=sprintf(sprintf("%%0.%df",signif[i]),qtn[1L,i]),
           adj=if(tmp[i,1L]<tmp[i,2L]) 0.75 else 0.25)
      text(x=tmp[i,2L],y=i+yoffset,
           labels=sprintf(sprintf("%%0.%df",signif[i]),qtn[2L,i]),
           adj=if(tmp[i,1L]<tmp[i,2L]) 0.25 else 0.75)
    }
    (qtn[2L,] %*%
        dat[[mod]]$B[X,dat$Z$term!="Constant"&dat$Z$term!="Phylogeny"] +
        dat[[mod]]$B[dat$X$term=="Constant",Z]) %>%
      abline(v=.,lty=3L)
  }
  return(invisible(list(qtn,tmp)))
}
##
Space_display <- function(dat, space_grid, col, mod, Z, digits, title,
                          xlab, ylab, xleg, yleg) {
  Bz <- dat[[mod]]$B[,Z]
  cte <-
    dat$X$data[,dat$X$term!="Constant"&dat$X$term!="Space"] %>%
    colMeans %*%
    Bz[dat$X$term!="Constant"&dat$X$term!="Space"] +
    dat[[mod]]$B[dat$X$term=="Constant",Z]
  dim(cte) <- NULL
  fit <-
    space_grid$emapscan %*%
    dat[[mod]]$B[dat$X$term=="Space",Z] + cte
  rng <-
    fit %>%
    {c(floor(10^digits*min(.)),ceiling(10^digits*max(.)))/10^digits}
  fit %>%
    matrix(length(space_grid$lat),length(space_grid$lon)) %>% t %>%
    image(
      z=.,x=space_grid$lon,y=space_grid$lat,asp=1,col=col,
      xlab=xlab,ylab=ylab,zlim=rng,axes=FALSE
    )
  axis(1L)
  axis(2L)
  Info$river %>% {points(x=.$lon, y=.$lat)}
  ##
  leg <- seq(rng[1L],rng[2L],10^-digits)
  legend(
    title=title, x=xleg,y=yleg, pch=22L, legend=sprintf("%.2f",leg),
    bg="white", ncol=2L, cex=1.25,
    pt.bg=col[1L+as.integer(seq(0,length(col)-1L,length.out=length(leg)))]
  )
}
##
phylo_display <- function(dat, phylo_grid, mod, X, labels, phylo.width=0.4,
                          cex.phylo=1, pch=21L, std=TRUE, ...) {
  fit <-
    phylo_grid$pem %*%
    (dat[[mod]]$B[X,dat$Z$term=="Phylogeny",drop=FALSE] %>% t) +
    rep(dat[[mod]]$B[X,dat$Z$term=="Constant"],each=nrow(phylo_grid$pem))
  stdev <- dat$X$data %>% apply(2L,function(x) sqrt(var(x)))
  bnd <- seq(phylo.width,1,length.out=length(X)+1L)
  ##
  par_safe <- par(no.readonly = TRUE)
  mar_tmp <- par_safe$mar
  mar_tmp[4L] <- 0
  par(mar=mar_tmp,fig=c(0,phylo.width,0,1))
  phylo_grid$tree %>% plot(cex=cex.phylo)
  mar_tmp <- par_safe$mar
  mar_tmp[2L] <- 0
  for(i in 1L:length(X)) {
    par(mar=mar_tmp,fig=c(bnd[i],bnd[i+1L],0,1),new=TRUE)
    if(std) {
      plot(x=fit[,X[i]]*stdev[X[i]],y=1L:nrow(fit),axes=FALSE,ylab="",xlab="",
           main=if(missing(labels)) X[i] else labels[i],xpd=TRUE,pch=pch,...)
      axis(1L,las=2L)
      abline(v=dat[[mod]]$B[X[i],dat$Z$term=="Constant"]*stdev[X[i]],lty=3L)
    } else {
      plot(x=fit[,X[i]],y=1L:nrow(fit),axes=FALSE,ylab="",xlab="",
           main=if(missing(labels)) X[i] else labels[i],xpd=TRUE,pch=pch,...)
      axis(1L,las=2L)
      abline(v=dat[[mod]]$B[X[i],dat$Z$term=="Constant"],lty=3L)
    }
  }
  par(par_safe)
  return(invisible(NULL))
}
##
