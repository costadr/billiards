orbita <- function(r,delta,teta0,alfa0,max_hits_circle,max_hits_square){
  #here we find the intercept of the particle with one of the line segments
  find_intercept<-function(){
    #intercept with line L1
    xe=-1
    dte=(xe-x0)/cos(mu0)
    if(dte>0){
      ye=y0+sin(mu0)*dte
      if(ye>-1 && ye<1){
        #reached line segment L1
        lines(c(x0,xe),c(y0,ye))
        x0 <<- xe
        y0 <<- ye
        mu0 <<- (-mu0+2*(3*pi/2)) %% (2*pi)
        return()
      }
    }
    #intercept with line L2
    ye=-1
    dte=(ye-y0)/sin(mu0)
    if(dte>0){
      xe=x0+cos(mu0)*dte
      if(xe>-1 && xe<1){
        #reached line segment L2
        lines(c(x0,xe),c(y0,ye))
        x0 <<- xe
        y0 <<- ye
        mu0 <<- (-mu0+2*(0)) %% (2*pi)
        return()
      }
    }
    #intercept with line L3
    xe=+1
    dte=(xe-x0)/cos(mu0)
    if(dte>0){
      ye=y0+sin(mu0)*dte
      if(ye>-1 && ye<1){
        #reached line segment L3
        lines(c(x0,xe),c(y0,ye))
        x0 <<- xe
        y0 <<- ye
        mu0 <<- (-mu0+2*(pi/2)) %% (2*pi)
        return()
      }
    }
    #intercept with line L4
    ye=+1
    dte=(ye-y0)/sin(mu0)
    if(dte>0){
      xe=x0+cos(mu0)*dte
      if(xe>-1 && xe<1){
        #reached line segment L4
        lines(c(x0,xe),c(y0,ye))
        x0 <<- xe
        y0 <<- ye
        mu0 <<- (-mu0+2*(pi)) %% (2*pi)
        return()
      }
    }
    print("Some error has occurred! It is necessary to find an intercept")
  }
  
  #drawing billiard boundary
  teta<-seq(0,2*pi,length.out=1000)
  x=delta+r*cos(teta)
  y=r*sin(teta)
  par(mar=c(2,2,0.1,0.1))
  plot(x,y,xlim=c(-1,1),ylim=c(-1,1),type='l',col='red')
  lines(c(-1,-1,+1,+1,-1),c(-1,+1,+1,-1,-1),col='red')
  
  #initial considerations
  x0=delta+r*cos(teta0)
  y0=r*sin(teta0)
  xl=-r*sin(teta0)
  yl=+r*cos(teta0)
  fi0=(pi+atan2(yl,xl)) %% (2*pi)
  mu0=(alfa0+fi0) %% (2*pi) #angle that gives the particle's direction
  
  nhits_circle=0
  for(ntimes in seq(1,max_hits_square)){
    #if a particle starts in the circle, it cannot hit the circle again
    #i.e., the particle needs to reach one of the line segments (L1, L2, L3 or L4)
    find_intercept()
    
    #now it is necessary to verify if the particle touches the circle
    b=2*((x0-delta)*cos(mu0)+y0*sin(mu0))
    c=(x0-delta)**2+y0**2-r**2
    del=b**2-4*c
    if(del>0){
      #particle touches the circle
      dte=(-b-sqrt(del))/2
      x1=x0+cos(mu0)*dte
      y1=y0+sin(mu0)*dte
      teta1=atan2(y1,x1-delta)
      xl=-r*sin(teta1)
      yl=+r*cos(teta1)
      fi1=(pi+atan2(yl,xl)) %% (2*pi)
      alfa1=(fi1-mu0) %% (pi)
      
      nhits_circle=nhits_circle+1
      if(nhits_circle==max_hits_circle){
        return()
      }
      
      lines(c(x0,x1),c(y0,y1))
      
      x0=x1
      y0=y1
      mu0=(alfa1+fi1) %% (2*pi)
    }
  }
}

orbita(r=0.6,delta=0.2,teta0=0.1,alfa0=0.1,max_hits_circle=2,max_hits_square=2)
