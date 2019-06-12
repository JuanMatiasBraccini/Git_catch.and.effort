VARIABLE.LIST=list(c("CPUE","FINYEAR"),c("CPUE","FINYEAR","new.level.vessel"),
          c("CPUE","FINYEAR","new.level.vessel","MONTH"),c("CPUE","FINYEAR","new.level.vessel","MONTH","Freo"),
          c("CPUE","FINYEAR","new.level.vessel","MONTH","Freo"),c("CPUE","FINYEAR","new.level.vessel","MONTH","Freo"),
          c("CPUE","FINYEAR","new.level.vessel","MONTH","Freo"),c("CPUE","FINYEAR","new.level.vessel","MONTH","Freo","zone"),
                   c("CPUE","FINYEAR","new.level.vessel","MONTH","Freo","zone"))




#K-fold test for model terms
fn.best.terms=function(DATA,VARIAB,FORMULA.pos,FORMULA.pos.log)
{
  # 1. Put data in proper shape
  DataFile=DATA[,match(VARIAB,names(DATA))]
  
  #drop levels not occurring in data
  for(f in 1:ncol(DataFile))
  {
    if(class(DataFile[,f])=='factor')DataFile[,f]=DataFile[,f][, drop=TRUE]
  }
  
  #select subsets
  id=round(nrow(DataFile)*training.data)
  
  list.ID=list()
  for(t in 1:K)list.ID[[t]]=sort(sample(1:nrow(DataFile),id))
  
  
  #2. loop for k-fold validation
  Correlation=vector('list',length=K)
    
    for(k in 1:K)
    {
      # 2.1. select data subsets
      All.data.train=DataFile[list.ID[[k]],]
      All.data.test=DataFile[-list.ID[[k]],]
      
      # create binary dataset
      Bi.train=All.data.train
      Bi.train[,1]=as.numeric(Bi.train[,1]>0)
      Bi.test=All.data.test
      Bi.test[,1]=as.numeric(Bi.test[,1]>0)
      
      # create positive dataset
      Pos.train=All.data.train[All.data.train[,1]>0,]
      Pos.test=All.data.test[All.data.test[,1]>0,]    
      
      
      
      #2.2. Fit different error structures
      
      #2.2.1   Binomial
      Binomial <- glm(FORMULA.pos, data=Bi.train, family="binomial", maxit=100)
      
      #2.2.2 Lognormal
      GLMlog <- glm(FORMULA.pos.log, data=Pos.train, family=gaussian, maxit=100)
      
      
      
      
      #2.3. Predict test.data
      
      IDD=which(All.data.test$CPUE>0)
      Pred.this=All.data.test[IDD,]
      
      #Binomial
      Bin.pred=predict(Binomial,newdata=Pred.this, type='response')
    
      
      #Lognormal
      biasCorr <- GLMlog$deviance/GLMlog$df.residual/2
      Log.pred=exp(predict(GLMlog,newdata=Pred.this,type='response')+biasCorr)
      Log.pred=Log.pred*Bin.pred
      
 
      
      
      #2.4. Measure fit
      #Correlation
      Correlation[[k]]=cor(Pred.this[,1],Log.pred,method='pearson')

    }
   
  return(list(Correlation=unlist(Correlation)))
}

Model.Term.Sel=vector('list',length=N.species)
for (j in 1:length(Species.list)))
  {
    Store.dummy=vector('list',length=length(form.bin))
    for(q in 1:length(form.bin))
    {
      Store.dummy[[q]]=fn.best.terms(DATA.list.LIVEWT.c[[j]],VARIABLE.LIST[[q]],form.bin[[q]],form.pos[[q]])
    }
    Model.Term.Sel[[j]]=Store.dummy
  }

a=matrix(unlist(Model.Term.Sel[[j]]),nrow=K)

for(q in 1:length(form.bin))
{
  mean(a[,q])
  sd(a[,q])
}