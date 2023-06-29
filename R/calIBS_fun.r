#' To predict the heterosis of a given breed with other breeds
#'
#' @param x old score
#'
#' @return new score
#' @export
#'

#' @examples calIBS(snp_lst,breeds)

        
calIBS<-function(snp_lst,breeds){
    for (i in 1:length(snp_lst)){
     assign(breeds[i],snp_lst[[i]])
     }
    object_name<-breeds
    names<-NULL
    dup<-0
    for (i in 1:length(object_name)){
        name1<-substr(object_name[i],1,1)
        if (name1 %in% names){
            dup<-dup+1
            names<-c(names,paste0(name1,dup))
            }else{
            names<-c(names,name1)
            }
        }
    names(snp_lst)<-object_name
    len<-length(snp_lst)
    nmarker<-nrow(snp_lst[[1]])
    p<-lapply(snp_lst,function(x)freq(x))
    names(p)<-object_name
    com<-combn(1:len,2)
    com_name<-combn(names,2)
    for (i in 1:ncol(com)){
        p[[len+i]]<-freq2(p[[com[1,i]]],p[[com[2,i]]])
        names(p)[len+i]<-paste0(com_name[1,i],com_name[2,i])
        }
    
    temp<-data.frame()
    size<-length(p)
    for (i in 1:size){
        for (j in 1:size){
            if(i==j){
                temp[i,j]=fibs(p[[i]],p[[j]])
                }
            if(i<j){
                temp[i,j]=fibs(p[[i]],p[[j]])
                temp[j,i]=temp[i,j]
                }
            }
        }   
    row.names(temp)<-names(p)
    names(temp)<-names(p)
    return(temp)
    }

freq<-function(geno){  
        n0 <- apply(geno<=0.5,1,sum,na.rm=T)
        n1 <- apply(geno>0.5 & geno<=1.5,1,sum,na.rm=T)
        n2 <- apply(geno>1.5 & geno<=2,1,sum,na.rm=T)

        n <- n0 + n1 + n2

        ## calculate allele frequencies
        p <- ((2*n0)+n1)/(2*n)
        q <- 1 - p
        freq.mat <- data.frame(p, q)
        return(freq.mat)    
        }

freq2<-function(a,b){
        ## calculate the probability of three genotype for each locus
        temp<-cbind(a,b)
        freq.mat<-matrix(unlist(
                    apply(temp,1,
                    function(x){
                        d<-x[1]*x[3]
                        e<-x[2]*x[4]
                        f<-x[1]*x[4]+x[2]*x[3]
                        return(list((2 * d + f) / 2,
                                    (2 * e + f) / 2
                            ))
                        })),ncol=2,byrow=T)
        return(freq.mat)
        }

fibs<-function(a,b){
        ## calculate the IBS between all the probable population
        temp<-cbind(a,b)
        return(mean(apply(temp,1,function(x)1-(x[1]*x[3]+x[2]*x[4]))))
        }
