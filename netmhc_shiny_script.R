library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(reshape2)
palette(brewer.pal(9,'Spectral'))



#Get the affinities for a 15mer for the alleles in a given file of allele percentages
#This will store the affinities in a mysql database if not already there and return
# a data frame of alleles, 15mers and percent rank affinities

get_affinities<-function(sequence15,allele_file,pool){
  
  conn<-poolCheckout(pool)
  dbBegin(conn)
  allele_list=read.csv(allele_file,stringsAsFactors = FALSE)$Allele
    
  query<-paste0("select * from affinity_storage where Peptide like '",sequence15,"' and Allele in ('",paste(allele_list,collapse = "','"),"');")
  dbres<-dbGetQuery(conn,query)
  
  if (dim(dbres)[1]!=length(allele_list)){
    fileConn<-file("temp.fasta")
    writeLines(c(">temp",sequence15), fileConn)
    close(fileConn)
  
   
    allele_list=read.csv(allele_file,stringsAsFactors = FALSE)$Allele
    #print(allele_list)
    bash_command = paste0('for allele in ',paste(allele_list,collapse = ' '),'; do /home/ubuntu/netMHCIIpan-3.2/netMHCIIpan -f temp.fasta -a $allele ; done')

    netmhc_ii_results<-system(bash_command,intern=TRUE)
    netmhc_ii_results2<-strsplit(netmhc_ii_results,'\\s+')
    print(netmhc_ii_results2)
    filtered_results<-data.frame('Allele'=NA,'Peptide'=NA,'percent_rank'=NA)
    i=1
    for (i in 1:length(netmhc_ii_results2)){
  
      try(   if(netmhc_ii_results2[[i]][3]!='Allele' & !(is.na(netmhc_ii_results2[[i]][3])) & length(netmhc_ii_results2[[i]])%in%12:13){
        filtered_results[i,]<-netmhc_ii_results2[[i]][c(3,4,11)]
        })
      i<-i+1
      }
    filtered_results<-filtered_results[complete.cases(filtered_results),]
    dbWriteTable(conn,value=filtered_results,name='affinity_storage',append=TRUE,row.names=FALSE)
    } 
 
  result<-dbGetQuery(conn,query)
  dbCommit(conn)
  poolReturn(conn)
  return(result)
}

# Take a sequence and create the heatmaps and promiscuity plots
# for each one

make_matrix<-function(full_seq,allele_file,pool){
  fifteen_mers<-character()
  allele_list=read.csv(allele_file,stringsAsFactors = FALSE)$Allele
  
  for (i in 1:(nchar(full_seq)-14)){
    fifteen_mers<-c(fifteen_mers,substr(full_seq,i,i+14))
    }

  affinity_matrix<-matrix(NA,nrow=length(allele_list),ncol=length(fifteen_mers),dimnames = list('Allele'=allele_list,"Sequence"=fifteen_mers))
    
  for (i in fifteen_mers){
    temp_affinity_scores<-get_affinities(i,allele_file,pool)
    print(temp_affinity_scores[1,])
    for (i in 1:dim(temp_affinity_scores)[1]){
      affinity_matrix[temp_affinity_scores[i,1],temp_affinity_scores[i,2]]<-as.numeric(temp_affinity_scores[i,3])
      }
    }
  return(t(affinity_matrix))
  }

plot_matrix_threshold<-function(affinity_matrix,threshold,legend_position){
  affinity_matrix[affinity_matrix>threshold]<-NA
  aff_mat_thr<-melt(affinity_matrix)
  aff_mat_thr$Allele<-factor(aff_mat_thr$Allele,levels=rev(levels(aff_mat_thr$Allele)))
  af_mat<-ggplot(aes(x=Sequence,y=Allele,fill=value),data=data.frame(aff_mat_thr))+geom_tile(colour=rgb(0,0,0,.05),size=.2)+scale_x_discrete(labels=unlist(substring(rownames(affinity_matrix),8,8)))
  af_mat<-af_mat+scale_fill_gradient(low='red',high='lightgoldenrod1',na.value='white')
  af_mat<-af_mat+theme(panel.background = element_blank(),legend.text = element_text(size=16),axis.text = element_text(size=16),axis.title = element_text(size=20),
                       legend.title = element_text(size=20),legend.position = legend_position)
  if (legend_position=='bottom'){
  af_mat<-af_mat+guides(fill=guide_colorbar(title='Percent Rank',barwidth = unit(20,'lines'),title.position='bottom',title.hjust = .5))
  af_mat<-af_mat+theme(axis.text=element_text(size=11))
  }
  if (legend_position=='right'){
    af_mat<-af_mat+guides(fill=guide_colorbar(title='Percent Rank',barheight = unit(20,'lines'),title.position='top',title.hjust = .5))
    af_mat<-af_mat+theme(axis.text=element_text(size=13))
    }
  af_mat<-af_mat+labs(x='')
  
  return(af_mat)
  }

plot_single_promiscuity<-function(prom_matrix,legend_logical,axis_sequence){
  dat<-melt(prom_matrix)
  dat$var2<-factor(dat$Var2)
  prom_plot<-ggplot(data=dat,aes(x=Var1,y=value,group=Var2,colour=Var2))+geom_line(show.legend = legend_logical,size=2,alpha=.4,position = position_dodge(width=.4))
  prom_plot<-prom_plot+theme_classic()
  prom_plot<-prom_plot+scale_y_continuous(breaks=c(0,.2,.4,.6,.8,1),labels=c('0','20','40','60','80','100'),limits = c(0,1))
  prom_plot<-prom_plot+scale_x_continuous(breaks=1:dim(prom_matrix)[1],labels=axis_sequence)
  prom_plot<-prom_plot+theme(legend.position = 'bottom',legend.text=element_text(size=16),axis.text=element_text(size=16),
                               axis.title=element_text(size=16))
  prom_plot<-prom_plot+labs(x='',y='Percent Allele Promiscuity')
  return(prom_plot)
  }

plot_multi_promiscuity<-function(multi_prom_melted){
  mpplot<-ggplot(data=multi_prom_melted,aes(x=Position,y=Promiscuity,group=Peptide,colour=Peptide))+geom_line(size=2,alpha=.4,position = position_dodge(width=.5),na.rm = TRUE)
  mpplot<-mpplot+theme_classic()+scale_x_continuous(breaks=1:nchar(m_align[1]),labels=unlist(strsplit(consensusString(m_align),'')))
  mpplot<-mpplot+theme(legend.position = 'bottom',legend.text=element_text(size=16),axis.text=element_text(size=16),
                       axis.title=element_text(size=16))
  mpplot<-mpplot+guides(colour=guide_legend(title='',nrow = as.integer(length(unique(multi_prom_melted$Peptide))/1))) 
  mpplot<-mpplot+labs(x='',y='Percent Allele Promiscuity')+scale_y_continuous(breaks=c(0,.2,.4,.6,.8,1),labels=c(0,20,40,60,80,100),limits=c(0,1))
  return(mpplot)
  }

