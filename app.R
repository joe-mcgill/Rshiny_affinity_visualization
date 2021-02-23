library(RMariaDB)
library(pool)
library(shinythemes)
#lapply( dbListConnections( dbDriver( drv = "MariaDB")), dbDisconnect)
pool <- dbPool(
  drv = RMariaDB::MariaDB(),
  dbname = "shiny_affinities",
  host = "shiny-aff-1.c38jawccscp1.us-east-1.rds.amazonaws.com",
  username = "admin_shiny_aff",
  password = "5Yn2*$pxqRnK#p%8"
)



library(shiny)
library(Biostrings)
library(muscle)
source('netmhc_shiny_script.R')
library(xtable)
#library(shinythemes)
setClass('peptide',slots=c(sequence='character',affinities_matrix='matrix',blanks='numeric'))



ui <- fluidPage(
  theme='shiny_aff.css',
  navbarPage(

  title="Affinity Visualization Comparison",
  tabPanel('Home',
	fluidPage(
		includeHTML('www/home.html')
	)
  ),
  
  tabPanel('Input Sequences',
    sidebarLayout(
      sidebarPanel(width=4,
        fileInput('fasta','Select file of sequences',multiple=FALSE,width='100%',buttonLabel='Select...'),
        radioButtons('align_choice','Align Sequences',choices=c('Yes','No'))
        ),
      mainPanel(
        tableOutput('aligned_seqs')
        )
      )
    ),
  
  tabPanel('Alleles to Use',
    sidebarLayout(
      sidebarPanel(width=4,
        radioButtons('which_allele','Select Allele File',choices=c('Use Precomputed Allele Frequencies','Enter a custom set of Frequencies'),selected=NULL),
        fileInput('hlafreq','Select HLA frequency File',multiple=FALSE,width='100%',buttonLabel = 'Select',placeholder='Choose File')
        ),
      mainPanel(
        tableOutput('hla')
        )
      )
    ),
   
  tabPanel('Individual Heatmaps',
    sidebarLayout(
      sidebarPanel(width=4,
        selectInput('peptide_heatmap','Peptide','raw_seqs'),
        sliderInput('aff_thresh_heatmap','Affinity Threshold',0,100,10)
        ),
      mainPanel(
        plotOutput('heatmap',height = 800)
        )
      )
    ),
  
  tabPanel('Individual Promiscuity Plots',
    sidebarLayout(
      sidebarPanel(width=4,
        selectInput('peptide_promiscuity','Peptide','raw_seqs'),
        sliderInput('aff_thresh_prom','Affinity Threshold',0,100,10),
        checkboxGroupInput('Population_prom','Choose Population of Interest')
        ),
      mainPanel(
        plotOutput('prom_plot',height = 800)
        )
      )
    ),
  
  tabPanel('Multiple Promiscuity Plots',
    sidebarLayout(
      sidebarPanel(width=4,
        sliderInput('aff_thresh_multi_prom','Affinity Threshold',0,100,10),
        selectInput('Population_multi_prom','Choose Population of Interest',choices='Select Population'),
        checkboxGroupInput('peptide_multi_promiscuity','Peptide(s) of Interest')
        ),
      mainPanel(
        plotOutput('multi_prom_plot',height = 800)
        )
      )
    ),
  tabPanel('Generate Report',
    sidebarLayout(
      sidebarPanel(width=4,
      textInput('report_name','Enter a file name with no extensions','value'='Report'),
      checkboxGroupInput("peptide_report","Peptide(s) of interest"),
      selectInput('Population_report','Select Population of Interest',''),
      sliderInput('aff_thresh_report','Affinity Threshold',0,100,10),
      downloadButton('generate_html_report','Generate HTML Report')
      #downloadButton('generate_pdf_report','Generate PDF Report')
      ),
	mainPanel()	  
    )
  )
)
)
server <- function(input, output,session) {

#pool <- dbPool(
#  drv = RMariaDB::MariaDB(),
#  dbname = "shiny_affinities",
#  host = "shiny-aff-1.c38jawccscp1.us-east-1.rds.amazonaws.com",
#  username = "admin_shiny_aff",
#  password = "5Yn2*$pxqRnK#p%8"
#)

observe('input_format')
  
  file_format='html'
infile='./Nworld.csv'
  
output$aligned_seqs<-renderTable({
  fasta_in<-input$fasta
  achoice<-input$align_choice
  validate(
    need(input$fasta != "",'Choose Fasta sequence on the left')
    )
  aas<-readAAStringSet(fasta_in$datapath)
  if (achoice=='Yes'){
    m_align<<-muscle(aas)@unmasked
    m_align_char<<-as.character(m_align)
    } else {
      m_align<<-as.character(aas)
      }
  aligned_seq_show<<-data.frame("Title"=names(m_align_char),'Aligned Sequence'=m_align_char,'Un-Aligned Sequence'=as.character(aas),stringsAsFactors = FALSE)
  print(aligned_seq_show)
  updateSelectInput(session,'peptide_heatmap',choices=names(m_align_char))
  updateSelectInput(session,'peptide_promiscuity',choices=names(m_align_char))
  updateCheckboxGroupInput(session,'peptide_multi_promiscuity',choices=names(m_align),selected = names(m_align))
  updateCheckboxGroupInput(session,'peptide_report',choices=names(m_align),selected = names(m_align))
  for (i in 1:dim(aligned_seq_show)[1]){
    seq=as.character(aligned_seq_show[i,3])
    
    name=as.character(aligned_seq_show[i,1])
    eval(parse(text=paste0(name,'<<-new("peptide",sequence=seq,affinities_matrix=make_matrix(seq,"Nworld.csv",pool))')))
    eval(parse(text=paste0(name,'@blanks<<-which(',name,'@sequence=="-")')))
    }
  if(dim(aligned_seq_show)[1]>1){
    aligned_seq_show<-data.frame("Title"=c(names(m_align_char),'Consensus'),'Aligned Sequence'=c(m_align_char,consensusString(m_align)),'Un-Aligned Sequence'=c(as.character(aas),consensusString(m_align)),stringsAsFactors = FALSE)
    rownames(aligned_seq_show)[(length(rownames(aligned_seq_show)))]<-'Consensus'
    
    }
 
  for (name in aligned_seq_show$Title){
    if (name != 'Consensus')  {
      eval(parse(text=paste0(name,"@blanks<<-gregexpr(pattern='-',aligned_seq_show[aligned_seq_show$Title=='",name,"',]$Aligned.Sequence,fixed=TRUE)[[1]]")))
     
      }
  }
  colnames(aligned_seq_show)<-c('','Aligned Sequence','Sequence')
 return(xtable(aligned_seq_show[,c(1,3,2)],auto=TRUE))
})

output$hla<-renderTable({
  
  
  if(input$which_allele=='Use Precomputed Allele Frequencies'){
  hla_in<-read.csv('Nworld.csv')
  } else {
    hla_in<-read.csv(input$hlafreq$datapath,sep=',')
  }
  population_names<-colnames(hla_in)
  updateCheckboxGroupInput(session,'Population_prom',choices=population_names[2:length(population_names)],selected='World')
  updateSelectInput(session,'Population_multi_prom',choices=population_names[2:length(population_names)],selected='World')
  updateSelectInput(session,'Population_report',choices=population_names[2:length(population_names)],selected='World')
  hla_in<<-hla_in
  print(hla_in)
  return(hla_in)
})



output$heatmap<-renderPlot({
  mat<-eval(parse(text=paste0(input$peptide_heatmap,'@affinities_matrix')))
  return(plot_matrix_threshold(mat,input$aff_thresh_heatmap,'bottom'))
})

output$prom_plot<-renderPlot({
  if(is.null(input$hlafreq)){
    allele_frequencies<-read.csv('Nworld.csv')
  } else {
    allele_frequencies<-read.csv(input$hlafreq$datapath,sep=',')
  }
  #allele_frequencies<-read.csv(input$hlafreq$datapath)
  allele_frequencies<-as.matrix(allele_frequencies[,colnames(allele_frequencies) %in% input$Population_prom])
  mat<-eval(parse(text=paste0(input$peptide_promiscuity,'@affinities_matrix')))
  mat<-mat<input$aff_thresh_prom
  mat<-matrix(as.numeric(mat),nrow=nrow(mat),byrow=FALSE)
  promlines=mat%*%(allele_frequencies/100)
  show_legend_logical<-length(input$Population)>1
  axis_seq<-unlist(strsplit(eval(parse(text=paste0(input$peptide_promiscuity,'@sequence'))),''))
  axis_seq<-axis_seq[8:(length(axis_seq)-7)]
  return(plot_single_promiscuity(promlines,show_legend_logical,axis_seq))
  })

output$multi_prom_plot<-renderPlot({
  peptides_for_multiplot<-input$peptide_multi_promiscuity
  
  if(is.null(input$hlafreq)){
    allele_frequencies<-read.csv('Nworld.csv',row.names = 1)
  } else {
    allele_frequencies<-read.csv(input$hlafreq$datapath,sep=',')
  }
  allele_frequencies<-as.matrix(allele_frequencies[,colnames(allele_frequencies) %in% input$Population_multi_prom])
  
  multi_promiscuity_melted<-data.frame('Peptide'=NA,'Position'=NA,'Promiscuity'=NA)
  
  for (name in aligned_seq_show$Title){
    
    if (name %in% peptides_for_multiplot){
      if (name != 'Consensus')  {
        mat<-eval(parse(text=paste0(name,'@affinities_matrix')))
        mat<-mat<input$aff_thresh_multi_prom
        mat<-matrix(as.numeric(mat),nrow=nrow(mat),byrow=FALSE)
        
        promlines=unlist(mat%*%as.matrix(allele_frequencies/100))
        #eval(parse(text=paste0(name,'@blanks')))
        if (eval(parse(text=paste0('max(',name,'@blanks)>0')))){
          for (i in eval(parse(text=paste0(name,'@blanks')))) {
            promlines<-append(promlines,NA,after=i-1)
            }
        }
        temp_frame=data.frame('Peptide'=rep(name,(nchar(m_align[1]))),'Position'=1:(nchar(m_align[1])),'Promiscuity'=c(rep(NA,7),promlines[1:(nchar(m_align[1])-14)],rep(NA,7)))
        multi_promiscuity_melted<-rbind(multi_promiscuity_melted,temp_frame)
        }
      }
    }
  multi_promiscuity_melted<-multi_promiscuity_melted[-is.na(multi_promiscuity_melted$Peptide),]
  return(plot_multi_promiscuity(multi_promiscuity_melted))
  })

output$generate_html_report<-downloadHandler(
  filename=paste0(input$report_name,'.html'),
  content = function(file) {
    tempReport <- file.path(tempdir(), 'Report_html.rmd')
    file.copy(paste0('Report_html.rmd'), tempReport, overwrite = TRUE)
    
    # Set up parameters to pass to Rmd document
    params <- list(report_peptides=input$peptide_report,
                   thresh=input$aff_thresh_report,
                   aligntable=aligned_seq_show[input$peptide_report,c(1,3,2)],
                   
                    hlatable=hla_in[,input$Population_report],
                   population=input$Population_report,
                   output_format='html')
    
    # Knit the document, passing in the `params` list, and eval it in a
    # child of the global environment (this isolates the code in the document
    # from the code in this app).
    rmarkdown::render(tempReport, output_file = file,
                      params = params,
                      envir = new.env(parent = globalenv()),
                      output_format = 'html_document')
    
  }
  )
  
output$generate_pdf_report<-downloadHandler(
  filename=paste0(input$report_name,'.pdf'),
  content = function(file) {
    
    tempReport <- file.path(tempdir(), paste0('Report_pdf.rmd'))
    file.copy(paste0('Report_pdf.rmd'), tempReport, overwrite = TRUE)
    
    # Set up parameters to pass to Rmd document
    params <- list(report_peptides=input$peptide_report,
                   thresh=input$aff_thresh_report,
                   aligntable=aligned_seq_show[input$peptide_report,c(1,3,2)],
                   
                   hlatable=hla_in[,input$Population_report],
                   population=input$Population_report,
                   output_format='pdf')
    
    # Knit the document, passing in the `params` list, and eval it in a
    # child of the global environment (this isolates the code in the document
    # from the code in this app).
    rmarkdown::render(tempReport, output_file = file,
                      params = params,
                      envir = new.env(parent = globalenv()),
                      output_format = 'pdf_document')
    
  }
)

}#server end bracket

shinyApp(ui = ui, server = server)

