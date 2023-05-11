 #' @export
 #'


 shinyAppsDB <- function(GAlist=NULL, what="DrugDBTables") {

  require(shiny)
  require(data.table)
  require(DT)
  require(ggplot2)
  require(httr)
  require(jsonlite)
  #require(stringr)

  if(what=="customGSEA") {

   # Define the UI
   ui <- fluidPage(
    titlePanel("Gene Set Enrichment Analysis of Gene Complexes"),
    sidebarLayout(
     sidebarPanel(
      selectInput("studyID", "Select Study:",
                  choices = GAlist$study$id,
                  selected="GWAS1"),
      textAreaInput("symbols", "Enter Gene Symbols (one per line or separated by spaces):"),
      numericInput("threshold", "Enter Interaction Threshold (0-999):", value = 900),
      numericInput("pval", "Enter P-value Threshold Enrichment Test:", value = 0.0001),
      br(),
      actionButton("submit", "Submit")
     ),
     mainPanel(
      tabsetPanel(
       tabPanel("Enrichment", DTOutput("hgtest")),
       tabPanel("Interactions", DTOutput("interactions")),
       tabPanel("Ensembl Gene IDs", DTOutput("features")),
       tabPanel("Summary Statistics", DTOutput("markers")),
       tabPanel("Drug Targets", DTOutput("targets")),
       tabPanel("Manhatten Plot", plotOutput("manhattenplot")),
       tabPanel("Studies", DTOutput("studies"), selected = TRUE)
      )
     )
    )
   )

   # Define the server logic
   server <- function(input, output, session) {

    # Prepare data frame for study information
    df_studies <- as.data.frame(GAlist$study)[,-c(2,11)]
    df_studies$neff <- as.integer(df_studies$neff)
    df_studies$reference <- createURL(url="https://pubmed.ncbi.nlm.nih.gov/",
                                      urlid=gsub("PMID:","",df_studies$reference))

    output$studies <- renderDT(server=FALSE,{
     createDT(df_studies)
    })

    # Get marker data if new study selected
    stat <- reactive({
     getMarkerStatDB(GAlist = GAlist, studyID = input$studyID)
    })


    # Display table when submit button is clicked
    observeEvent(input$submit, {

     # Check input symbols
     symbols <- gsub("[^[:alnum:]]", " ", input$symbols)
     symbols <- strsplit(symbols, "\\s+")[[1]]

     # Prepare data frame for gene symbol information
     df_features <- as.data.frame(t(sapply(symbols, function(x){getGeneDB(x)})))
     colnames(df_features) <- c("Symbol", "Ensembl Gene Id", "Description", "Type", "Chr", "Start", "End")

     interactions <- getInteractionsDB(ids=symbols, species="9606", threshold=input$threshold)
     interactions[,1] <- gsub("9606.","",interactions[,1], fixed=TRUE)
     interactions[,2] <- gsub("9606.","",interactions[,2], fixed=TRUE)
     colnames(interactions)[1:5] <- c("Prot1","Prot2","Sym1","Sym2","TaxID")
     ensp <- unique(c(interactions[,1],interactions[,2]))

     # Extract gene-marker sets
     marker_sets_indices <- getMarkerSetsDB(GAlist = GAlist, feature = "Genes")
     marker_sets_indices <- qgg:::mapSets(sets=marker_sets_indices, rsids=stat()$rsids, index=TRUE)
     pgenes <- sapply(marker_sets_indices, function(x) {min(stat()$p[x])})


     # Identify feature genes
     feature_genes <- unique(na.omit(unlist(GAlist$gsets$ensp2ensg[ensp])))
     feature_genes <- unique(c(feature_genes,na.omit(df_features[,2])))
     feature_genes <- feature_genes[feature_genes%in%names(pgenes)]
     feature_markers <- unique(unlist(marker_sets_indices[feature_genes]))

     # Get selected genes based on p-value threshold
     selected_genes <- names(pgenes)[pgenes < input$pval]

     # Calculate number of genes in each feature and number of associated genes based on p-value threshold
     number_genes_feature <- length(feature_genes)
     number_associated_genes_feature <- sum(feature_genes%in%selected_genes)

     # Calculate hypergeometric test p-values
     phgt <- hgtestDB(p = pgenes, sets = list(feature_genes), threshold = input$pval)

     # Calculate enrichment factor
     ef <- (number_associated_genes_feature/number_genes_feature)/
      (length(selected_genes)/length(pgenes))

     # Prepare data frame for hg test results
     df_hg <- data.frame(feature = "PPI",
                         ng = number_genes_feature,
                         nag = number_associated_genes_feature,
                         ef=ef,
                         phgt = phgt)
     colnames(df_hg) <- c("Feature", "Number of Genes",
                          "Number of Associated Genes",
                          "Enrichment Factor",
                          "P-value")

     # Render feature table
     output$hgtest <- renderDT(server=FALSE,{
      createDT(df_hg)%>%
       formatRound("Enrichment Factor", digits=2)%>%
       formatSignif('P-value',3)
     })

     # Render protein interaction data
     output$interactions <- renderDT(server=FALSE,{
      createDT(interactions)
     })

     # Prepare data frame for marker information
     df_markers <- stat()[feature_markers,]
     rsids <- df_markers$rsids
     df_markers$rsids <- createURL(url="https://www.ncbi.nlm.nih.gov/snp/",
                                   urlid=rsids)
     df_markers$gwas <- createURL(url="https://www.ebi.ac.uk/gwas/variants/",
                                  urlid=rsids)

     # Prepare data frame for feature information
     urlid <- df_features[,2]
     df_features[,2] <- createURL(url="https://www.ensembl.org/Homo_sapiens/Gene/Summary?g=",
                                  urlid=urlid)
     df_features$open_target_gene <- createURL(url="https://genetics.opentargets.org/gene/",
                                               urlid=urlid)
     df_features$open_target_target <- createURL(url="https://platform.opentargets.org/target/",
                                                 urlid=urlid)

     # Prepare data frame for feature drug information
     df_targets <- GAlist$targets
     df_targets$ATC <- createURL(url="https://www.whocc.no/atc_ddd_index/?code=",
                                 urlid=df_targets$ATC)
     df_feature_targets <- NULL
     has_targets <- df_targets$Target%in%feature_genes
     if(any(has_targets)) {
      df_feature_targets <- data.frame(df_targets[df_targets$Target%in%feature_genes,],P=NA)
      rws <- match(df_feature_targets$Target, names(pgenes))
      df_feature_targets$P[!is.na(rws)] <- pgenes[rws[!is.na(rws)]]
     }

     # Render marker table
     output$markers <- renderDT(server=FALSE,{
      createDT(df_markers)%>%
       formatRound("n", digits=0)
     })

     # Render feature table
     output$features <- renderDT(server=FALSE,{
      createDT(df_features)
     })

     # Render targets table
     output$targets <- renderDT(server=FALSE,{
      createDT(df_feature_targets)
     })

     # Create qqplot
     output$qqplot <- renderPlot({
      qqplotDB(p=pgenes, main="Markers")
     })

     # Create manhatten plot
     output$manhattenplot <- renderPlot({
      mhplotDB(p=df_markers$p, main="Markers in selected features")
     })

    })
   }

  }

  if(what=="geneDrugs") {
   # Define the UI
   ui <- fluidPage(
    titlePanel("Drug Target Analysis of Genomic Features"),
    sidebarLayout(
     sidebarPanel(
      selectInput("studyID", "Select Study:",
                  choices = GAlist$study$id,
                  selected="GWAS1"),
      textAreaInput("symbols", "Enter Gene Symbols (one per line or separated by spaces):"),
      br(),
      actionButton("submit", "Submit")
     ),
     mainPanel(
      tabsetPanel(
       tabPanel("Ensembl Gene IDs", DTOutput("features")),
       tabPanel("Summary Statistics", DTOutput("markers")),
       tabPanel("Drug Targets", DTOutput("targets")),
       tabPanel("Manhatten Plot", plotOutput("manhattenplot")),
       tabPanel("Studies", DTOutput("studies"), selected = TRUE)
      )
     )
    )
   )

   # Define the server logic
   server <- function(input, output, session) {

    # # Connect to the Ensembl database
    # ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
    #
    # # Define the function to translate gene symbols to Ensembl IDs
    # symbol_to_id <- function(symbols) {
    #  query <- getBM(attributes=c("ensembl_gene_id", "external_gene_name"),
    #                 filters="external_gene_name",
    #                 values=symbols,
    #                 mart=ensembl)
    #  return(query)
    # }
    #
    # # Call the function when the user inputs a list of gene symbols
    # output$output <- renderTable({
    #  if (nchar(input$symbols) > 0) {
    #   symbol_list <- unlist(strsplit(input$symbols, "[\\s\n]+"))
    #   result <- symbol_to_id(symbol_list)
    #   if (nrow(result) == 0) {
    #    "No results found."
    #   } else {
    #    result
    #   }
    # }
    # })

    # Define the function to translate gene symbols to Ensembl IDs
    symbol_to_id <- function(symbol) {
     base_url <- "https://rest.ensembl.org"
     url <- paste0(base_url, "/lookup/symbol/homo_sapiens/", symbol)
     response <- GET(url = url, content_type("application/json"))
     gene <- fromJSON(content(response, "text"), flatten = TRUE)
     if(!is.null(gene$error)) gene_information <- c(symbol, rep(NA, 6))
     if(is.null(gene$error)) gene_information <- c(symbol, gene$id, gene$description, gene$biotype, gene$seq_region_name, gene$start, gene$end)
     return(gene_information)
    }


    # Prepare data frame for study information
    output$studies <- renderDT(server=FALSE,{
     studies <- as.data.frame(GAlist$study)[,-c(2,11)]
     studies$neff <- as.integer(studies$neff)
     pmid <- gsub("PMID:","",studies$reference)
     pmid_link <- paste0("https://pubmed.ncbi.nlm.nih.gov/", pmid, "/")
     html_code <- paste0("<a href='", pmid_link, "' target='_blank'>", pmid, "</a>")
     studies$reference <- html_code

     datatable(studies,
               extensions = "Buttons",
               options = list(paging = TRUE,
                              scrollX=TRUE,
                              searching = TRUE,
                              ordering = TRUE,
                              dom = 'Bfrtip',
                              buttons = c('csv', 'excel'),
                              pageLength=10,
                              lengthMenu=c(3,5,10) ),
               rownames= FALSE,
               escape = FALSE)
    })

    stat <- reactive({
     getMarkerStatDB(GAlist = GAlist, studyID = input$studyID)
    })


    # Display table when submit button is clicked
    observeEvent(input$submit, {

     symbol_list <- unlist(strsplit(input$symbols, "[\\s\n]+"))
     #result <- symbol_to_id(symbol_list)
     result <- as.data.frame(t(sapply(symbol_list, function(x){symbol_to_id(x)})))
     colnames(result) <- c("Symbol", "Ensembl Gene Id", "Description", "Type", "Chr", "Start", "End")

     feature_genes <- result[,2]

     message("Extract gene-marker sets")

     # Extract gene-marker sets
     marker_sets_indices <- getMarkerSetsDB(GAlist = GAlist, feature = "Genes")
     marker_sets_indices <- qgg:::mapSets(sets=marker_sets_indices, rsids=stat()$rsids, index=TRUE)
     pgenes <- sapply(marker_sets_indices, function(x) {min(stat()$p[x])})

     # Create data frame for markers
     selected_markers <- unique(unlist(marker_sets_indices[feature_genes]))

     df_markers <- stat()[selected_markers,]
     rsid <- df_markers$rsids
     gwas_url <- paste0("https://www.ebi.ac.uk/gwas/variants/", rsid)
     dbsnp_url <- paste0("https://www.ncbi.nlm.nih.gov/snp/", rsid)
     html_code <- paste0("<a href='", dbsnp_url, "' target='_blank'>", rsid, "</a>")
     df_markers$rsids <- html_code

     html_code <- paste0("<a href='", gwas_url, "' target='_blank'>", rsid, "</a>")
     df_markers$gwas <- html_code



     # Create data frame for features
     df_features <- result

     gene_id <- df_features[,2]
     gene_link <- paste0("https://www.ensembl.org/Homo_sapiens/Gene/Summary?g=", gene_id)
     html_code <- paste0("<a href='", gene_link, "' target='_blank'>", gene_id, "</a>")
     df_features[,2] <- html_code

     ot_gene_link <- paste0("https://genetics.opentargets.org/gene/", gene_id)
     html_code <- paste0("<a href='", ot_gene_link, "' target='_blank'>", gene_id, "</a>")
     df_features$open_target_gene <- html_code
     ot_target_link <- paste0("https://platform.opentargets.org/target/", gene_id)
     html_code <- paste0("<a href='", ot_target_link, "' target='_blank'>", gene_id, "</a>")
     df_features$open_target_target <- html_code


     # Create data frame for feature drug targets
     df_targets <- GAlist$targets
     atc_code <- df_targets$ATC
     atc_url <- paste0("https://www.whocc.no/atc_ddd_index/?code=", atc_code)
     html_code <- paste0("<a href='", atc_url, "' target='_blank'>", atc_code, "</a>")
     df_targets$ATC <- html_code
     df_feature_targets <- NULL
     has_targets <- df_targets$Target%in%feature_genes
     if(any(has_targets)) {
      df_feature_targets <- data.frame(df_targets[df_targets$Target%in%feature_genes,],P=NA)
      rws <- match(df_feature_targets$Target, names(pgenes))
      df_feature_targets$P[!is.na(rws)] <- pgenes[rws[!is.na(rws)]]
     }

     # Render marker table
     output$markers <- renderDT(server=FALSE,{
      datatable(df_markers,
                escape = FALSE,
                extensions = "Buttons",
                options = list(paging = TRUE,
                               scrollX=TRUE,
                               searching = TRUE,
                               ordering = TRUE,
                               dom = 'Bfrtip',
                               buttons = c('csv', 'excel'),
                               pageLength=10,
                               lengthMenu=c(3,5,10) ),
                rownames= FALSE)%>%
       formatRound("n", digits=0)
     })

     # Render feature table
     output$features <- renderDT(server=FALSE,{
      datatable(df_features, extensions = "Buttons",
                escape = FALSE,
                options = list(paging = TRUE,
                               scrollX=TRUE,
                               searching = TRUE,
                               ordering = TRUE,
                               dom = 'Bfrtip',
                               buttons = c('csv', 'excel'),
                               pageLength=10,
                               lengthMenu=c(3,5,10) ),
                rownames= FALSE)
     })

     # Render target table
     output$targets <- renderDT(server=FALSE,{
      datatable(df_feature_targets,
                extensions = "Buttons",
                escape = FALSE,
                options = list(paging = TRUE,
                               scrollX=TRUE,
                               searching = TRUE,
                               ordering = TRUE,
                               dom = 'Bfrtip',
                               buttons = c('csv', 'excel'),
                               pageLength=10,
                               lengthMenu=c(3,5,10) ),
                rownames= FALSE)
     })

     # Create qqplot
     output$qqplot <- renderPlot({
      qqplotDB(p=pgenes, main="Markers")
     })


     # Create manhatten plot
     output$manhattenplot <- renderPlot({
      ggplot(data = stat()[selected_markers,],
             aes(x = 1:length(selected_markers),
                 y = -log10(p))) +
       geom_point(size = 2) +
       geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
       xlab("Markers") +
       ylab("-log10(p-value)") +
       ggtitle("Association of Markers with Disease") +
       theme_minimal()
     })
    })
   }

  }

  if(what=="featureDRUGS") {

   # Define UI
   ui <- fluidPage(
    titlePanel("Drug Target Analysis of Genomic Features"),
    sidebarLayout(
     sidebarPanel(
      selectInput("studyID", "Select Study:",
                  choices = GAlist$study$id,
                  selected="GWAS1"),
      selectInput("feature", "Select Feature:", choices = c("Genes",
                                                            "Pathways",
                                                            "GO",
                                                            "DrugGenes",
                                                            "ProteinComplexes",
                                                            "ChemicalComplexes")),
      textInput("featureID", "Select Feature ID:", value = "R-HSA-202427"),
      br(),
      actionButton("submit", "Submit")
     ),
     mainPanel(
      tabsetPanel(
       tabPanel("Ensembl Gene IDs", DTOutput("features")),
       tabPanel("Summary Statistics", DTOutput("markers")),
       tabPanel("Drug Targets", DTOutput("targets")),
       tabPanel("Manhatten Plot", plotOutput("manhattenplot")),
       tabPanel("Studies", DTOutput("studies"), selected = TRUE)
      )
      #h4("Selected Inputs"),
      #verbatimTextOutput("selected_inputs_study"),
      #verbatimTextOutput("selected_inputs_feature")
     )
    )
   )


   # Define server
   server <- function(input, output, session) {


    # Prepare data frame for study information
    output$studies <- renderDT(server=FALSE,{
     studies <- as.data.frame(GAlist$study)[,-c(2,11)]
     studies$neff <- as.integer(studies$neff)
     pmid <- gsub("PMID:","",studies$reference)
     pmid_link <- paste0("https://pubmed.ncbi.nlm.nih.gov/", pmid, "/")
     html_code <- paste0("<a href='", pmid_link, "' target='_blank'>", pmid, "</a>")
     studies$reference <- html_code

     datatable(studies,
               extensions = "Buttons",
               options = list(paging = TRUE,
                              scrollX=TRUE,
                              searching = TRUE,
                              ordering = TRUE,
                              dom = 'Bfrtip',
                              buttons = c('csv', 'excel'),
                              pageLength=10,
                              lengthMenu=c(3,5,10) ),
               rownames= FALSE,
               escape = FALSE)
    })

    observeEvent(input$studyID, {
     # Perform actions when the studyID is changed
     message("StudyID event observed")
     output$selected_inputs_study <- renderText({
      paste0("Selected Study: ", input$studyID)
     })
    })

    stat <- reactive({
     getMarkerStatDB(GAlist = GAlist, studyID = input$studyID)
    })


    observeEvent(input$feature, {
     output$selected_inputs_feature <- renderText({
      paste0("Selected Feature: ", input$feature)
     })

     message("New feature selected")
     if(input$feature=="Genes") feature_value <- "ENSG00000149257"
     if(input$feature=="Pathways") feature_value <- "R-HSA-202427"
     if(input$feature=="GO") feature_value <- "GO:0000002"
     if(input$feature=="DrugGenes") feature_value <- "REPAGLINIDE"
     if(input$feature=="ProteinComplexes") feature_value <- "ENSP00000000233"
     if(input$feature=="ChemicalComplexes") feature_value <- "CIDm00000001"
     updateTextInput(session, "featureID", value = feature_value )
    })


    # Display table when submit button is clicked
    observeEvent(input$submit, {

     message("Extract gene-marker sets")

     # Extract gene-marker sets
     marker_sets_indices <- getMarkerSetsDB(GAlist = GAlist, feature = "Genes")
     marker_sets_indices <- qgg:::mapSets(sets=marker_sets_indices, rsids=stat()$rsids, index=TRUE)
     pgenes <- sapply(marker_sets_indices, function(x) {min(stat()$p[x])})

     if (input$feature=="Genes") feature_genes <- input$featureID

     message("Extract feature-gene sets")

     # Extract feature-gene sets
     if(!input$feature=="Genes") {
      feature <- input$feature
      featureID <- input$featureID
      #if(feature=="ProteinComplexes") feature <- "ProteinComplexes2Genes"
      #if(feature=="ChemicalComplexes") feature <- "ChemicalComplexes2Genes"
      feature_genes <- getSetsDB(GAlist = GAlist, feature = feature)
      feature_genes <- qgg:::mapSets(sets=feature_genes, rsids=names(marker_sets_indices), index=FALSE)
      feature_genes <- unlist(feature_genes[featureID])
     }


     # Create data frame for markers
     selected_markers <- unique(unlist(marker_sets_indices[feature_genes]))

     df_markers <- stat()[selected_markers,]
     rsid <- df_markers$rsids
     gwas_url <- paste0("https://www.ebi.ac.uk/gwas/variants/", rsid)
     dbsnp_url <- paste0("https://www.ncbi.nlm.nih.gov/snp/", rsid)
     html_code <- paste0("<a href='", dbsnp_url, "' target='_blank'>", rsid, "</a>")
     df_markers$rsids <- html_code

     html_code <- paste0("<a href='", gwas_url, "' target='_blank'>", rsid, "</a>")
     df_markers$gwas <- html_code



     # Create data frame for features
     df_features <- data.frame(gene_id=feature_genes)

     gene_id <- df_features$gene_id
     gene_link <- paste0("https://www.ensembl.org/Homo_sapiens/Gene/Summary?g=", gene_id)
     html_code <- paste0("<a href='", gene_link, "' target='_blank'>", gene_id, "</a>")
     df_features$gene_id <- html_code

     ot_gene_link <- paste0("https://genetics.opentargets.org/gene/", gene_id)
     html_code <- paste0("<a href='", ot_gene_link, "' target='_blank'>", gene_id, "</a>")
     df_features$open_target_gene <- html_code
     ot_target_link <- paste0("https://platform.opentargets.org/target/", gene_id)
     html_code <- paste0("<a href='", ot_target_link, "' target='_blank'>", gene_id, "</a>")
     df_features$open_target_target <- html_code


     # Create data frame for feature drug targets
     df_targets <- GAlist$targets
     atc_code <- df_targets$ATC
     atc_url <- paste0("https://www.whocc.no/atc_ddd_index/?code=", atc_code)
     html_code <- paste0("<a href='", atc_url, "' target='_blank'>", atc_code, "</a>")
     df_targets$ATC <- html_code
     df_feature_targets <- NULL
     has_targets <- df_targets$Target%in%feature_genes
     if(any(has_targets)) {
      df_feature_targets <- data.frame(df_targets[df_targets$Target%in%feature_genes,],P=NA)
      rws <- match(df_feature_targets$Target, names(pgenes))
      df_feature_targets$P[!is.na(rws)] <- pgenes[rws[!is.na(rws)]]
     }

     # Render marker table
     output$markers <- renderDT(server=FALSE,{
      datatable(df_markers,
                escape = FALSE,
                extensions = "Buttons",
                options = list(paging = TRUE,
                               scrollX=TRUE,
                               searching = TRUE,
                               ordering = TRUE,
                               dom = 'Bfrtip',
                               buttons = c('csv', 'excel'),
                               pageLength=10,
                               lengthMenu=c(3,5,10) ),
                rownames= FALSE)%>%
       formatRound("n", digits=0)
     })

     # Render feature table
     output$features <- renderDT(server=FALSE,{
      datatable(df_features, extensions = "Buttons",
                escape = FALSE,
                options = list(paging = TRUE,
                               scrollX=TRUE,
                               searching = TRUE,
                               ordering = TRUE,
                               dom = 'Bfrtip',
                               buttons = c('csv', 'excel'),
                               pageLength=10,
                               lengthMenu=c(3,5,10) ),
                rownames= FALSE)
     })

     # Render target table
     output$targets <- renderDT(server=FALSE,{
      datatable(df_feature_targets,
                extensions = "Buttons",
                escape = FALSE,
                options = list(paging = TRUE,
                               scrollX=TRUE,
                               searching = TRUE,
                               ordering = TRUE,
                               dom = 'Bfrtip',
                               buttons = c('csv', 'excel'),
                               pageLength=10,
                               lengthMenu=c(3,5,10) ),
                rownames= FALSE)
     })

     # Create qqplot
     output$qqplot <- renderPlot({
      qqplotDB(p=pgenes, main="Markers")
     })


     # Create manhatten plot
     output$manhattenplot <- renderPlot({
      ggplot(data = stat()[selected_markers,],
             aes(x = 1:length(selected_markers),
                 y = -log10(p))) +
       geom_point(size = 2) +
       geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
       xlab("Markers") +
       ylab("-log10(p-value)") +
       ggtitle("Association of Markers with Disease") +
       theme_minimal()
     })


    })

   }


  }

  if(what=="featureGSEA") {
   # Define UI
   ui <- fluidPage(
    titlePanel("Gene Set Enrichment Analysis of Genomic Features"),
    sidebarLayout(
     sidebarPanel(
      selectInput("studyID", "Select Study:",
                  choices = GAlist$study$id,
                  selected="GWAS1"),
      selectInput("feature", "Select Feature:", choices = c("Pathways",
                                                            "GO",
                                                            "DrugGenes",
                                                            "ProteinComplexes",
                                                            "ChemicalComplexes")),
      numericInput("threshold", "Enter P-value Threshold:", value = 10e-5),
      br(),
      actionButton("submit", "Submit")
     ),
     mainPanel(
      tabsetPanel(
       tabPanel("Table", DTOutput("result")),
       tabPanel("GSEAplot", plotOutput("gseaplot")),
       tabPanel("QQplot", plotOutput("qqplot")),
       tabPanel("Studies", DTOutput("studies"), selected = TRUE)
      )
     )
    )
   )


   # Define server
   server <- function(input, output, session) {


    # Prepare data frame for study information
    output$studies <- renderDT(server=FALSE,{
     studies <- as.data.frame(GAlist$study)[,-c(2,11)]
     studies$neff <- as.integer(studies$neff)
     pmid <- gsub("PMID:","",studies$reference)
     pmid_link <- paste0("https://pubmed.ncbi.nlm.nih.gov/", pmid, "/")
     html_code <- paste0("<a href='", pmid_link, "' target='_blank'>", pmid, "</a>")
     studies$reference <- html_code

     datatable(studies,
               extensions = "Buttons",
               options = list(paging = TRUE,
                              scrollX=TRUE,
                              searching = TRUE,
                              ordering = TRUE,
                              dom = 'Bfrtip',
                              buttons = c('csv', 'excel'),
                              pageLength=10,
                              lengthMenu=c(3,5,10) ),
               rownames= FALSE,
               escape = FALSE)
    })

    # Perform actions when the studyID is changed
    observeEvent(input$studyID, {
     message("StudyID event observed")
     output$selected_inputs_study <- renderText({
      paste0("Selected Study: ", input$studyID)
     })
    })

    stat <- reactive({
     getMarkerStatDB(GAlist = GAlist, studyID = input$studyID)
    })


    # Display table when submit button is clicked
    observeEvent(input$submit, {

     message("Extract gene-marker sets")

     # Extract gene-marker sets
     marker_sets_indices <- getMarkerSetsDB(GAlist = GAlist, feature = "Genes")
     marker_sets_indices <- qgg:::mapSets(sets=marker_sets_indices, rsids=stat()$rsids, index=TRUE)
     pgenes <- sapply(marker_sets_indices, function(x) {min(stat()$p[x])})

     message("Extract feature-gene sets")

     # Extract feature-gene sets
     feature_genes <- getSetsDB(GAlist = GAlist, feature = input$feature)
     # if(input$feature=="ProteinComplexes") {
     #  min_combined_score <- 900
     #  min_interactions <- 5
     #
     #  file_string <- file.path(GAlist$dirs["gsets"],"9606.protein.links.v11.5.txt.gz")
     #  string <- fread(file_string, data.table=FALSE)
     #  string  <- string[string$combined_score>=min_combined_score,]
     #  string <- string[string$protein2%in%names(GAlist$gsets$ensp2ensg),]
     #  string <- split( string$protein2,f=as.factor(string$protein1))
     #  feature_genes <- string[sapply(string ,length)>=min_interactions]
     #  feature_genes <- lapply(feature_genes,function(x){na.omit(unlist(GAlist$gsets$ensp2ensg[x]))})
     # }
     # if(input$feature=="ChemicalComplexes") {
     #  min_combined_score <- 700
     #  min_interactions <- 5
     #  file_stitch <- "http://stitch.embl.de/download/protein_chemical.links.v5.0/9606.protein_chemical.links.v5.0.tsv.gz"
     #  stitch <- fread(file_stitch, data.table=FALSE)
     #  stitch$protein <- gsub("9606.","",stitch$protein)
     #  stitch <- stitch[stitch$protein%in%names(GAlist$gsets$ensp2ensg),]
     #  stitch  <- stitch[stitch$combined_score>=min_combined_score,]
     #  stitch <- split( stitch$protein,f=as.factor(stitch$chemical))
     #  feature_genes  <- stitch[sapply(stitch ,length)>=min_interactions]
     #  feature_genes <- lapply(feature_genes,function(x){na.omit(unlist(GAlist$gsets$ensp2ensg[x]))})
     # }

     feature_genes <- qgg:::mapSets(sets=feature_genes, rsids=names(pgenes), index=FALSE)

     # Get selected genes based on p-value threshold
     selected_genes <- names(pgenes)[pgenes < input$threshold]

     # Calculate number of genes in each feature and number of associated genes based on p-value threshold
     number_genes_feature <- sapply(feature_genes, function(x){length(x)})
     number_associated_genes_feature <- sapply(feature_genes, function(x){sum(x %in% selected_genes)})

     # Calculate hypergeometric test p-values
     phgt <- hgtestDB(p = pgenes, sets = feature_genes, threshold = input$threshold)

     # Calculate enrichment factor
     ef <- (number_associated_genes_feature[names(phgt)]/number_genes_feature[names(phgt)])/
      (length(selected_genes)/length(pgenes))

     # Create data frame for table
     df <- data.frame(feature = names(phgt),
                      ng = number_genes_feature[names(phgt)],
                      nag = number_associated_genes_feature[names(phgt)],
                      ef=ef,
                      phgt = phgt)
     colnames(df) <- c("Feature", "Number of Genes",
                       "Number of Associated Genes",
                       "Enrichment Factor",
                       "P-value")
     if(input$feature=="DrugGenes") {
      df$ATC <- rep("Unknown",nrow(df))
      has_atc <- match(tolower(df$Feature),tolower(GAlist$atc$name))
      df$ATC[!is.na(has_atc)] <- as.character(GAlist$atc$code[has_atc[!is.na(has_atc)]])
      atc_code <- df$ATC
      atc_url <- paste0("https://www.whocc.no/atc_ddd_index/?code=", atc_code)
      html_code <- paste0("<a href='", atc_url, "' target='_blank'>", atc_code, "</a>")
      df$ATC <- html_code
     }
     if(!input$feature=="DrugGenes") {
      feature_id <- df$Feature
      if(input$feature=="Pathways") feature_link <- paste0("https://reactome.org/content/detail/", feature_id)
      if(input$feature=="GO") feature_link <- paste0("https://www.ebi.ac.uk/QuickGO/term/", feature_id)
      if(input$feature=="ProteinComplexes") feature_link <- paste0("https://string-db.org/network/9606.", feature_id)
      if(input$feature=="ChemicalComplexes") feature_link <- paste0("https://pubchem.ncbi.nlm.nih.gov/compound/",
                                                                    substring(feature_id,5,nchar(as.character(feature_id))))
      html_code <- paste0("<a href='", feature_link, "' target='_blank'>", feature_id, "</a>")
      df$Feature <- html_code
     }

     # Render feature table
     output$result <- renderDT(server=FALSE,{
      datatable(df, extensions = "Buttons",
                escape = FALSE,
                options = list(paging = TRUE,
                               scrollX=TRUE,
                               searching = TRUE,
                               ordering = TRUE,
                               dom = 'Bfrtip',
                               buttons = c('csv', 'excel'),
                               pageLength=10,
                               lengthMenu=c(3,5,10) ),
                rownames= FALSE)%>%
       formatRound("Enrichment Factor", digits=2)%>%
       formatSignif('P-value',3)
     })

     # Create plot of -log10 p-values
     output$qqplot <- renderPlot({
      qqplotDB(p=phgt, main="Features")
     })

     # Create plot of -log10 p-values
     output$gseaplot <- renderPlot({
      plot(y=-log10(df[,5]),x=df[,4], frame.plot=FALSE,
           xlab="Enrichment factor", ylab="Enrichment -log10(P)")
     })

    })

   }

  }

  if(what=="DrugDBTables") {


   shinyApp(
    ui <- fluidPage(
     sidebarLayout(
      sidebarPanel(
       selectInput("what", "Select a file to load:",
                   choices = c("Drug Labels PharmGKB",
                               "Drug Labels Genes PharmGKB",
                               "Variants PharmGKB",
                               "Automated Annotations PharmGKB",
                               "Drugs DGI",
                               "Genes DGI",
                               "Categories DGI",
                               "Interactions DGI",
                               "Direct Associations OpenTarget",
                               "Direct Associations By Data Type OpenTarget",
                               "Direct Associations By Data Source OpenTarget"))
      ),
      mainPanel(
       DTOutput("data")
      )
     )
    )
    ,
    server <- function(input, output) {

     observeEvent(input$what, {
      if (input$what == "Drug Labels PharmGKB") {
       drugdb_file <- file.path(GAlist$dirs["pharmgkb"], "drugLabels.tsv")
      } else if (input$what == "Drug Labels Genes PharmGKB") {
       drugdb_file <- file.path(GAlist$dirs["pharmgkb"], "drugLabels.byGene.tsv")
      } else if (input$what == "Variants PharmGKB") {
       drugdb_file <- file.path(GAlist$dirs["pharmgkb"], "clinicalVariants.tsv")
      } else if (input$what == "Automated Annotations PharmGKB") {
       drugdb_file <- file.path(GAlist$dirs["pharmgkb"], "automated_annotations.tsv")
      } else if (input$what == "Drugs DGI") {
       drugdb_file <- file.path(GAlist$dirs["dgidb"], "drugs.tsv")
      } else if (input$what == "Genes DGI") {
       drugdb_file <- file.path(GAlist$dirs["dgidb"], "genes.tsv")
      } else if (input$what == "Categories DGI") {
       drugdb_file <- file.path(GAlist$dirs["dgidb"], "categories.tsv")
      } else if (input$what == "Interactions DGI") {
       drugdb_file <- file.path(GAlist$dirs["dgidb"], "interactions.tsv")
      } else if (input$what == "Direct Associations OpenTarget") {
       drugdb_file <- file.path(GAlist$dirs["opentargets"], "associationByOverallDirect.tsv")
      } else if (input$what == "Direct Associations By Data Type OpenTarget") {
       drugdb_file <- file.path(GAlist$dirs["opentargets"], "associationByDatatypeDirect.tsv")
      } else if (input$what == "Direct Associations By Data Source OpenTarget") {
       drugdb_file <- file.path(GAlist$dirs["opentargets"], "associationByDatasourceDirect.tsv")
      }
      df <- fread(drugdb_file, quote = "", data.table = FALSE)
      output$data <- renderDT(df,
                              options = list(pageLength = 50, dom = 'Bfrtip',
                                             buttons = c('copy', 'csv')),
                              extensions = "Buttons")
     })
    }

   )
  }

 #  library(shiny)
 #
 #  # Define UI
 #  ui <- fluidPage(
 #   titlePanel("Update and Process Genetic Association Study Summary Statistics"),
 #   sidebarLayout(
 #    sidebarPanel(
 #     selectInput("source", "Source of the Summary Statistics:",
 #                 choices = c("Source 1", "Source 2", "Source 3")),
 #     textInput("trait", "Trait Being Studied:"),
 #     selectInput("type", "Type of Trait:",
 #                 choices = c("Binary", "Quantitative")),
 #     textInput("gender", "Gender of Study Participants:"),
 #     textInput("ancestry", "Ancestry of Study Participants:"),
 #     textInput("build", "Genome Build Used in the Study:"),
 #     numericInput("n", "Sample Size:", value = 0, min = 0),
 #     numericInput("ncase", "Number of Cases:", value = 0, min = 0),
 #     numericInput("ncontrol", "Number of Controls:", value = 0, min = 0),
 #     textInput("reference", "Reference for the Study:"),
 #     textInput("comments", "Additional Comments:"),
 #     checkboxInput("writeStatDB", "Perform Quality Control and Write the Processed Summary Statistics to a File", value = TRUE),
 #     numericInput("excludeMAFDIFF", "Threshold to Exclude SNPs with a Difference in MAF Between Cases and Controls:", value = 0.05, min = 0, max = 1, step = 0.01),
 #     actionButton("run", "Run updateStatDB()")
 #    ),
 #    mainPanel(
 #     verbatimTextOutput("result")
 #    )
 #   )
 #  )
 #
 #  # Define server
 #  server <- function(input, output) {
 #
 #   # Define reactive object to store user inputs
 #   inputs <- reactive({
 #    list(
 #     GAlist = GAlist,
 #     stat = NULL,
 #     source = input$source,
 #     trait = input$trait,
 #     type = input$type,
 #     gender = input$gender,
 #     ancestry = input$ancestry,
 #     build = input$build,
 #     n = input$n,
 #     ncase = input$ncase,
 #     ncontrol = input$ncontrol,
 #     reference = input$reference,
 #     comments = input$comments,
 #     writeStatDB = input$writeStatDB,
 #     excludeMAFDIFF = input$excludeMAFDIFF
 #    )
 #   })
 #
 #   # Define action to run the updateStatDB() function
 #   observeEvent(input$run, {
 #    result <- tryCatch({
 #     updateStatDB(inputs())
 #    }, error = function(e) {
 #     paste("Error:", e$message)
 #    })
 #    output$result <- renderPrint(result)
 #   })
 #  }
 #
   # Run the app
   shinyApp(ui = ui, server = server)

 }
