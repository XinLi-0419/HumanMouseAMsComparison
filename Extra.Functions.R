# hGene/mGene homologue function
mouse_human_genes <- read.csv("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",sep = "\t")
mouse <- split.data.frame(mouse_human_genes,mouse_human_genes$Common.Organism.Name)[[2]]
human <- split.data.frame(mouse_human_genes,mouse_human_genes$Common.Organism.Name)[[1]]
mouse <- mouse[, c(1,4)]
human <- human[, c(1,4)]
mh_data <- merge.data.frame(mouse, human, by = "DB.Class.Key", all.y = TRUE)
colnames(mh_data) <- c("DB.Class.Key", "Symbol.Mouse", "Symbol.Human")
HtoM <- function(Human.Gene) {if (Human.Gene %in% mh_data$Symbol.Human) {mh_data$Symbol.Mouse[which(mh_data$Symbol.Human == Human.Gene)]} else {Human.Gene}}
MtoH <- function(Mouse.Gene) {if (Mouse.Gene %in% mh_data$Symbol.Mouse) {mh_data$Symbol.Human[which(mh_data$Symbol.Mouse == Mouse.Gene)]} else {Mouse.Gene}}
HtoM.melt <- function(Human.Gene) {Mouse.Gene <- c() 
for (i in 1:length(Human.Gene)) {Mouse.Gene <- c(Mouse.Gene, HtoM(Human.Gene[i]))} 
return(Mouse.Gene)}
MtoH.melt <- function(Mouse.Gene) {Human.Gene <- c() 
for (i in 1:length(Mouse.Gene)) {Human.Gene <- c(Human.Gene, MtoH(Mouse.Gene[i]))} 
return(Human.Gene)}

# topGOtable.targetGO # instead of return top GOs, return specific GOs
topGOtable.targetGO <- function(DEgenes,                  # Differentially expressed genes
                       BGgenes,                 # background genes
                       ontology = "BP",            # could use also "MF"
                       annot = annFUN.org,       # parameters for creating topGO object
                       mapping = "org.Mm.eg.db",
                       geneID = "symbol",       # could also beID = "entrez"
                       topTablerows = 200,
                       fullNamesInRows = TRUE,
                       addGeneToTerms = TRUE,
                       plotGraph = FALSE, 
                       plotNodes = 10,
                       writeOutput = FALSE, 
                       outputFile = "",
                       topGO_method2 = "elim",
                       TargetGO,
                       do_padj = FALSE) {
  # checking the additional topGO_method2
  topgo_methods <- c("elim", "weight", "weight01", "lea", "parentchild")
  if (!(topGO_method2 %in% topgo_methods))
    stop("Please provide one of the following topGO methods in addition to the classic method:\n",
         paste(topgo_methods, collapse = ", "))
  # creating the vectors
  DEgenes_input <- factor(as.integer(BGgenes %in% DEgenes))
  names(DEgenes_input) <- BGgenes
  # instantiating the topGOdata object
  GOdata <- new("topGOdata",
                ontology = ontology,
                allGenes = DEgenes_input,
                nodeSize = 10,
                annot = annot,
                mapping = mapping,
                ID = geneID)
  # performing the test(s)
  result_method2 <- runTest(GOdata, algorithm = topGO_method2, statistic = "fisher")
  resultClassic <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  sTab <- GenTable(GOdata,
                   p.value_method2 = result_method2,
                   p.value_classic = resultClassic,
                   orderBy = "p.value_method2",
                   ranksOf = "p.value_classic",
                   topNodes = length(score(resultClassic)))
  
  names(sTab)[which(names(sTab) == "p.value_method2")] <- paste0("p.value_", topGO_method2)
  
  sTab[["p.value_classic"]] <- as.numeric(sTab[["p.value_classic"]])
  sTab[[paste0("p.value_", topGO_method2)]] <- as.numeric(sTab[[paste0("p.value_", topGO_method2)]])
  
  # if FDR, then apply it here
  if (do_padj)
    sTab[[paste0("padj_BY_", topGO_method2)]] <- p.adjust(sTab[[paste0("p.value_", topGO_method2)]], method = "BY")
  
  # subset to specified number of rows
  sTab <- sTab[seq_len(nrow(sTab)), ]
  
  if (fullNamesInRows) {
    sTab$Term <- sapply(sTab$GO.ID, function(go) {
      Term(GOTERM[[go]])
    })
  }
  
  if (addGeneToTerms) {
    # adapted from an elegant one liner found here: https://support.bioconductor.org/p/65856/
    SignificantGenes <- sigGenes(GOdata)
    sTab$genes <- sapply(sTab$GO.ID, function(x) {
      genes <- genesInTerm(GOdata, x)
      tmp <- genes[[1]][genes[[1]] %in% SignificantGenes]
    })
    # coerce the list to a comma separated vector
    sTab$genes <- unlist(lapply(sTab$genes, function(arg) paste(arg, collapse = ",")))
  }
  
  sTab <- sTab[sTab$GO.ID %in% TargetGO, ]
  
  # write all entries of the table
  if (writeOutput) write.table(sTab, file = outputFile, sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
  if (plotGraph) showSigOfNodes(GOdata, topGO::score(result_method2), firstSigNodes = plotNodes, useInfo = "all")
  #   if(outputToLatex) sTabSig <- xtable(apply(sTabSig[1:15,], 2, as.character)) # take a smaller subset
  
  # and returns the significant ones # or all, like here
  return(sTab)
}
