#' Create the ASOspotter Shiny app
#' 
#' @param vcf A path to a VCF file.
#' @param bam A path to a BAM file.
#' @return A shiny app to interactively view variants.
#' 
#' @importFrom vcfR read.vcfR
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom GenomicAlignments readGAlignments 
#' @importFrom Rsamtools ScanBamParam
#' @importFrom igvShiny renderIgvShiny igvShinyOutput showGenomicRegion loadBamTrackFromLocalData 
#'   parseAndValidateGenomeSpec loadVcfTrack
#' @importFrom DT renderDT datatable
#' @importFrom VariantAnnotation readVcf
#' @importFrom shiny shinyApp addResourcePath observeEvent reactiveValues
#' 
#' @author Jared Andrews
#' @export
#' @examples
#' library(ASOspotter)
#' vcf <- system.file("extdata", "NA12878_HG001.hg38.benchmark.subset.vcf.gz", package = "ASOspotter")
#' bam <- system.file("extdata", "NA12878.sub.sorted.bam", package = "ASOspotter")
#' \dontrun{
#' ASOspotter(vcf, bam)
#' }
ASOspotter <- function(vcf, bam, genome = "hg38", track_dir = "./tracks") {

    # Make local dir for track subsets.
    if(!dir.exists(track_dir))
        dir.create(track_dir)
    addResourcePath("tracks", track_dir)

    # Load the VCF
    vcf.records <- read.vcfR(vcf)
    vcf.records <- cbind(as.data.frame(getFIX(vcf.records)), as.data.frame(vcf.records@gt))

    # Get initial locus
    initialLocus <- paste0(vcf.records[1, "CHROM"], ":", 
                           as.numeric(vcf.records[1, "POS"]) - 250, "-", as.numeric(vcf.records[1, "POS"]) + 250)

    # Set igvShiny options
    opts <- parseAndValidateGenomeSpec(genomeName = genome,  initialLocus = initialLocus)

    # Create the UI
    ui <- .create_ui()
    
    server <- function(input, output, session) {
        
        robjects <- reactiveValues(
            locus = initialLocus,
            locus_gr = GRanges(seqnames = vcf.records[1, "CHROM"], 
                               ranges = IRanges(start = as.numeric(vcf.records[1, "POS"]) - 250, 
                                                end = as.numeric(vcf.records[1, "POS"]) + 250))
        )

        # Observe for variant table click
        observeEvent(input$variants_rows_selected, {
            chr <- vcf.records[input$variants_rows_selected, "CHROM"]
            start <- as.numeric(vcf.records[input$variants_rows_selected, "POS"]) - isolate(input$wing_size)
            end <- as.numeric(vcf.records[input$variants_rows_selected, "POS"]) + isolate(input$wing_size)
            robjects$locus <- paste0(chr, ":", start, "-", end)
            robjects$locus_gr <- GRanges(seqnames = chr, ranges = IRanges(start = start, end = end))

            # Get the bam data
            reads <- readGAlignments(bam, param = ScanBamParam(which = robjects$locus_gr, what = "seq"))
            showGenomicRegion(session, id = "igv", robjects$locus)
            loadBamTrackFromLocalData(session, id = "igv", data = reads, trackName = "Alignments")

            # And variants.
            vc <- readVcf(vcf, genome, param = robjects$locus_gr)
            loadVcfTrack(session, id = "igv", trackName = "Variants", vcfData = vc)
        })

        # Create the IGV viewer
        output$igv <- renderIgvShiny({
            igvShiny(opts)
        })

        output$variants <- renderDT({
            datatable(vcf.records, selection = "single", rownames = FALSE)
        })
    }

    shinyApp(ui, server)
}