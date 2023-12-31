#' Create the ASOspotter Shiny app
#'
#' @param vcf Character scalar of path to a VCF file.
#' @param bams A named list of character scalars of paths to BAM files. 
#'   Element names will be used as track names. 
#' @param beds An optional named list of character scalars of paths to BED files. 
#'   Element names will be used as track names. 
#' @param genome The genome to use. 
#'   See \code{\link[igvShiny]{currently.supported.stock.genomes()}} for all supported genomes.
#' @return A shiny app to interactively view variants.
#'
#' @importFrom vcfR read.vcfR getFIX
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom GenomicAlignments readGAlignments
#' @importFrom Rsamtools ScanBamParam
#' @importFrom igvShiny igvShiny renderIgvShiny igvShinyOutput showGenomicRegion loadBamTrackFromLocalData
#'   parseAndValidateGenomeSpec loadVcfTrack loadBedTrack
#' @importFrom DT renderDT datatable formatStyle
#' @importFrom VariantAnnotation readVcf
#' @importFrom shiny shinyApp addResourcePath observeEvent reactiveValues isolate
#'
#' @author Jared Andrews
#' @export
#' @examples
#' library(ASOspotter)
#' vcf <- system.file("extdata", "NA12878_HG001.hg38.benchmark.subset.vcf.gz", package = "ASOspotter")
#' bam <- system.file("extdata", "NA12878.sub.sorted.bam", package = "ASOspotter")
#' bed <- system.file("extdata", "hg38_cgi.bed", package = "ASOspotter")
#' beds <- list("CpG Islands" = bed)
#' bams <- list("NA12878" = bam, "NA12878_again" = bam)
#' \dontrun{
#' ASOspotter(vcf = vcf, bams = bams, beds = beds)
#' }
ASOspotter <- function(vcf, bams, beds = NULL, genome = "hg38") {
    # Make local dir for track subsets.
    if (!dir.exists("tracks")) {
        dir.create("tracks")
    }
    addResourcePath("tracks", "tracks")

    # Load the VCF
    vcf.records <- read.vcfR(vcf)
    vcf.records <- cbind(as.data.frame(getFIX(vcf.records)), as.data.frame(vcf.records@gt))

    # Get initial locus
    initialLocus <- paste0(
        vcf.records[1, "CHROM"], ":",
        as.numeric(vcf.records[1, "POS"]) - 250, "-", as.numeric(vcf.records[1, "POS"]) + 250
    )

    # Set igvShiny options
    opts <- parseAndValidateGenomeSpec(genomeName = genome, initialLocus = initialLocus)

    # Create the UI
    ui <- .create_ui()

    server <- function(input, output, session) {
        robjects <- reactiveValues(
            locus = initialLocus,
            locus_gr = GRanges(
                seqnames = vcf.records[1, "CHROM"],
                ranges = IRanges(
                    start = as.numeric(vcf.records[1, "POS"]) - 50,
                    end = as.numeric(vcf.records[1, "POS"]) + 50
                )
            )
        )

        # Observe for variant table click
        observeEvent(input$variants_rows_selected, {
            chr <- vcf.records[input$variants_rows_selected, "CHROM"]
            start <- as.numeric(vcf.records[input$variants_rows_selected, "POS"]) - isolate(input$wing_size)
            end <- as.numeric(vcf.records[input$variants_rows_selected, "POS"]) + isolate(input$wing_size)
            robjects$locus <- paste0(chr, ":", start, "-", end)
            robjects$locus_gr <- GRanges(seqnames = chr, ranges = IRanges(start = start, end = end))

            showGenomicRegion(session, id = "igv", robjects$locus)

            # Add/plot the tracks
            for (b in seq_along(bams)) {
                reads <- readGAlignments(bams[[b]], param = ScanBamParam(which = robjects$locus_gr, what = "seq"))
                loadBamTrackFromLocalData(session, id = "igv", data = reads, trackName = names(bams)[b])
            }

            for (be in seq_along(beds)) {
                bed <- read.table(beds[[be]], sep = "\t", header = FALSE, stringsAsFactors = FALSE)
                colnames(bed)[1:3] <- c("chr", "start", "end")
                loadBedTrack(session, id = "igv", tbl = bed, trackName = names(beds)[be])
            }

            # And variants.
            vc <- readVcf(vcf, genome, param = robjects$locus_gr)
            loadVcfTrack(session, id = "igv", trackName = "Variants", vcfData = vc)
        })

        # Create the IGV viewer
        output$igv <- renderIgvShiny({
            igvShiny(opts)
        })

        output$variants <- renderDT({
            datatable(vcf.records, 
                filter = "top",
                extensions = c("Buttons"), 
                selection = "single", 
                rownames = FALSE,
                options = list(
                    search = list(regex = TRUE),
                    pageLength = 10,
                    dom = "Blfrtip",
                    buttons = c("copy", "csv", "excel", "pdf", "print")
                )
            ) %>% formatStyle(0, target = "row", lineHeight = "50%")
        })
    }

    shinyApp(ui, server)
}
