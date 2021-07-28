

options(shiny.reactlog=TRUE) 
#setwd("D:/GCRF_UoG/Vicky_JCR_Shiny/Tcells_amp_adopted")
source("tcell_libs_raw_dash.R", local = TRUE)
##Creating the ui object for the user interface

ui <- tagList(
    dashboardPage(title = "Malaria Host Atlas",
                  dashboardHeader(title = "Malaria Host Atlas",
                                  tags$li(class = "dropdown",
                                          tags$style(".main-header {max-height: 90px}"),
                                          tags$style(".main-header .logo {height: 90px}"),
                                          tags$a(href='https://www.gla.ac.uk/', target="_blank",
                                                 tags$img(src='uog_logo.png', height = 60, width = 180))
                                  ),
                                  titleWidth = 380),
                  dashboardSidebar(
                      tags$style(".left-side, .main-sidebar {padding-top: 90px}"),
                      width = 250,
                      sidebarMenu(
                          id = 'sidebar',
                          menuItem("About", tabName = "About", icon = icon("door-open")),
                          menuItem("Cluster exploration", tabName = "cluster_res", icon = icon("puzzle-piece"),
                                                           menuSubItem("UMAP", tabName = "all_cluster_res"),
                                                           menuSubItem("Cluster markers", tabName = "grps_cluster_res")
                                                          
                          ),
                          menuItem("Differential expression (DE)", tabName = "ra_de", icon = icon("balance-scale"),
                                                           menuItem("Gene view", tabName = "vis", icon = icon("braille"),
                                                                                            menuSubItem("Single gene view", tabName = "ge_vis_gene"),
                                                                                            menuSubItem("Multiple gene view", tabName = "m_ge_vis_gene") 
                                                                                          
                                                           ),
                                                           menuItem("Cell population view", tabName = "ge_vis_cell", icon = icon("braille")) 
                                                           )
                          
                          
                          
                      )
                  ),
                  dashboardBody(
                      tags$head(
                          tags$link(rel = "stylesheet", type = "text/css", href = "tcell_amp_custom_dash.css"),
                          includeHTML("tcell_js.htm")
                      ),
                      # Boxes need to be put in a row (or column)
                      tabItems(
                          # First tab content
                          # source(file.path("tcell_introduction_page_dash.R"), local = TRUE)$value,
                          tabItem(tabName = "About",
                                  source(file.path("tcell_introduction_page_dash.R"), local = TRUE)$value),
                          
                          tabItem(tabName = "all_cluster_res",
                                  
                                  wellPanel(fluidRow(
                                      
                                      column(
                                          12,
                                          tabBox(
                                              title = "UMAP",
                                              # The id lets us use input$tabset1 on the server to find the current tab
                                              #id = "vis_de", 
                                              #height = "250px",
                                              width = 12,
                                              tabPanel("All cells", 
                                                       wellPanel(
                                                           sliderInput(inputId = "clusters_res", label = strong("Louvain algorithm resolution"), value = res1, min = res1, max = res2, step = diff_res, round = F),
                                                           withSpinner(plotOutput("labelled_umap", width = "600px", height = "400px")),
                                                           wellPanel(
                                                               h4("Download specifications"),
                                                               flowLayout(numericInput("lumap_height", "Plot height (cm):", value = 14),
                                                                          numericInput("lumap_width", "Plot width (cm):", value = 14),
                                                                          radioButtons("lumap_format", "File format:", choices = list("PNG" = ".png", "PDF" = ".pdf", "JPEG" = ".jpeg"), selected = ".png", inline = T)),
                                                               downloadButton('dwnl_lumap','Download Plot')
                                                           ),
                                                           # )
                                                           
                                                           tags$hr(),
                                                           uiOutput("cluster_annot"),
                                                           wellPanel(style = "background:#385A4F",
                                                                     tags$hr(),
                                                                     tags$p(style = "font-family:Arial;color:white",
                                                                            paste("The resolution is fixed to 0.51.")
                                                                            
                                                                            
                                                                     )
                                                           )
                                                           
                                                       )
                                              ),
                                              tabPanel("Sample comparison", withSpinner(plotOutput("all_groups", height = "400px")),
                                                       
                                                       wellPanel(
                                                           h4("Download specifications"),
                                                           flowLayout(numericInput("grps_height", "Plot height (cm):", value = 7),
                                                                      numericInput("grps_width", "Plot width (cm):", value = 30),
                                                                      radioButtons("grps_format", "File format:", choices = list("PNG" = ".png", "PDF" = ".pdf", "JPEG" = ".jpeg"), selected = ".png", inline = T)),
                                                           downloadButton('dwnl_grps','Download Plot')
                                                       )
                                              )
                                              
                                          )
                                          
                                      ),
                                      
                                      
                                  )
                                  )
                                  
                                  
                          ),
                          tabItem(tabName = "grps_cluster_res",
                                  
                                  wellPanel(fluidRow(
                                      
                                      column(
                                          12, 
                                          tabBox(
                                              title = "Cluster markers",
                                              # The id lets us use input$tabset1 on the server to find the current tab
                                              id = "vis_de", 
                                              #height = "250px",
                                              width = 12,
                                              tabPanel("Marker table", 
                                                       wellPanel(
                                                           # box(
                                                           # width = NULL,
                                                           # solidHeader = TRUE,
                                                           uiOutput("dyn_clusters"),
                                                           # uiOutput("topclgenes"),
                                                           withSpinner(DTOutput("top_conserved_genes"))
                                                       )
                                                       
                                                       # )
                                              ),
                                              tabPanel("Marker feature plots", 
                                                       uiOutput("top_markers_umap"),
                                                       withSpinner(plotOutput("conserved_markers_umap")),
                                                       wellPanel(
                                                           h4("Download specifications"),
                                                           flowLayout(numericInput("markers_height", "Plot height (cm):", value = 20),
                                                                      numericInput("markers_width", "Plot width (cm):", value = 20),
                                                                      radioButtons("markers_format", "File format:", choices = list("PNG" = ".png", "PDF" = ".pdf", "JPEG" = ".jpeg"), selected = ".png", inline = T)),
                                                           downloadButton('dwnl_markers','Download Plot')
                                                       )
                                              )
                                          )
                                          
                                          
                                      )
                                  ),
                                  uiOutput("box_2_2")
                                  )
                                  
                                  
                          ),
                          tabItem(tabName = "ge_vis_gene",
                                  #fluidRow(
                                  wellPanel(
                                      fluidRow(
                                          column(
                                              12,
                                              h4("Single gene DE visualization"),
                                              selectizeInput(inputId = "de_genes", label = strong("Choose gene:"),choices = NULL, multiple = F)
                                          )
                                      ),
                                      fluidRow(
                                          tabsetPanel(
                                              tabPanel("Violin plot", withSpinner(plotOutput("de_stim_vs_ctrl_vp")),
                                                       wellPanel(
                                                           h4("Download specifications"),
                                                           flowLayout(numericInput("violp_height", "Plot height (cm):", value = 12),
                                                                      numericInput("violp_width", "Plot width (cm):", value = 30),
                                                                      radioButtons("violp_format", "File format:", choices = list("PNG" = ".png", "PDF" = ".pdf", "JPEG" = ".jpeg"), selected = ".png", inline = T)),
                                                           downloadButton('dwnl_violp','Download Plot')
                                                       )
                                              ),
                                              tabPanel("UMAP feature plot", withSpinner(plotOutput("de_stim_vs_ctrl_um")),
                                                       wellPanel(
                                                           h4("Download specifications"),
                                                           flowLayout(numericInput("featurep_height", "Plot height (cm):", value = 10),
                                                                      numericInput("featurep_width", "Plot width (cm):", value = 30),
                                                                      radioButtons("featurep_format", "File format:", choices = list("PNG" = ".png", "PDF" = ".pdf", "JPEG" = ".jpeg"), selected = ".png", inline = T)),
                                                           downloadButton('dwnl_featurep','Download Plot')
                                                       )
                                              )
                                          ),
                                          
                                          
                                      ),
                                      tags$hr(),
                                      uiOutput("box_1_1")
                                      
                                  )
                          ),
                          tabItem(tabName = "m_ge_vis_gene",
                                  fluidRow(
                                      wellPanel(
                                          h4("Dotplot for multiple gene DE visualization"),
                                          selectizeInput(inputId = "select_markers_dotplot", label = strong("Choose gene:"), choices = NULL, multiple = T),
                                          withSpinner(plotOutput("marker_dotplot")),
                                          wellPanel(
                                              h4("Download specifications"),
                                              flowLayout(numericInput("dotp_height", "Plot height (cm):", value = 10),
                                                         numericInput("dotp_width", "Plot width (cm):", value = 30),
                                                         radioButtons("dotp_format", "File format:", choices = list("PNG" = ".png", "PDF" = ".pdf", "JPEG" = ".jpeg"), selected = ".png", inline = T)),
                                              downloadButton('dwnl_dotp','Download Plot')
                                          ),
                                          tags$hr(),
                                          uiOutput("box_1_2")
                                          
                                          # )
                                      )
                                  )
                          ),
                          tabItem(tabName = "ge_vis_cell",
                                  
                                  wellPanel(
                                      
                                      fluidRow(
                                          column(
                                              4,
                                              uiOutput("cluster_ids"),
                                              ##Comparing conditions new addition
                                              selectInput(inputId = "ra_conds", label = strong("Choose conditions to compare:"),choices = conds)
                                          )
                                      ),
                                      fluidRow(
                                          
                                          column(6,
                                                 wellPanel(
                                                     h4("DE scatterplot"),
                                                     fluidRow(
                                                         withSpinner(plotOutput("cell_type_plot", click = clickOpts(id ="plot_click"))),
                                                         dataTableOutput("click_info"),
                                                         wellPanel(
                                                             h4("Download specifications"),
                                                             flowLayout(numericInput("scatter_height", "Plot height (cm):", value = 14),
                                                                        numericInput("scatter_width", "Plot width (cm):", value = 14),
                                                                        radioButtons("scatter_format", "File format:", choices = list("PNG" = ".png", "PDF" = ".pdf", "JPEG" = ".jpeg"), selected = ".png", inline = T)),
                                                             downloadButton('dwnl_scatter','Download Plot')
                                                         ),
                                                         
                                                         
                                                         
                                                     ),
                                                     tags$hr(),
                                                     uiOutput("box_1_3a")
                                                     
                                                 )
                                          ),
                                          
                                          column(
                                              6,
                                              wellPanel(
                                                  h4("DE table"),
                                                  
                                                  # sliderInput(inputId = "top_genes", label = strong("Number of top DE genes:"), value = 100, min = 1, max = dim(cluster)[1], step = 1),
                                                  # uiOutput("topdegenes"),
                                                  withSpinner(DTOutput("top_de_genes")),
                                                  tags$hr(),
                                                  uiOutput("box_1_3b")
                                                  
                                                  # )
                                              )
                                          )
                                      )
                                      
                                  )
                                  
                          )
                      )
                  )
    ),
    includeHTML("tcell_footer_dash.htm")
)

# alveri = readRDS("alveri.rds")

##server function to compute the outputs
server = function(input, output, session) {
    
    
    updateSelectizeInput(session = session, inputId = 'de_genes', choices = all_genes_common_in_all_groups, selected = choice_gene, server = TRUE)
    updateSelectizeInput(session = session, inputId = 'select_markers_dotplot', choices = all_genes_common_in_all_groups, selected = fav_genes, server = TRUE)
    
    ##Selecting the Seurat object
    umap_clusters = reactive({
        
        tcells_combined_umap_list_res_skinny[[(input$clusters_res + 0.49)]]
    })
    
    ######################
    ### Labelling clusters under subtitle 2.3
    ######################
    ##Generating dynamic fields for labelling UMAP clusters and initializing the fields with placeholders
    
    ##Preparing and plotting UMAP cluster markers for annotating cell types
    umap_cluster_modified_rna = reactive({
        umap_cluster_modified_ul = umap_clusters()
        DefaultAssay(umap_cluster_modified_ul) = "RNA"
        umap_cluster_modified_ul
    })
    
    output$cluster_annot <- renderUI({
        
        if(length(unique(umap_cluster_modified_rna()[["seurat_clusters"]][["seurat_clusters"]])) == length(cluster_names)){
            do.call(flowLayout, 
                    lapply(0:(length(unique(umap_cluster_modified_rna()[["seurat_clusters"]][["seurat_clusters"]]))-1), function(x) {
                        textInput(inputId = paste("labeller", x), label = strong(paste("Input cluster", x, "label")), value = cluster_names[x+1])
                    })
            )
            
        } else {
            do.call(flowLayout,
                    lapply(0:(length(unique(umap_cluster_modified_rna()[["seurat_clusters"]][["seurat_clusters"]]))-1), function(x) {
                        textInput(inputId = paste("labeller", x), label = strong(paste("Input cluster", x, "label")), value = x)
                    })
            )
        }
    })
    
    ##Storing names decided on by the researchers based on optimal clustering to start off the differential expression visualization
    annotation = reactiveValues(annot = cluster_names)
    
    ##Observer to allow updating of cluster names dynamically as typed in
    observe({
        
        req(unlist(lapply(0:(length(unique(umap_cluster_modified_rna()[["seurat_clusters"]][["seurat_clusters"]]))-1), function(x) {
            new_cluster_name = input[[paste("labeller",x)]]
        })))
        annotation$annot = unlist(lapply(0:(length(unique(umap_cluster_modified_rna()[["seurat_clusters"]][["seurat_clusters"]]))-1), function(x) {
            new_cluster_name = input[[paste("labeller",x)]]
        }))
    })
    
    ##Dynamic input for selecting celltypes (clusters) for diffential expression visualization
    output$cluster_ids <- renderUI({
        umap_names = annotation$annot
        #FindMarkers(stim_markers(), ident.1 = paste(input$select_cell_type, "KO", sep = "_"), ident.2 = paste(input$select_cell_type, "WT", sep = "_"), verbose = FALSE)
        
        if(length(unique(umap_cluster_modified_rna()[["seurat_clusters"]][["seurat_clusters"]])) == length(umap_names)){
            umap_names = annotation$annot
            selectInput(inputId = "select_cell_type", label = strong(paste("Select cell population to compare gene expression between",conditions[1],"and", conditions[2],":")),choices = umap_names, multiple = F)
            
        } else {
            selectInput(inputId = "select_cell_type", label = strong("Select cell population to compare gene expression across conditions:"),choices = unique(umap_cluster_modified_rna()[["seurat_clusters"]][["seurat_clusters"]]), multiple = F)
        }
        
    })
    
    ##Renaming clusters
    umap_cluster_modified_ren_reo = reactive({
        umap_names = annotation$annot
        umap_cluster_modified1 = umap_cluster_modified_rna()
        if(length(unique(umap_cluster_modified1$seurat_clusters)) == length(umap_names)){
            names(umap_names) <- levels(umap_cluster_modified1)
            umap_cluster_modified1 <- RenameIdents(umap_cluster_modified1, umap_names)
            
        } else {
            umap_cluster_modified1
        }
        
    })
    
    ##Plotting labelled umap
    labelled_umap_r = reactive({
        
        if(length(unique(umap_cluster_modified_ren_reo()[["seurat_clusters"]][["seurat_clusters"]])) == length(cluster.colours)){
            DimPlot(umap_cluster_modified_ren_reo(), order = T,  pt.size = 1, label = TRUE, label.size = 6)#, cols = cluster.colours)
            
        } else {
            DimPlot(umap_cluster_modified_ren_reo(), order = T, pt.size = 1, label = TRUE, label.size = 6)
        }
        
    })
    output$labelled_umap = renderPlot({
        labelled_umap_r()
    })
    
    
    output$dwnl_lumap <- downloadHandler(
        filename = function(){paste("labelled_umap",input$lumap_format,sep="")},
        content = function(file){
            ggsave(file,plot=labelled_umap_r(), width = input$lumap_width, height = input$lumap_height, units = "cm", dpi = 300)
        },
        contentType = "image"
    )
    
    ######################
    ###END Labelling clusters under subtitle 2.3
    ######################
    
    ######################
    ### Selection of seurat object and plotting UMAP plots under subtitle 2.1
    ######################
    
    ##Plotting UMAP plots for clustering
    umap_p_split = reactive({
        withProgress(message = 'Plotting', detail = 'Please wait...', value = 0.8, {
            if(length(unique(umap_cluster_modified_ren_reo()[["seurat_clusters"]][["seurat_clusters"]])) == length(cluster.colours)){
                DimPlot(umap_cluster_modified_ren_reo(),reduction = "umap", split.by = "sample",  pt.size = 1, label = TRUE, label.size = 6)#, cols = cluster.colours)
                
            } else {
                DimPlot(umap_cluster_modified_ren_reo(), reduction = "umap", split.by = "sample", pt.size = 1, label = TRUE, label.size = 6)
            }
            
        })
    })
    
    output$all_groups = renderPlot({
        
        umap_p_split()    
    })
    
    output$dwnl_grps <- downloadHandler(
        filename = function(){paste("umap_split_by_groups",input$grps_format,sep="")},
        content = function(file){
            ggsave(file,plot=umap_p_split(), width = input$grps_width, height = input$grps_height, units = "cm", dpi = 300)
        },
        contentType = "image"
    )
    ######################
    ### END Selection of seurat object and plotting UMAP plots under subtitle 2.1
    ######################
    
    ######################
    ### Cluster markers plots and tables under subtitle 2.2
    ######################
    ##Dynamic input field for selecting cluster to plot table of markers
    output$dyn_clusters <- renderUI({
        selectInput(inputId = "marker_genes_cluster", label = strong("Choose cluster to display markers for"), choices = unique(umap_clusters()[["seurat_clusters"]][["seurat_clusters"]]), multiple = F)
    })
    
    ##Dynamic input for selecting celltypes (clusters) for diffential expression visualization
    output$dyn_clusters <- renderUI({
        umap_names = annotation$annot
        #FindMarkers(stim_markers(), ident.1 = paste(input$select_cell_type, "KO", sep = "_"), ident.2 = paste(input$select_cell_type, "WT", sep = "_"), verbose = FALSE)
        
        if(length(unique(umap_cluster_modified_rna()[["seurat_clusters"]][["seurat_clusters"]])) == length(umap_names)){
            umap_names = annotation$annot
            selectInput(inputId = "marker_genes_cluster", label = strong("Choose cluster to display markers for"), choices = umap_names, multiple = F)
            
        } else {
            selectInput(inputId = "marker_genes_cluster", label = strong("Choose cluster to display markers for"), choices = unique(umap_cluster_modified_rna()[["seurat_clusters"]][["seurat_clusters"]]), multiple = F)
        }
        
    })
    
    ##Displaying table of cluster markers for annotating cell types
    cluster_markers = reactive({
        req(input$marker_genes_cluster)
        umap_names = annotation$annot
        
        if(length(unique(umap_cluster_modified_rna()[["seurat_clusters"]][["seurat_clusters"]])) == length(umap_names)){
            withProgress(message = 'Tabulating', {
                setProgress(detail = 'Please wait...', value = 0.4)
                marker_tb =     tcells_combined_clusters_tables_res[[(input$clusters_res + 0.49)]][[as.numeric(match(input$marker_genes_cluster,umap_names))]] %>% rownames_to_column(var = 'genes') %>% dplyr::left_join(x = ., y = uniprot_info, by = c("genes" = "Gene")) %>% dplyr::distinct(., genes, .keep_all = T) %>% mutate_at(vars(matches("p_val|pval") ), ~formatC(., format = "e", digits = 2)) %>% dplyr::select(genes,`Protein name`, PlasmoDB, contains(c("avg", "adj"))) %>% mutate_if(is.numeric, ~sprintf("%.3f", .))
                setProgress(detail = 'Please wait...', value = 0.8)
                return(marker_tb)
            })
        } else {
            # ra_macrophage_combined_de_tables_full[[1]][[(as.numeric(input$select_cell_type) + 1)]]
            withProgress(message = 'Tabulating', {
                setProgress(detail = 'Please wait...', value = 0.4)
                marker_tb =     tcells_combined_clusters_tables_res[[(input$clusters_res + 0.49)]][[(as.numeric(input$marker_genes_cluster) + 1)]] %>% rownames_to_column(var = 'genes') %>% dplyr::left_join(x = ., y = uniprot_info, by = c("genes" = "Gene")) %>% dplyr::distinct(., genes, .keep_all = T)  %>% mutate_at(vars(matches("p_val|pval") ), ~formatC(., format = "e", digits = 2)) %>% dplyr::select(genes,`Protein name`, PlasmoDB, contains(c("avg", "adj"))) %>% mutate_if(is.numeric, ~sprintf("%.3f", .))
                setProgress(detail = 'Please wait...', value = 0.8)
                return(marker_tb)
            })            
        }
        
    })
    
    
    output$top_conserved_genes = DT::renderDataTable({
        # req(input$topclgenes_i)
        #numeric_cols =  colnames(data.frame(cluster_markers()))[which_numeric_cols(data.frame(cluster_markers()))]
        
        #   # Javascript-enabled table.
        #   datatable(
        DT::datatable(
            cluster_markers(), 
            selection = "single",
            escape = F,
            rownames = FALSE,
            filter = list(position = "top", plain = TRUE),
            style = "default",
            extensions = c("Buttons","Scroller"),
            options = list(
                deferRender = TRUE,
                scrollY = 350,
                #scroller = TRUE,
                #lengthMenu = FALSE,
                autoWidth = FALSE,
                dom = "Blfrtip",
                buttons = 
                    list(list(
                        extend = "collection",
                        buttons = c("csv", "pdf"),
                        text = "Download"
                    )),  # end of buttons customization
                
                # customize the length menu
                lengthMenu = list( c(10, 20, -1), c(10, 20, "All")), # end of lengthMenu customization
                pageLength = 10
            ), fillContainer = TRUE
        )
    }, server = TRUE)
    
    
    output$top_markers_umap <- renderUI({
        req(cluster_markers())
        selectInput(inputId = "select_markers_umap", label = strong("Select marker to visualize in clusters:"), choices = cluster_markers()[,1], multiple = T, selected = head(cluster_markers()[,1], n=4))
        
    })
    conserved_markers_umap_r = reactive({
        req(input$select_markers_umap)
        req(umap_cluster_modified_ren_reo())
        withProgress(message = 'Plotting',{
            shiny::setProgress( detail = 'Please wait...', value = 0.4)
            fp_conserved = FeaturePlot(umap_cluster_modified_rna(), features = input$select_markers_umap, min.cutoff = "q9")
            shiny::setProgress( detail = 'Please wait...', value = 0.8)
            return(fp_conserved)
        })
    })
    
    output$conserved_markers_umap = renderPlot({
        req(input$select_markers_umap)
        req(umap_cluster_modified_ren_reo())
        conserved_markers_umap_r()
    })
    
    
    output$dwnl_markers <- downloadHandler(
        filename = function(){paste(input$select_markers_umap, "_feature_plot",input$markers_format,sep="")},
        content = function(file){
            ggsave(file,plot=conserved_markers_umap_r(), width = input$markers_width, height = input$markers_height, units = "cm", dpi = 300)
        },
        contentType = "image"
    )
    
    ##information box
    output$box_2_2 =renderUI({
        wellPanel(style = "background:#385A4F",
                  tags$hr(),
                  tags$p(style = "font-family:Arial;color:white",
                         paste("Listing top cluster marker which can subsequently be used in labelling the cluster.Here, markers are genes highly expressed in a cluster as compared to all other clusters in both", conditions[1], "and", conditions[2],".")
                  ),
                  tags$hr()
        )
    })
    
    ######################
    ### END Cluster markers plots and tables under subtitle 2.2
    ######################
    
    
    ######################
    ### Differential expression dynamic using UMAP and Violin plots under subtitle 1.1
    ######################
    stim_markers = reactive({
        
        umap_cluster_modified = umap_cluster_modified_ren_reo()
        umap_cluster_modified$celltype.group <- interaction(Idents(umap_cluster_modified), umap_cluster_modified$group, sep = "_")
        umap_cluster_modified$celltype <- Idents(umap_cluster_modified)
        Idents(umap_cluster_modified) <- "celltype.group"
        umap_cluster_modified
    })
    
    
    #Functions to update differentially expressed genes
    de_stim_vs_ctrl_um_r = eventReactive(input$de_genes,{
        withProgress(message = 'Plotting',{
            shiny::setProgress(detail = 'Please wait...', value = 0.5)
            fp_umap = FeaturePlot(stim_markers(), features = input$de_genes, split.by = "group", max.cutoff = 3,cols = c("grey", "red"))
            shiny::setProgress(detail = 'Please wait...', value = 0.8)
            return(fp_umap)
        })
    })
    
    output$de_stim_vs_ctrl_um = renderPlot({
        req(input$de_genes)
        req(de_stim_vs_ctrl_um_r())
        withProgress(message = 'Plotting', detail = 'Please wait...', value = 0.9, {
            
            de_stim_vs_ctrl_um_r()
        })
    })
    
    output$dwnl_featurep <- downloadHandler(
        filename = function(){paste(input$de_genes, "_feature_plot",input$featurep_format,sep="")},
        content = function(file){
            ggsave(file,plot=de_stim_vs_ctrl_um_r(), width = input$featurep_width, height = input$featurep_height, units = "cm", dpi = 300)
        },
        contentType = "image"
    )
    
    
    de_stim_vs_ctrl_vp_r = reactive({
        # plots <- VlnPlot(stim_markers(), features = input$de_genes, split.by = "group", group.by = "celltype", pt.size = 0, combine = FALSE, multi.group = T, cols = group.cols, assay = "RNA") 
        # for(i in 1:length(plots)) {
        #     plots[[i]] <- plots[[i]] + stat_summary(fun.y= median, geom='point', size = 2, colour = "black", position = position_dodge(0.9)) + scale_fill_manual(values=group.cols) 
        # }
        # CombinePlots(plots)
        VlnPlot(stim_markers(), features = input$de_genes, split.by = "group", group.by = "celltype", pt.size = 0, combine = FALSE, cols = group.cols, assay = "RNA") 
        
    })
    
    output$de_stim_vs_ctrl_vp = renderPlot({
        req(input$de_genes)
        req(stim_markers())
        withProgress(message = 'Plotting', detail = 'Please wait...', value = 0.7, {
            
            de_stim_vs_ctrl_vp_r()
        })
    })
    
    output$dwnl_violp <- downloadHandler(
        filename = function(){paste(input$de_genes,"_violin_plot",input$violp_format,sep="")},
        content = function(file){
            ggsave(file,plot=de_stim_vs_ctrl_vp_r(), width = input$violp_width, height = input$violp_height, units = "cm", dpi = 300)
        },
        contentType = "image"
    )
    
    ##Information box
    output$box_1_1 <- renderUI({
        wellPanel(style = "background:#385A4F",
                  tags$hr(),
                  tags$p(style = "font-family:Arial;color:white",
                         
                         paste("Comparison of", input$de_genes, "expression between", cond," across all clusters using violin plots and umap feature plots.")
                         
                  ),
                  tags$hr()
        )  })
    
    ######################
    ### END Differential expression dynamic using UMAP and Violin plots under subtitle 1.1
    ######################
    
    
    ######################
    ### Differential expression using dotplot under subtitle 1.2
    ######################
    ##Dotplot for DE comparison between KO and WT across cell types
    marker_dotplot_r = eventReactive(input$select_markers_dotplot,{
        req(input$select_markers_dotplot)
        req(umap_cluster_modified_ren_reo())
        withProgress(message = 'Plotting',{
            shiny::setProgress(detail = 'Please wait...', value = 0.3)
            
            umap_cluster_modified_ren_reo = umap_cluster_modified_ren_reo()
            umap_cluster_modified_ren_reo@meta.data$grp_od <- umap_cluster_modified_ren_reo@meta.data$group
            umap_cluster_modified_ren_reo@meta.data <- umap_cluster_modified_ren_reo@meta.data %>% mutate(grp_od = case_when(grp_od == "Healthy" ~ 4,grp_od == "UPA" ~ 3,grp_od == "Naive RA" ~ 2,grp_od == "Resistant RA" ~ 1,grp_od == "Remission RA" ~ 0))
            shiny::setProgress(detail = 'Please wait...', value = 0.5)
            
            umap_cluster_modified_ren_reo@meta.data <- dplyr::arrange(umap_cluster_modified_ren_reo@meta.data, umap_cluster_modified_ren_reo@meta.data$grp_od)
            
            ##or **(not the '-' sign)**
            umap_cluster_modified_ren_reo@meta.data <- dplyr::arrange(umap_cluster_modified_ren_reo@meta.data, -umap_cluster_modified_ren_reo@meta.data$grp_od)
            dp = DotPlot(umap_cluster_modified_ren_reo, features = input$select_markers_dotplot, cols = group.cols, dot.scale = 6, split.by = "group") + RotatedAxis()
            shiny::setProgress(detail = 'Please wait...', value = 0.7)
            return(dp)
            
        })
    })
    
    output$marker_dotplot = renderPlot({
        marker_dotplot_r()
    })
    
    
    output$dwnl_dotp <- downloadHandler(
        filename = function(){paste(input$select_markers_dotplot,"dotplot",input$dotp_format,sep="_")},
        content = function(file){
            ggsave(file,plot=marker_dotplot_r(), width = input$dotp_width, height = input$dotp_height, units = "cm", dpi = 300)
        },
        contentType = "image"
    )
    
    ##Information box
    output$box_1_2 <- renderUI({
        wellPanel(style = "background:#385A4F",
                  tags$hr(),
                  tags$p(style = "font-family:Arial;color:white",
                         
                         paste("Comparison of gene expression between", conditions[1], "and", conditions[2], "cells across clusters using a dotplot. The genes are on the y-axis and the clusters on the x-axis. Red is for", conditions[1], "and Blue for", conditions[2], "cells with the increase in intensity of the respective colour (from grey to blue/red) correlating with the increase in the average level of gene expression across all cells in the cluster. The size of the dot corresponds to the percentage of cells in the cluster expressing the gene.")
                         
                  ),
                  tags$hr()
        )})
    ######################
    ###END Differential expression using dotplot under subtitle 1.2
    ######################
    
    
    ######################
    ### Differential expression using ggplot and tables under subtitle 1.3
    ######################
    
    ##Retrieving table for DE expression from precomputed list
    genes_in_de_order = reactive({
        umap_names = annotation$annot
        
        if(length(unique(umap_cluster_modified_rna()[["seurat_clusters"]][["seurat_clusters"]])) == length(umap_names)){
            tcells_combined_de_tables[[(input$clusters_res + 0.49)]][[as.numeric(match(input$select_cell_type,umap_names))]][[as.numeric(match(input$ra_conds,conds))]]
            
        } else {
            # ra_macrophage_combined_de_tables_full[[1]][[(as.numeric(input$select_cell_type) + 1)]]
            tcells_combined_de_tables[[(input$clusters_res + 0.49)]][[(as.numeric(input$select_cell_type) + 1)]][[as.numeric(match(input$ra_conds,conds))]]
            
        }
        
    })
    
    ##Retrieving table for DE expression from precomputed list
    top_de_g = reactive({
        req(input$select_cell_type, input$ra_conds)
        withProgress(message = 'Tabulating',{
            setProgress(detail = 'Please wait...', value = 0.4)
            t_d_g = genes_in_de_order()  %>% rownames_to_column(var = 'genes') %>% filter(p_val_adj <= 0.05) %>% mutate_at(vars(matches("p_val|pval") ), ~formatC(., format = "e", digits = 2)) %>% mutate_if(is.numeric, ~sprintf("%.3f", .))
            setProgress(detail = 'Please wait...', value = 0.8)
            return(t_d_g)
        })
    })
    
    output$top_de_genes = DT::renderDataTable({
        DT::datatable(
            top_de_g(),
            selection = "single",
            rownames = FALSE,
            filter = list(position = "top", plain = TRUE),
            style = "default",
            extensions = c("Buttons","Scroller"),
            options = list(
                deferRender = TRUE,
                scrollY = 350,
                #scroller = TRUE,
                #lengthMenu = FALSE,
                autoWidth = FALSE,
                dom = "Blfrtip",
                buttons =
                    list(list(
                        extend = "collection",
                        buttons = c("csv", "pdf"),
                        text = "Download"
                    )),  # end of buttons customization
                
                # customize the length menu
                lengthMenu = list( c(10, 20, -1), c(10, 20, "All")), # end of lengthMenu customization
                pageLength = 10
            )
        )
        # DT::formatSignif(columns = numeric_cols, digits = 3)
    }, server = TRUE)
    
    
    ##Allowing for download of DE table
    output$downloadData <- downloadHandler(
        filename = function() {
            paste("differentially_expressed_genes_in",input$select_cell_type,input$ra_conds, ".csv", sep = "_")
            
        },
        content = function(file) {
            write.csv(genes_in_de_order() %>% rownames_to_column(var = 'genes') %>% filter(p_val_adj <= 0.05)  %>% mutate_at(vars(matches("p_val|pval") ), ~formatC(., format = "e", digits = 2)), file) %>% mutate_if(is.numeric, ~sprintf("%.3f", .))
        }
    )
    
    ##Retrieving table for DE scatterplotfrom precomputed list
    cell_type_de = reactive({
        req(input$select_cell_type, input$ra_conds) 
        umap_names = annotation$annot
        
        if(length(unique(umap_cluster_modified_rna()[["seurat_clusters"]][["seurat_clusters"]])) == length(umap_names)){
            tcells_combined_de_ggplots_table[[(input$clusters_res + 0.49)]][[as.numeric(match(input$select_cell_type,umap_names))]]
            
        } else {
            tcells_combined_de_ggplots_table[[(input$clusters_res + 0.49)]][[(as.numeric(input$select_cell_type) + 1)]]
        }
        
    })
    
    
    ##preparing ggplot for average DE expression for genes above
    cell_type_de_plot_no_grb = reactive({
        req(input$select_cell_type, input$ra_conds)
        theme_set(theme_cowplot())
        ggplot(data=cell_type_de(), aes_string(paste("`",unlist(str_split(input$ra_conds[1], " VS "))[1],"`", sep=""), paste("`",unlist(str_split(input$ra_conds[1], " VS "))[2],"`", sep=""))) + geom_point() + ggtitle(input$select_cell_type) + theme_bw()
        
        
    })
    cell_type_de_plot = reactive({
        req(input$select_cell_type, input$ra_conds)
        #theme_set(theme_cowplot())
        grob <- grobTree(textGrob("Click on points to diplay more information about the gene", x=0.1,  y=0.95, hjust=0,
                                  gp=gpar(col="red", fontsize=13, fontface="italic")))
        cell_type_de_plot_no_grb() + annotation_custom(grob)
        
        
    })
    
    ##plotting ggplot for average DE expression for genes above
    output$cell_type_plot = renderPlot({
        #theme_set(theme_cowplot())
        cell_type_de_plot()
    })
    
    output$dwnl_scatter <- downloadHandler(
        filename = function(){paste("scatter_plot_of_average_expression_in_", input$ra_conds,"_among_", input$select_cell_type, input$scatter_format,sep="")},
        content = function(file){
            ggsave(file,plot=cell_type_de_plot_no_grb(), width = input$scatter_width, height = input$scatter_height, units = "cm", dpi = 300)
        },
        contentType = "image"
    )
    
    ##Displaying further details upon clicking points
    displayed_text <- reactive({
        req(input$plot_click)
        nearPoints(cell_type_de(), input$plot_click)
        
    })
    
    ##Displaying table with gene details upon click of point in DE scatterplot
    output$click_info <- renderDataTable({
        req(displayed_text())
        displayed_text()
        # DT::datatable(displayed_text(),
        #               extensions=c('Scroller'),
        #               options = list(dom = 'Bfrtip',
        #                              scroller = TRUE,
        #                              scrollX=TRUE)) 
    }, escape = F)
    
    ##Information box
    output$box_1_3a <- renderUI({
        wellPanel(style = "background:#385A4F",
                  tags$hr(),
                  tags$p(style = "font-family:Arial;color:white",
                         paste("Comparison of average gene expression between", input$ra_conds,"in", input$select_cell_type, "cells using a scatter plot.")
                  ),
                  tags$hr()
        )
    })
    
    output$box_1_3b <- renderUI({
        
        wellPanel(style = "background:#385A4F",
                  tags$hr(),
                  tags$p(style = "font-family:Arial;color:white",
                         paste("Listing differentially expressed genes (adjusted P value <0.05) between" , input$ra_conds, "in", input$select_cell_type,"cells.")
                  ),
                  tags$hr()
        )
    })
    ######################
    ### END Differential expression using ggplot and tables under subtitle 1.3
    ######################
    
    
}


shinyApp(ui = ui, server = server)


