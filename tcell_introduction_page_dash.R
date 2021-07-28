HTML("
     
<div>
  <h1>The Malaria Host Atlas</h1>
  <p>
    Malaria remains a global health problem with over 400,000 deaths annually. Plasmodium parasites, the causative agents of malaria, replicate asexually in red blood cells (RBCs) of their vertebrate host, while a subset differentiates into sexual stages (gametocytes) for mosquito transmission. Parasite replication and gametocyte maturation in the erythropoietic niches of the bone marrow and spleen contribute to pathogenesis and drive transmission, but the mechanisms underlying this organ enrichment remain unknown. We performed a comprehensive single cell analysis of rodent P. berghei in spleen, bone marrow and blood to define parasite phenotypes specific to those niches. After quality control, we obtained over 19,000 host single cell transcriptomes from spleen and bone marrow and over 33,000 P. berghei transcriptomes from all three organs, allowing us to investigate host response to infection as well as parasite-specific adaptation to host organ and host cell. Single cell RNA-seq analysis of host and parasite cells revealed an interferon-driven host response to infection as well as transcriptional adaptations of Plasmodium to host organ and RBC maturation status. Our data provides a thorough characterisation of host-parasite interactions in erythropoietic niches and defines host cell maturation state as the key driver of parasite adaptation.
    
  <p>
  Here, you can explore the parasite transcriptomes derived from blood, spleen, and bone marrow. The data set includes the complete asexual replication cycle as well as parasite cells differentiating into male and female sexual transmission forms (gametocytes). The results are displayed in 2 sections

<div>
  <b>Clustering</B> - Allows cluster visualization and exploration of top cluster markers. Each cluster was annotated as a parasite stage based on comparison with the malaria cell atlas - see: 

  
    <a href=' https://www.pnas.org/content/early/2020/08/25/1921930117' target='_blank'>BioRxiv: Host cell maturation modulates parasite invasion and sexual differentiation in <i>Plasmodium</i></a></h4>
    <p>Hentzschel, Franziska; Gibbins, Matthew; Attipa, Charalampos; Beraldi, Dario; Moxon, Chris; Otto, Thomas & Marti, Matthias
    </p>
  </div>
  
  
  <div>
   <b>Differential expression (DE)</b> - Comparison of gene expression of different parasite stages across between different organs (spleen, blood, bone marrow), using the settings specified in our publication (i.e the first 24 principal components and 0.51 as the clustering resolution). Note that the differential gene expression analysis was performed using a different tool in our publication, thus, results might vary slightly.  
  

  </div>
  
  
  <div style='border:thin solid black ; padding:30px ; margin:30px'>
    <figure>  
      <h4></h4>
      <img src='CA_figure1.png' alt='Results' style='width:100%;'>
      <figcaption>
        scRNAseq analysis of P. berghei and host cells from spleen, bone marrow and blood. 
        A) Cartoon depicting the experimental strategy used to enrich infected RBCs from spleen, bone marrow and blood. 
        Three days post infection, we enriched for P. berghei-infected splenic, bone marrow and blood cells by flow sorting, 
        labelled surface CD71 and CD44 expression with barcoded antibodies (Cellular Indexing of Transcriptomes and Epitopes 
        by Sequencing, CITE-seq18) and analysed host and parasite transcriptomes by droplet-based scRNA-seq. 
        B) Annotated UMAP of P. berghei cells. C) UMAP of P. berghei cells separated according to host organ. 
        D) Heatmap depicting average expression levels of genes in late trophozoites that are differentially expressed 
        between organs.
        
      </figcaption>
    </figure>
  </div>
  
</div>")
