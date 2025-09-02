retrieve_genesets = function(gobject,GeneSet){
  gene_level_reference = c("low","med","high")
  names(gene_level_reference) = sapply(gene_level_reference,FUN=function(x){toupper(substring(x,1,1))})
  
  names(gene_level_reference)
  if(grepl(x=GeneSet,pattern="hvg")){
    gs = gobject@gene_metadata[get(GeneSet)=="yes",gene_ID]
  }else if(grepl(x=GeneSet,pattern="svg")){
    gs = gobject@gene_metadata[get(GeneSet)==1,gene_ID]
       if(!is.null(gobject@custom_expr)){gs = intersect(gs,rownames(gobject@custom_expr))}       
  }else if(grepl(x=GeneSet,pattern="non_union")){
    gs = gobject@gene_metadata[get(GeneSet)=="yes",gene_ID]
  }else if(grepl(x=GeneSet,pattern="union")){
    levels = strsplit(GeneSet,split="_")[[1]][2]
    HVG_method = strsplit(GeneSet,split="_")[[1]][3]
    SVG_method = strsplit(GeneSet,split="_")[[1]][4]
    HVG_level = paste(HVG_method,"hvg",gene_level_reference[substring(levels,1,1)],sep="_")
    SVG_level = paste(SVG_method,"svg",gene_level_reference[substring(levels,2,2)],sep="_")
    
    hvgs = gobject@gene_metadata[get(HVG_level)=="yes",gene_ID]
    svgs = gobject@gene_metadata[get(SVG_level)==1,gene_ID]
       if(!is.null(gobject@custom_expr)){svgs = intersect(svgs,rownames(gobject@custom_expr))} 
    gs = union(hvgs,svgs)
  }else if(grepl(x=GeneSet,pattern="all")){       
    if(!is.null(gobject@custom_expr)){
            gs=rownames(gobject@custom_expr)
          }else{
            gs=rownames(gobject@norm_expr)
          }
  }
  
  return(gs)
}