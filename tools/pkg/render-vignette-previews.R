#!/usr/bin/env Rscript
# Render installed guides and website-only articles, with local help and an index.
# --index-only refreshes navigation from existing HTML without executing examples.
root <- normalizePath('.')
args <- commandArgs(trailingOnly=TRUE)
index.only <- '--index-only' %in% args
args <- args[args != '--index-only']
stopifnot(length(args) <= 1L)
out <- if(length(args)) normalizePath(args[1],mustWork=FALSE) else file.path(root,'output','vignette-previews')
for(sub in c('doc','html','articles')) dir.create(file.path(out,sub),recursive=TRUE,showWarnings=FALSE)
config <- yaml::read_yaml('_pkgdown.yml')
article.names <- unique(unlist(lapply(config$articles, `[[`, 'contents')))
website.articles <- article.names[startsWith(article.names, 'articles/')]
if (!index.only) {
options(rgl.useNULL=TRUE, grip.local.preview=TRUE)
if (file.exists('../ivue/DESCRIPTION')) pkgload::load_all('../ivue',quiet=TRUE)
if (pkgbuild::needs_compile(root)) pkgbuild::compile_dll(root,debug=FALSE)
pkgload::load_all(root,quiet=TRUE,export_all=FALSE,helpers=FALSE,recompile=FALSE)
for(input in list.files('vignettes',pattern='[.]Rmd$',full.names=TRUE)) {
  message('Rendering ',basename(input))
  rmarkdown::render(input,output_dir=file.path(out,'doc'),
    envir=new.env(parent=globalenv()),quiet=TRUE)
}
for (name in website.articles) {
  status <- system2(file.path(R.home('bin'),'Rscript'),
    c('tools/pkg/render-representation-vignette.R', shQuote(basename(name)),
      shQuote(file.path(out,'articles'))))
  if (status != 0L) stop('Website article render failed: ',name)
}
for(input in list.files('man',pattern='[.]Rd$',full.names=TRUE)) {
  tools::Rd2HTML(input,out=file.path(out,'html',sub('[.]Rd$','.html',basename(input))),package='grip')
}
file.copy(file.path(R.home('doc'),'html','R.css'),file.path(out,'html','R.css'),overwrite=TRUE)
}
redirect <- function(from,to) {
  writeLines(sprintf('<!doctype html><html lang="en"><meta charset="utf-8"><meta http-equiv="refresh" content="0; url=%s"><title>Moved</title><a href="%s">Continue to the current vignette</a></html>',to,to),from)
}
for(file in setdiff(list.files(file.path(out,'doc'),pattern='[.]html$'),
                    paste0(basename(website.articles),'.html')))
  redirect(file.path(out,'articles',file),paste0('../doc/',file))
for(entry in config$redirects) {
  destination <- sub('^articles/','../doc/',entry[2])
  redirect(file.path(out,entry[1]),destination)
}
escape <- function(s) as.character(htmltools::htmlEscape(s))
sections <- vapply(config$articles,function(section) {
  items <- vapply(section$contents,function(name) {
    meta <- rmarkdown::yaml_front_matter(file.path('vignettes',paste0(name,'.Rmd')))
    website.only <- startsWith(name,'articles/')
    target <- if (website.only) paste0(name,'.html') else paste0('doc/',name,'.html')
    if (!file.exists(file.path(out,target))) stop('Missing preview: ',target,
      '. Render it before refreshing the index.')
    sprintf('<li><a href="%s">%s</a>%s</li>',target,escape(meta$title),
      if (website.only) ' <small>(website article)</small>' else '')
  },'')
  paste0('<section><h2>',escape(section$title),'</h2><ul>',paste(items,collapse=''),' </ul></section>')
},'')
dir.create(file.path(out,'figures'),showWarnings=FALSE)
ids <- c('dwt_1005','1138_bus','carpet','porous')
labels <- c('Structural mesh','Power network','Recursive carpet','Porous volume')
for(id in ids) file.copy(file.path('inst','extdata','documentation-display',paste0(id,'.png')),
  file.path(out,'figures',paste0(id,'.png')),overwrite=TRUE)
cards <- paste(vapply(seq_along(ids),function(i) sprintf(
  '<a href="doc/synthetic-graph-families.html#graph-%s" style="display:block;text-decoration:none"><img src="figures/%s.png" alt="%s metric-MDS layout" style="width:100%%;height:auto"><strong>%s</strong></a>',
  ids[i],ids[i],labels[i],labels[i]),''),collapse='')
writeLines(paste0('<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>grip vignettes</title><style>body{max-width:1000px;margin:48px auto;padding:0 24px;font:18px/1.6 system-ui;color:#243347;background:#fcfcfa}a{color:#1769aa}h1{font-size:2.8rem;line-height:1.1}h2{font-size:1.3rem}main{display:grid;grid-template-columns:repeat(auto-fit,minmax(300px,1fr));gap:20px}section{padding:18px 24px;background:#edf3f8;border-radius:10px}li{margin:10px 0}.version{color:#637487;font-size:.85em}</style><h1>Explore graphs with grip</h1><p>Start with a drawing, find your way around the functions, then explore methods and applications. Browse installed vignettes and interactive website articles below. All links open local previews; website-only articles are marked.</p><div style="display:grid;grid-template-columns:repeat(auto-fit,minmax(200px,1fr));gap:15px;margin:25px 0">',cards,'</div><main>',paste(sections,collapse=''),'</main><p class="version">Base commit ',system2('git',c('rev-parse','--short','HEAD'),stdout=TRUE),' (see source checksums in build-manifest.tsv) • grip ',read.dcf('DESCRIPTION')[1,'Version'],' • Index refreshed ',format(Sys.time(),tz='UTC',usetz=TRUE),'</p></html>'),file.path(out,'index.html'))
# Keep numerical/render provenance intact during an index-only refresh.
if (!index.only) {
# Evidence used by this render: checksums distinguish local edits from the base commit.
inputs <- c('DESCRIPTION','_pkgdown.yml','tools/pkg/render-vignette-previews.R',
  list.files('R',pattern='[.]R$',full.names=TRUE),
  list.files('vignettes',pattern='[.]Rmd$',full.names=TRUE,recursive=TRUE),
  list.files('inst/scripts',pattern='[.]R(md)?$',full.names=TRUE),
  list.files('inst/extdata',pattern='[.](rds|png)$',recursive=TRUE,full.names=TRUE))
write.table(data.frame(path=inputs,md5=unname(tools::md5sum(inputs))),
  file.path(out,'build-manifest.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
}
message(if (index.only) 'Refreshed navigation in ' else 'Rendered guides and local help in ',out)
