#!/usr/bin/env Rscript
args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=1) stop("Usage: build-report.R EXPERIMENT_OUTPUT")
output<-normalizePath(args[1])
script<-sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)[1])
source_dir<-dirname(normalizePath(script))
summary<-read.csv(file.path(output,"summary.csv"))
graphs<-read.csv(file.path(output,"graph-results.csv"))
fits<-read.csv(file.path(output,"optimizer-status.csv"))
audit<-read.csv(file.path(output,"extended-budget-audit.csv"))
validation<-readRDS(file.path(output,"postrun-validation.rds"))
stopifnot(length(unique(graphs$replicate[graphs$n==250]))==100,
          length(unique(graphs$replicate[graphs$n==500]))==100)
unique_graphs<-graphs[!duplicated(graphs[,c("n","replicate","k")]),]
texnum<-function(x) {
  s<-formatC(x,format="e",digits=2); parts<-strsplit(s,"e")[[1]]
  sprintf("$%s\\times10^{%d}$",parts[1],as.integer(parts[2]))
}
values<-list(samplesTotal=nrow(unique(graphs[,c("n","replicate")])),graphsTotal=nrow(unique_graphs),
  obstructionsTotal=sum(unique_graphs$obstruction),edgeFits=sum(fits$optimizer=="edge-KK"),
  pathFits=sum(fits$optimizer=="L-BFGS-B"),pathConverged=sum(fits$convergence==0,na.rm=TRUE),
  pathLimited=sum(fits$convergence==1,na.rm=TRUE),pathOther=sum(!is.na(fits$convergence)&!fits$convergence%in%c(0,1)),
  pathCheckError=texnum(validation$max_independent_R_path_score_difference),
  nearTieGraphs=validation$near_tie_graphs,nearTieRelative=texnum(validation$max_relative_distance_difference),
  strictScoreChange=texnum(validation$max_strict_distance_score_difference))
writeLines(vapply(names(values),function(name) sprintf("\\newcommand{\\%s}{%s}",name,values[[name]]),character(1)),
  file.path(output,"report-numbers.tex"))
writeLines(sprintf("\\renewcommand{\\reportbuilddatetime}{%s}",format(Sys.time(),"%Y-%m-%d %H:%M:%S %Z",tz="America/New_York")),
  file.path(output,"report-build-info.tex"))
rule_tex<-c(k_conn="$k_{\\rm conn}$",`k_conn+2`="$k_{\\rm conn}+2$",`2*k_conn`="$2k_{\\rm conn}$")
lines<-c("\\begin{table}[H]\\centering\\small",
  "\\caption{Best-found path error and accompanying edge error, in percent, over 100 samples per row.}",
  "\\begin{tabular}{rlrrr}\\toprule",
  "$n$ & Neighborhood & Path median [IQR] & Edge median [IQR] & Obstructions\\\\\\midrule")
for(n in c(250,500)) for(rule in names(rule_tex)) {
  d<-summary[summary$n==n&summary$rule==rule&summary$method=="Best path candidate",]
  p<-d[d$metric=="path_rel",];e<-d[d$metric=="edge_rel",]
  g<-graphs[graphs$n==n&graphs$rule==rule,]
  lines<-c(lines,sprintf("%d & %s & %.3f [%.3f, %.3f] & %.2f [%.2f, %.2f] & %d/%d\\\\",
    n,rule_tex[[rule]],100*p$median,100*p$q25,100*p$q75,100*e$median,100*e$q25,100*e$q75,sum(g$obstruction),nrow(g)))
}
writeLines(c(lines,"\\bottomrule\\end{tabular}\\end{table}"),file.path(output,"report-summary-table.tex"))
lines<-c("\\begin{table}[H]\\centering\\small",
  "\\caption{Selected-case additional-budget check. Errors are percentages. Status 0: numerical stopping rule met; 1: iteration limit.}",
  "\\begin{tabular}{rlrrrrr}\\toprule",
  "$n$ & Neighborhood & Sample & Path before & Path after & Edge after & Status\\\\\\midrule")
for(i in seq_len(nrow(audit))) {
  a<-audit[i,]
  lines<-c(lines,sprintf("%d & %s & %d & %.4g & %.4g & %.3g & %d\\\\",a$n,rule_tex[[a$rule]],a$replicate,
    100*a$before_path,100*a$after_path,100*a$after_edge,a$convergence))
}
writeLines(c(lines,"\\bottomrule\\end{tabular}\\end{table}"),file.path(output,"report-audit-table.tex"))
old<-setwd(output);on.exit(setwd(old))
for(i in 1:2) {
  status<-system2("pdflatex",c("-interaction=nonstopmode","-halt-on-error",shQuote(file.path(source_dir,"report.tex"))),
    stdout="report-build.log",stderr="report-build.err")
  if(status!=0) stop("PDF build failed; inspect report-build.log")
}
cat(file.path(output,"report.pdf"),"\n")
