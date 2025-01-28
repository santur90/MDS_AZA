# input.rnk
sed -n '2,$p' input.txt | perl -alne '@x=split(/\./,$F[0]);print join("\t",$x[0],$F[4])' | grep -v "NA" > input.rnk
# GSEA
gsea-cli.sh GSEAPreranked -nperm 1000 -rnk input.rnk -rpt_label x -plot_top_x y -out ./ -collapse No_Collapse -gmx z.gmt -set_min 15 -set_max 500 -zip_report false
