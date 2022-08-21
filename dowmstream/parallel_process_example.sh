script=RNA-seq_different_expression_and_data_visualization.R

parallel --link  Rscript $script -s hsa -o org.Hs.eg.db \
-g group.example -f 0.5849 -a {1} -b {2} \
::: control control control ip6 ip6 ins \
::: ip6 ins ip6_ins ins ip6_ins ip6_ins


