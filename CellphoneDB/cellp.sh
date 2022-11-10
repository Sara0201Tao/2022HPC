source activate cellp


# cellp2, 表达阈值0.15
cellphonedb method statistical_analysis NaiveGood_meta_cellp2.csv NaiveGood_Ncount.csv --counts-data gene_name --iterations 1000 --threshold 0.15 --pvalue 0.05 --threads 25 --project-name NaiveGood_cellp2_015
cellphonedb method statistical_analysis NaiveBad_meta_cellp2.csv NaiveBad_Ncount.csv --counts-data gene_name --iterations 1000 --threshold 0.15 --pvalue 0.05 --threads 25 --project-name NaiveBad_cellp2_015
cellphonedb method statistical_analysis TreatBad_meta_cellp2.csv TreatBad_Ncount.csv --counts-data gene_name --iterations 1000 --threshold 0.15 --pvalue 0.05 --threads 25 --project-name TreatBad_cellp2_015

cellphonedb plot heatmap_plot NaiveGood_meta_cellp2.csv --output-path ./out/NaiveGood_cellp2_015 --pvalues-path ./out/NaiveGood_cellp2_015/pvalues.txt
cellphonedb plot heatmap_plot NaiveBad_meta_cellp2.csv --output-path ./out/NaiveBad_cellp2_015 --pvalues-path ./out/NaiveBad_cellp2_015/pvalues.txt
cellphonedb plot heatmap_plot TreatBad_meta_cellp2.csv --output-path ./out/TreatBad_cellp2_015 --pvalues-path ./out/TreatBad_cellp2_015/pvalues.txt
