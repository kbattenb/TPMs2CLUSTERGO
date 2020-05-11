#!/usr/bin/perl

#########################################
#                                       #
#     Written by Kai Battenberg         #
#     Plant Symbiosis Research Team     #
#                                       #
#########################################

use strict;
use warnings;
use Getopt::Long;

#####SCRIPT DESCRIPTION#####
#Script "TPMs2CLUSTERGO.pl" given a set of TPM files, runs PCA, differential gene expression, hierarchical clustering, and GO enrichment analysis.
###########



######Options#####
print "\n";
print "checking options.\n";

#script directory
my $script_dir = "/Users/kai/Desktop/TPMs2CLUSTERGO/TPMs2CLUSTERGO";
my $pca_r = $script_dir."/PCA.R";
my $dgea_r = $script_dir."/DifferentialGeneExpressionAnalysis.R";
my $ca_r = $script_dir."/ClusteringAnalysis.R";
my $go_r = $script_dir."/GOEnrichmentAnalysis.R";
my $gomodule_r = $script_dir."/GOModulePlot.R";
my $sum_r = $script_dir."/Summurize.R";

#getting tool versions
my $r_version = `Rscript --version 2>&1 | rev | cut -f2 -d' ' | rev | perl -pe chomp`;

#setting the default values.
my $experiment_name = ""; #experiment name
my $treatment_list = ""; #file that includes all treatments and all respective samples.
my $exp_design = ""; #file that includes the experimental design matrix of all samples.
my $tpm_dir = ""; #directory that includes all the files with TPM data.
my $gene_annotation = ""; #gene annotation file with GO terms, Pfam, KEGG ID, etc.
my $out_dir = "";

my $fdr = 0.05; #Benjamini-Hochberg false discovery rate

my $module_size = 30; #minimum number of genes per cluster. Decrease to have finer clusters.
my $deepsplit = 4; #sensitivity to splitting clusters. Increase to have finer clusers. (integers 0-4)
my $cutheight = "AUTO"; #maximum joining height allowed. Decrease to have finer clusters. (0 to 1 in increments of 0.05 or AUTO)
my $cluster_sign = "signed"; #type of correlation distance. ('signed' or 'unsigned')

my $go_min_count = 3;

my $log = "TPMs2CLUSTERGO_log.txt";
my $threads = 2; #number of threads given
##########



#####Checking Options#####
#making the options into external arguments.
GetOptions (
	'experiment_name=s' => \$experiment_name,
	'treatment_list=s' => \$treatment_list,
	'exp_design=s' => \$exp_design,
	'tpm_dir=s' => \$tpm_dir,
	'gene_annotation=s' => \$gene_annotation,
	'out_dir=s' => \$out_dir,
	'fdr=s' => \$fdr,
	'module_size=s' => \$module_size,
	'deepsplit=s' => \$deepsplit,
	'cutheight=s' => \$cutheight,
	'cluster_sign=s' => \$cluster_sign,
	'go_min_count=s' => \$go_min_count,
	'threads=s' => \$threads
	);

#checking for required options.
if (!$experiment_name) {
	die "USAGE: option --experiment_name <EXPERIMENT NAME> is required.\n";
}
elsif (!$treatment_list) {
	die "USAGE: option --treatment_list <TREATMENT LIST FILE> is required.\n";
}
elsif (!$exp_design) {
	die "USAGE: option --exp_design <EXPERIMENTAL DESIGN FILE> is required.\n";
}
elsif (!$tpm_dir) {
	die "USAGE: option --tpm_dir <TPM DIR> is required.\n";
}
elsif (!$gene_annotation) {
	die "USAGE: option --gene_annotation <GENE ANNOTATION> is required.\n";
}
elsif (!$out_dir) {
	die "USAGE: option --out_dir is required.\n";
}

#checking annotation file
if (-e $treatment_list) {
	print " Treatment list file checked.\n";
}
else {
	die "Error: Selected treatment list file $treatment_list does not exist.\n";
}

#checking option tpm_dir
if (-e $tpm_dir && -d $tpm_dir) {
	print " TPM directory checked.\n";
}
else {
	die "Error: Selected TPM directory $tpm_dir does not exist.\n";
}

#checking option out_dir
if (-e $out_dir && -d $out_dir) {
	print " Output directory checked.\n";
}
else {
	die "Error: Selected output directory $out_dir does not exist.\n";
}
chdir $out_dir;

#checking annotation file
if (-e $gene_annotation) {
	print " Gene annotation file checked.\n";
}
else {
	die "Error: Selected gene annotation file $gene_annotation does not exist.\n";
}

#checking FDR
if ( $fdr !~ /^-?\d+\.?\d*$/ ) {
	die "Error: option --fdr needs to be numeric\n";
}
elsif ($fdr >= 1 || $fdr <= 0 ) {
	die "Error: option --fdr needs to be within 0 to 1.\n";
}

#checking minimum module size
if ( $module_size !~ /^-?\d+\.?\d*$/ ) {
	die "Error: option --module_size needs to be numeric\n";
}
elsif ( $module_size < 1 || $module_size !~ /^-?\d+$/ ) {
	die "Error: option --module_size needs to be an integer greater than 0\n";
}

#checking deep split
if ( $deepsplit !~ /^-?\d+\.?\d*$/ ) {
	die "Error: option --deepsplit needs to be numeric\n";
}
elsif ( $deepsplit < 0 || $deepsplit > 4 || $deepsplit !~ /^-?\d+$/ ) {
	die "Error: option --deepsplit needs to be an integer between 0 and 4\n";
}

#checking deep split
if ( $deepsplit !~ /^-?\d+\.?\d*$/ ) {
	die "Error: option --deepsplit needs to be numeric\n";
}
elsif ( $deepsplit < 0 || $deepsplit > 4 || $deepsplit !~ /^-?\d+$/ ) {
	die "Error: option --deepsplit needs to be an integer between 0 and 4\n";
}

#checking deep split
if ( $deepsplit !~ /^-?\d+\.?\d*$/ ) {
	die "Error: option --deepsplit needs to be numeric\n";
}
elsif ( $deepsplit < 0 || $deepsplit > 4 || $deepsplit !~ /^-?\d+$/ ) {
	die "Error: option --deepsplit needs to be an integer between 0 and 4\n";
}

#checking deep split
if ( $cutheight !~ /^-?\d+\.?\d*$/ && $cutheight !~ m/^AUTO$/ ) {
	die "Error: option --cutheight needs to be numeric or AUTO\n";
}
elsif ( $cutheight !~ m/^AUTO$/ ) {
	if ( $cutheight < 0 || $cutheight > 1 ) {
		die "Error: option --cutheight needs to be within 0 and 1\n";
	}
	elsif ( $cutheight*100 % 5 > 0 ) {
		die "Error: option --cutheight needs to be divisible by 0.05\n";
	}
}

#checking cluster sign
if ($cluster_sign !~ m/^signed$/ && $cluster_sign !~ m/^unsigned$/ ) {
	die "Error: option --cluster_sign needs to be 'signed' or 'unsigned'\n";
}

#checking minimum gene count for GO enrichment
if ( $go_min_count !~ /^-?\d+\.?\d*$/ ) {
	die "Error: option --go_min_count needs to be numeric\n";
}
elsif ( $go_min_count < 1 || $go_min_count !~ /^-?\d+$/ ) {
	die "Error: option --go_min_count needs to be a positive integer\n";
}

print "check complete\n\n";
##########



#####Setting output file names#####
#output directory
my $out = $out_dir."/".$experiment_name."_TPMs2CLUSTERGO";

#file names
$log = $experiment_name."_".$log;
open (LOG, ">$log") or die "cannot open $log.\n";

my $input_folder = "01_INPUT";
my $input_tpm_file = $input_folder."/".$experiment_name."_TPMs.txt";
my $input_count_file = $input_folder."/".$experiment_name."_ExpReadCounts.txt";
my $input_treatment_file = $input_folder."/".$experiment_name."_Treatment.txt";
my $input_design_file = $input_folder."/".$experiment_name."_Design.txt";
my $input_annotation_file = $input_folder."/".$experiment_name."_GeneAnnotation.txt";

my $analyses_folder = "02_Analyses";
my $pca_1v2_file = $analyses_folder."/".$experiment_name."_PCA_PC1vsPC2.pdf";
my $pca_2v3_file = $analyses_folder."/".$experiment_name."_PCA_PC2vsPC3.pdf";
my $pca_contribution_file = $analyses_folder."/".$experiment_name."_PCA_PercentContribution.pdf";
my $pca_contributor_file = $analyses_folder."/".$experiment_name."_PCA_Top10ContributingGenes.txt";
my $dgea_threshold_file = $analyses_folder."/".$experiment_name."_DGEA_BestQuantileThreshold.pdf";
my $dgea_deg_file = $analyses_folder."/".$experiment_name."_DGEA_DifferentiallyExpressedGenes.txt";
#my $dgea_maplot_file = $analyses_folder."/".$experiment_name."_DGEA_MAplot.pdf";
my $ca_cor_gene_file = $analyses_folder."/".$experiment_name."_CLUSTER_CorrGene.txt";
my $ca_cor_sample_file = $analyses_folder."/".$experiment_name."_CLUSTER_CorrSample.txt";
my $ca_hierar_gene_file = $analyses_folder."/".$experiment_name."_CLUSTER_HierarchyGene.tre";
my $ca_hierar_sample_file = $analyses_folder."/".$experiment_name."_CLUSTER_HierarchySample.tre";
my $ca_modules_file = $analyses_folder."/".$experiment_name."_CLUSTER_Modules.txt";
my $ca_heatmapscore_file = $analyses_folder."/".$experiment_name."_CLUSTER_HeatmapScore.txt";
my $ca_heatmap_all_file = $analyses_folder."/".$experiment_name."_CLUSTER_HeatmapAll.pdf";
my $ca_heatmap_best_file = $analyses_folder."/".$experiment_name."_CLUSTER_HeatmapBest.pdf";
my $go_enrichment_file = $analyses_folder."/".$experiment_name."_GO_EnrichedTerms.txt";
my $go_module_file = $analyses_folder."/".$experiment_name."_GO_Module.pdf";

my $sum_file = $experiment_name."_summary.txt";
##########



#####STEP1: Make input file folder#####
print "Step-1: Making input file folder\n";
system "rm -rf $input_folder";
system "mkdir $input_folder";

#copying raw input file
print " generating TPM and readcount files with treatment info.\n";

#sample to treatment
print "  reading in treatments and sample names for each treatment.\n";
my %treatments;
my %tpms;
my %readcounts;
open (TREATMENT, "<", $treatment_list) or die "cannot open $treatment_list.\n";
while (my $line = <TREATMENT>) {
	my $data = $line;
	chomp $data;
	if ($data !~ m/\:/ || $data =~ m/ /) {
		die "Error: $treatment_list needs to have the following format for each treatment: <TREATMENT>:<SAMPLE1>,<SAMPLE2>,<SAMPLE3>\n";
	}
	
	my @data = split (/:/, $data);
	my $treatment = shift @data;
	my @samples = split (/\,/, $data[0]);
	
	$treatments{$treatment} = \@samples;
	
	foreach my $sample (@samples) {
		my %tpm;
		$tpms{$sample} = \%tpm;
		my %readcount;
		$readcounts{$sample} = \%readcount;
	}
}
close (TREATMENT);
my $treatments = join (",", sort keys %treatments);
my $samples = join (",", sort keys %tpms);

#tpm/readcount for each sample
print "  reading in TPM data from each sample.\n";
my $tpmfiles = `ls $tpm_dir | grep '_TPM_PerGene'`;
chomp $tpmfiles;
my @tpmfiles = split (/\n/, $tpmfiles);

foreach my $tpmfile (@tpmfiles) {
	my $path = $tpm_dir."/".$tpmfile;
	my $samplename = (split (/_/, $tpmfile))[0];
	if (!exists $tpms{$samplename}) {
		die "Error: TPM file for $samplename ($path) is not included in the treatment list ($treatment_list).\n";
	}
	
	print "   reading in from $path\n";
	my %tpm = %{ $tpms{$samplename} };
	my %readcount = %{ $readcounts{$samplename} };
	open (TPM, "<", $path) or die "cannot open $path.\n";
	my $linecount = 0;
	while (my $line = <TPM>) {
		$linecount++;
		if ($linecount > 1) {
			my $data = $line;
			chomp $data;
			my @data = split (/\t/, $data);
			$tpm{$data[0]} = $data[2];
			$readcount{$data[0]} = $data[1];
		}
	}
	close (TPM);
	$tpms{$samplename} = \%tpm;
	$readcounts{$samplename} = \%readcount;
}
foreach my $samplename (sort keys %tpms) {
	my %tpm = %{ $tpms{$samplename} };
	my $genecount = keys %tpm;
	if ($genecount < 1) {
		die "Error: TPM/readcount file is missing in $tpm_dir for $samplename\n";
	}
}

#tpm for each gene
my %gene2tpms;
my %gene2readcounts;
my $samplecount = 0;
foreach my $sample (keys %tpms) {
	$samplecount++;
	
	my %gene2tpm = %{ $tpms{$sample} };
	my %gene2readcount = %{ $readcounts{$sample} };
	
	foreach my $gene (keys %gene2tpm) {
		my %tpms;
		my %readcounts;
		if ($samplecount > 1 && !exists $gene2tpm{$gene}) {
			die "Error: $gene found in $sample is absent in at least one other sample.\n";
		}
		elsif (exists $gene2tpms{$gene}) {
			%tpms = %{ $gene2tpms{$gene} };
			%readcounts = %{ $gene2readcounts{$gene} };
		}
		$tpms{$sample} = $gene2tpm{$gene};
		$gene2tpms{$gene} = \%tpms;
		$readcounts{$sample} = $gene2readcount{$gene};
		$gene2readcounts{$gene} = \%readcounts;
	}
}

#print out tpm data
print "  printing out TPM input file.\n";

open (TPMS, ">", $input_tpm_file) or die "cannot open $input_tpm_file.\n";
my $line1_tpm = "Treatment";
my $line2_tpm = "SampleID";
my @sampleorder_tpm;
foreach my $treatment (sort keys %treatments) {
	my @samples = @{ $treatments{$treatment} };
	foreach my $sample (sort @samples) {
		push (@sampleorder_tpm, $sample);
		$line1_tpm = $line1_tpm."\t$treatment";
		$line2_tpm = $line2_tpm."\t$sample";
	}
}
print TPMS "$line1_tpm\n";
print TPMS "$line2_tpm\n";

foreach my $gene (sort keys %gene2tpms) {
	my $line3plus = $gene;
	my %tpms = %{ $gene2tpms{$gene} };
	foreach my $sample (@sampleorder_tpm) {
		$line3plus = $line3plus."\t".$tpms{$sample};
	}
	print TPMS "$line3plus\n";
}
close (TPMS);

#print out readcount data
print "  printing out ExpectedReadCount input file.\n";

open (READCOUNT, ">", $input_count_file) or die "cannot open $input_count_file.\n";
my $line1_count = "Treatment";
my $line2_count = "SampleID";
my @sampleorder_count;
foreach my $treatment (sort keys %treatments) {
	my @samples = @{ $treatments{$treatment} };
	foreach my $sample (sort @samples) {
		push (@sampleorder_count, $sample);
		$line1_count = $line1_count."\t$treatment";
		$line2_count = $line2_count."\t$sample";
	}
}
print READCOUNT "$line1_count\n";
print READCOUNT "$line2_count\n";

foreach my $gene (sort keys %gene2readcounts) {
	my $line3plus = $gene;
	my %readcounts = %{ $gene2readcounts{$gene} };
	foreach my $sample (@sampleorder_count) {
		$line3plus = $line3plus."\t".$readcounts{$sample};
	}
	print READCOUNT "$line3plus\n";
}
close (READCOUNT);

#check experimental design file
print "  checking consistency across treatment file and design file.\n";

my %indesign;
my $linecount = 0;
open (DESIGN, "<", $exp_design);
while (my $line = <DESIGN>) {
	$linecount++;
	if ($linecount > 1) {
		my $sampledata = $line;
		chomp $sampledata;
		if ($sampledata !~ m/\t/) {
			die "ERROR: Design file needs to be in tab delimited text file.\n";
		}
		my @sampledata = split (/\t/, $sampledata);
		my $samplename = $sampledata[0];
		shift @sampledata;
		$indesign{$samplename} = join (",", @sampledata);
		if (!exists $tpms{$samplename}) {
			die "ERROR: $samplename is present in design file but absent among TPM files.";
		}
		
	}
}
close (DESIGN);
foreach my $samplename (keys %tpms) {
	if (!exists $indesign{$samplename}) {
		die "ERROR: $samplename is present among TPM files but absent in design file.";
	}
}

foreach my $treatment (keys %treatments) {
	my @samplelist = @{ $treatments{$treatment} };
	my $ref = "";
	foreach my $sample (@samplelist) {
		my $conditions = $indesign{$sample};
		if ($ref eq "") {
			$ref = $conditions;
		}
		elsif ($ref ne $conditions) {
			die "ERROR: within design file, $sample has conditions different from other samples within $treatment\n";
		}
	}
}

#copy input files
system "cp $treatment_list $input_treatment_file";
system "cp $exp_design $input_design_file";
system "cp $gene_annotation $input_annotation_file";

#log
print LOG "Step-1:\n";
print LOG " Experiment name:\t$experiment_name\n";
print LOG " Experimental design:\n";
print LOG "  TreatmentName\tTreatmentDesign\tSampleList\n";
foreach my $treatment (keys %treatments) {
	my @samplelist = @{ $treatments{$treatment} };
	my $samplelist = join (",", @samplelist);
	my $design = $indesign{ $samplelist[0] };
	print LOG "  $treatment\t$design\t$samplelist\n";
}
print "STEP1 complete\n\n";
##########



#####STEP2: Run Analyses#####
print "Step-2: Running variouse analyses\n";
print LOG "Step-2:\n";
system "rm -rf $analyses_folder";
system "mkdir $analyses_folder";

###Running PCA###
print " Begin PCA\n";
my $r_pca_out = `Rscript $pca_r $input_tpm_file 2>&1 | grep \"^\\[1\\]\" | perl -pe chomp | sed 's/\"//g' | sed 's/\\[1\\] /@/g' | sed 's/@//' | perl -pe chomp`;
my @r_pca_out = split (/\@/, $r_pca_out);

system "mv PCA_PC1vPC2.pdf $pca_1v2_file";
system "mv PCA_PC2vPC3.pdf $pca_2v3_file";
system "mv PCA_PercentContribution.pdf $pca_contribution_file";
system "mv PCA_Top10Contributors.txt $pca_contributor_file";

#output
print "  Top 25% of genes with the highest covariance was subjected to PCA.\n";

#log
print LOG " (PCA)\n";
print LOG "  R version:\t$r_version\n";
print LOG "  R packages:\tNo special packages used.\n";
print LOG "  Subject:\tTop 25% of genes with the highest covariance\n";

print " PCA complete\n\n";
######

###Running Differential Gene Expression Analysis###
print " Begin Differential Gene Expression Analysis\n";

my $r_dgea_out = `Rscript $dgea_r $input_count_file $input_design_file $fdr 2>&1 | grep \"^\\[1\\]\" | perl -pe chomp | sed 's/\"//g' | sed 's/\\[1\\] /@/g' | sed 's/@//' | perl -pe chomp`;
my @r_dgea_out = split (/\@/, $r_dgea_out);

system "mv DGEA_BestQuantile_threshold.pdf $dgea_threshold_file";
system "mv DGEA_DEGs.txt $dgea_deg_file";
#system "mv DGEA_MAplot.pdf $dgea_maplot_file";

#extract data
my $filter_threshold = (split (/ /, $r_dgea_out[2]))[-1];
my $deg_count = (split (/ /, $r_dgea_out[3]))[-1];

#output
print "  All genes were subject to Differential gene expression analysis.\n";
print "  Genes with the lowest CPM at $filter_threshold quantile were filtered out for Differential Gene Expression Analysis.\n";
print "  Total of $deg_count genes with adjusted p-value less than $fdr were kept.\n";

#log
print LOG " (Differential Gene Expression Analysis)\n";
print LOG "  R version:\t$r_version\n";
print LOG "  R packages:\t$r_dgea_out[0], $r_dgea_out[1]\n";
print LOG "  Subject:\tAll genes\n";
print LOG "  Lowest CPM filtering threshold:\t$filter_threshold\n";
print LOG "  Benjamini-Hochberg FDR:\t$fdr\n";
print LOG "  Differentially expressed gene count:\t$deg_count\n";

print " Differential Gene Expression Analysis complete\n\n";
######

###Running Clustering Analysis###
print " Begin Clustering Analysis\n";
my $r_cluster_out = `Rscript $ca_r $input_tpm_file $dgea_deg_file $cluster_sign $module_size $deepsplit $cutheight 2>&1 | grep \"^\\[1\\]\" | perl -pe chomp | sed 's/\"//g' | sed 's/\\[1\\] /@/g' | sed 's/@//' | perl -pe chomp`;
my @r_cluster_out = split (/\@/, $r_cluster_out);

system "mv CLUSTER_corrilation_gene.txt $ca_cor_gene_file";
system "mv CLUSTER_corrilation_sample.txt $ca_cor_sample_file";
system "mv CLUSTER_hierarchy_gene.tre $ca_hierar_gene_file";
system "mv CLUSTER_hierarchy_sample.tre $ca_hierar_sample_file";
system "mv CLUSTER_modules.txt $ca_modules_file";
system "mv CLUSTER_Heatmap.txt $ca_heatmapscore_file";
system "mv CLUSTER_Heatmap_ALL.pdf $ca_heatmap_all_file";
system "mv CLUSTER_Heatmap_BEST.pdf $ca_heatmap_best_file";
system "rm Rplots.pdf";

#extract data
my $dist_method = (split (/ /, $r_cluster_out[1]))[2];
my $agglomeration_method = (split (/ /, $r_cluster_out[2]))[4];
my $treecut_method = (split (/ /, $r_cluster_out[3]))[2];
my $best_height = (split (/ /, $r_cluster_out[4]))[3];
my $module_count = (split (/ /, $r_cluster_out[5]))[3];

#output
print "  All differentially expressed genes were subject to clustering analysis.\n";
print "  Genes were clustered according $agglomeration_method based on $cluster_sign $dist_method distances.\n";
if ($cutheight =~ m/AUTO/) {
	print "  The best cut height was set automatically to $best_height, the lowest height setting among those that resulted in the largest number of modules.\n";
}
else {
	print "  The best cut height was set to $best_height.\n";
}
print "  $module_count modules were found\n";

#log
print LOG " (Clustering Analysis)\n";
print LOG "  R version:\t$r_version\n";
print LOG "  R packages:\t$r_cluster_out[0]\n";
print LOG "  Subject:\tAll differentially expressed genes\n";
print LOG "  Distance method:\t$cluster_sign $dist_method\n";
print LOG "  Clustering method:\t$agglomeration_method\n";
print LOG "  dynamicTreeCut options:\n";
print LOG "   Tree cut method:\t$treecut_method\n";
print LOG "   Minimum Module Size:\t$module_size\n";
print LOG "   DeepSplit:\t$deepsplit\n";
if ($cutheight =~ m/AUTO/) {
	print LOG "   Cut height:\t$best_height(Auto detected)\n";
}
else {
	print LOG "   Cut height:\t$best_height\n";
}
print LOG "  Module count:\t$module_count\n";

print " Clustering Analysis complete\n\n";
######

###Running GO Enrichment Analysis###
print " Begin GO Enrichment Analysis\n";

#GO enrichment
my $r_go_out = `Rscript $go_r $ca_modules_file $input_annotation_file $best_height $fdr $go_min_count 2>&1 | grep \"^\\[1\\]\" | perl -pe chomp | sed 's/\"//g' | sed 's/\\[1\\] /@/g' | sed 's/@//' | perl -pe chomp`;
my @r_go_out = split (/\@/, $r_go_out);
my $topgo_algorithm = (split (/ /, $r_go_out[1]))[3];
my $topgo_statistic = (split (/ /, $r_go_out[2]))[3];
system "mv GO_EnrichedPerModule.txt $go_enrichment_file";

#GO module plot
my $r_gomodule_out = `Rscript $gomodule_r $input_tpm_file $ca_modules_file $go_enrichment_file $best_height 2>&1`;
system "mv GO_module_summary.pdf $go_module_file";

#output
print "  All differentially expressed genes clustered to a module at the best height was subject to GO enrichment analysis.\n";
print "  GO terms enriched beyond adjusted p-value less than $fdr with $go_min_count or more representative genes were identified.\n";

#log
print LOG " (GO Enrichment Analysis)\n";
print LOG "  R version:\t$r_version\n";
print LOG "  R packages:\t$r_go_out[0]\n";
print LOG "  Subject:\tAll clustered genes\n";
print LOG "  topGO algorithm:\t$topgo_algorithm\n";
print LOG "  topGO statistic:\t$topgo_statistic\n";
print LOG "  Benjamini-Hochberg FDR:\t$fdr\n";
print LOG "  Minimum gene count:$go_min_count\n";

print " GO Enrichment Analysis complete\n\n";
######

###Summarizing###
print " Begin summarzigin\n";
my $r_sum_out = `Rscript $sum_r $input_tpm_file $input_design_file $dgea_deg_file $ca_modules_file $input_annotation_file $best_height 2>&1`;
system "mv summary.txt $sum_file";

print " Summarizing complete\n\n";
###
##########



#####Final STEP: Sorting files#####
print "Final step: Organizing output.\n";
system "rm -rf $out";
system "mkdir $out";

system "mv $input_folder $out";
system "mv $analyses_folder $out";
system "mv $sum_file $out";

print "Process complete\n";
print LOG "Process complete\n";
close (LOG);
system "mv $log $out";
##########
__END__
