This package contains the source code for isoform quantification with RNA-seq data, along with a few C++ library dependencies. 

Users can download the required boost library from: http://archive.gersteinlab.org/boost/

To compile the source code, go to the solve directory and type make from command line. Then in .bashrc file, add "export LD_LIBRARY_PATH= PATH_TO_PACKAGE/gsl/lib/:PATH_TO_PACKAGE/cppunit/lib/" and "export PATH=PATH_TO_PACKAGE/solve/bin:$PATH"

To run the pipeline, type solve from command line.And the parameters required are shown: 
solve
        log_level(0,1,2,...) proj_name out_prefix
        isoform_format isoforms_path g2i_format gene2isoform_path gene_begin_idx gene_end_idx
  (read_format read_type expected_read_length reads_path total_read_bases)+

log_level determines how much information to output during running; proj_name is the name of the project; out_prefix is the prefix for the output file; isoform_format is the annotation file for isoforms, choices are UCSC_GENE_TXT,GENELETS_GFF3,LH_GENE_TXT,UCSC_GFF,WORMBASE_GFF2; isoforms_path is the path to the isoform annotation file; g2i_format is the format of gene to isoform mapping file, choices are UCSC_GENE2ISOFORM,WORMBASE_GENE2ISOFORMS; gene2isoform_path is the path to the gene to isoform mapping file; gene_begin_idx is the number of the first of the gene to be quantified and gene_end_idx is the number of the last gene to be analyzed; read_format is the format of alinged reads, choices are UCSC_GFF,MRF_SINGLE,UCSC_BED,WORMBASE_GFF3; read_type is the type of reads, choices are SHORT_READ,MEDIUM_READ; expected_read_length is the average read length; reads_path is the path to the file of alignment; total_read_bases is the total number of bases in the alignment file.    

An example for running:
solve 3 042310_sample1_C1 /042310_sample1/bowtie/isoform/C1_gff/ LH_GENE_TXT /annotation/UCSC_knownGene_mm9_nh.interval UCSC_GENE2ISOFORM /annotation/UCSC_knownGene_mm9_nh.isoformMap 27000 27500 MRF_SINGLE SHORT_READ 60 /042310_sample1/bowtie/C1.mrf 2197588920 >/042310_sample1/bowtie/isoform/042310_sample1_C1-27000-27500.isoquants
