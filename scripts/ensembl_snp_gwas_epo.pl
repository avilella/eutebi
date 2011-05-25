#!/usr/bin/perl
use strict;
use Getopt::Long;
use Bio::EnsEMBL::Registry;

my ($inputfile,$debug,$flanks,$query,$target);
$flanks = 2500;
$query  = 'homo_sapiens';
$target = 'rattus_norvegicus';
GetOptions(
	   'i|input|inputfile:s' => \$inputfile,
	   'f|flanks:s' => \$flanks,
           'd|debug:s' => \$debug,
           'q|query:s'  => \$query,
           't|target:s' => \$target,
          );

Bio::EnsEMBL::Registry->load_registry_from_db
  (-host=>"ensembldb.ensembl.org",
   -user=>"anonymous",
   -db_version=>'62');
Bio::EnsEMBL::Registry->no_version_check(1) unless ($debug);

# variation
my $var_adaptor = Bio::EnsEMBL::Registry->get_adaptor($query,"variation","variation");

# compara
my $comparaDBA = Bio::EnsEMBL::Registry->get_DBAdaptor('compara', 'compara');
my $mlss_adaptor = $comparaDBA->get_MethodLinkSpeciesSetAdaptor;
my $gdb_adaptor = $comparaDBA->get_GenomeDBAdaptor;
my $query_genome_db  = $gdb_adaptor->fetch_by_name_assembly($query);
my $slice_adaptor = Bio::EnsEMBL::Registry->get_adaptor($query_genome_db->name, 'core','Slice');
my $target_genome_db = $gdb_adaptor->fetch_by_name_assembly($target);
throw("no slice") if (!$slice_adaptor);
throw("no mlss") if (!$mlss_adaptor);

my $genomic_align_tree_adaptor =
Bio::EnsEMBL::Registry->get_adaptor('Multi', 'compara','GenomicAlignTree');
throw("Cannot get GenomicAlignTree adaptors") if (!$genomic_align_tree_adaptor);

my $method_link_species_set = $mlss_adaptor->fetch_by_method_link_type_species_set_name("EPO", "mammals");
throw("Cannot get GenomicAlignTree") if (!$method_link_species_set);

my $query_short  = $query_genome_db->short_name;
my $target_short = $target_genome_db->short_name;
my @targets;
push @targets, $target; #TODO multiple species

# header
print "rsid\tphenotype\tassociated_gene\tsource\texternal_ref\tvfid\tquery\tquery_seq_region_name\tquery_strand\tquery_start\tquery_end\ttarget\ttarget_seq_region_name\ttarget_strand\ttarget_start\ttarget_end\n";

open FILE,"$inputfile" or die $!;
while (<FILE>) {
  chomp $_;
  next if ($_ =~ /^variation/);
  my ($rsid,$phenotype,$associated_gene,$source,$external_ref) = split("\t",$_);
  my $this_rsid;
  $this_rsid->{rsid}{$rsid}{phenotype} = $phenotype;
  $this_rsid->{rsid}{$rsid}{associated_gene} = $associated_gene;
  $this_rsid->{rsid}{$rsid}{source} = $source;
  $this_rsid->{rsid}{$rsid}{external_ref} = $external_ref;
  my $var = $var_adaptor->fetch_by_name($rsid);
  if ($var->is_failed()) {    my $desc = $var->failed_description();    warn("# $desc. Skip\n");    next;  }
  if ($var->is_somatic()) {    warn("# Somatic mutation. Skip\n");    next;  }

  foreach my $vf (@{$var->get_all_VariationFeatures()}) {
    my $slice = $vf->feature_Slice();
    my $query_slice = $slice->expand($flanks,$flanks);
    $this_rsid->{rsid}{$rsid}{vfid}{$vf->dbID}{seq_region_name} = $query_slice->seq_region_name;
    $this_rsid->{rsid}{$rsid}{vfid}{$vf->dbID}{strand} = $query_slice->strand;
    $this_rsid->{rsid}{$rsid}{vfid}{$vf->dbID}{start} = $query_slice->start;
    $this_rsid->{rsid}{$rsid}{vfid}{$vf->dbID}{end} = $query_slice->end;

    my $genomic_align_trees = $genomic_align_tree_adaptor->fetch_all_by_MethodLinkSpeciesSet_Slice($method_link_species_set,$query_slice);
    next if (!$genomic_align_trees);

    foreach my $this_genomic_align_tree (@$genomic_align_trees) {
      my $this_restricted_genomic_align_tree = $this_genomic_align_tree->restrict_between_reference_positions($query_slice->start,$query_slice->end);
      my $genomic_align_nodes =
        $this_restricted_genomic_align_tree->get_all_sorted_genomic_align_nodes();
      my $other_name = '';
      foreach my $genomic_align_node (@$genomic_align_nodes) {
        my $name = $genomic_align_node->name;
        if ($genomic_align_node->genomic_align_group->genome_db->name eq $query_genome_db->name) {
          $other_name = $name;
        } elsif ($genomic_align_node->genomic_align_group->genome_db->name eq $target_genome_db->name) {
          my $genomic_align_group = $genomic_align_node->genomic_align_group();
          my $genomic_aligns = $genomic_align_group->get_all_GenomicAligns();
          foreach my $genomic_align (@$genomic_aligns) {
            my $target_slice = $genomic_align->get_Slice;
            next unless defined($target_slice);
            my $length = sprintf("%08d",$target_slice->length);
            my $target_seq_region_name = $target_slice->seq_region_name;
            my $target_strand = $target_slice->strand;
            my $target_start = $target_slice->start;
            my $target_end = $target_slice->end;
            $this_rsid->{rsid}{$rsid}{vfid}{$vf->dbID}{target}{$target}{seq_region_name} = $target_seq_region_name;
            $this_rsid->{rsid}{$rsid}{vfid}{$vf->dbID}{target}{$target}{strand} = $target_strand;
            $this_rsid->{rsid}{$rsid}{vfid}{$vf->dbID}{target}{$target}{start} = $target_start;
            $this_rsid->{rsid}{$rsid}{vfid}{$vf->dbID}{target}{$target}{end} = $target_end;
          }
        }
      }
    }
  }
  foreach my $rsid (keys %{$this_rsid->{rsid}}) {
    foreach my $vfid (keys %{$this_rsid->{rsid}{$rsid}{vfid}}) {
      foreach my $target (@targets) {
        my $phenotype = $this_rsid->{rsid}{$rsid}{phenotype};
        my $associated_gene = $this_rsid->{rsid}{$rsid}{associated_gene};
        my $source = $this_rsid->{rsid}{$rsid}{source};
        my $external_ref = $this_rsid->{rsid}{$rsid}{external_ref};
        my $query_seq_region_name = $this_rsid->{rsid}{$rsid}{vfid}{$vfid}{seq_region_name};
        my $query_strand = $this_rsid->{rsid}{$rsid}{vfid}{$vfid}{strand};
        my $query_start           = $this_rsid->{rsid}{$rsid}{vfid}{$vfid}{start};
        my $query_end             = $this_rsid->{rsid}{$rsid}{vfid}{$vfid}{end};
        my $target_seq_region_name = $this_rsid->{rsid}{$rsid}{vfid}{$vfid}{target}{$target}{seq_region_name} || 'na';
        my $target_strand = $this_rsid->{rsid}{$rsid}{vfid}{$vfid}{target}{$target}{strand} || 'na';
        my $target_start = $this_rsid->{rsid}{$rsid}{vfid}{$vfid}{target}{$target}{start} || 'na';
        my $target_end = $this_rsid->{rsid}{$rsid}{vfid}{$vfid}{target}{$target}{end} || 'na';
        $DB::single=$debug;1;
        print "$rsid\t$phenotype\t$associated_gene\t$source\t$external_ref\t$vfid\t$query\t$query_seq_region_name\t$query_strand\t$query_start\t$query_end\t$target\t$target_seq_region_name\t$target_strand\t$target_start\t$target_end\n";
      }
    }
  }
}

close FILE;

1;

# Mapping of human GWAS positive regions to orthologous regions in the
# rat genome

# Identification of orthologous regions of interest within the rat
# genome

# Analysis of the multispecies alignments (MSA) provided by the
# ENSEMBL project. The MSA available through ENSEMBL use the
# Enrado-Pecan-Ortheous (EPO) pipeline, which incorporates a
# consistency-based multiple alignment algorithm in combination with
# the ability to properly align duplicated regions of the genome and
# accurately model the indel structure with the genomes. The EPO
# pipeline is currently the most comprehensive mapping of orthology
# among the sequenced mammalian genomes. Using these results we are
# able to pinpoint the regions in the rat genome that are orthologous
# to positive GWAS regions in the human genome and identified QTL
# regions in the mouse genome whether or not these regions contain
# protein coding genes."

# Get all of the GWAS SNPs from variation and we map them to rat.

# Take both the 5kb and 100 bp ("region" and "pinpointed area")
# centred on the human GWAS hit and get their orthologous positions in
# rat.

# STEPS:

# Install the Ensembl APIs:
# http://www.ensembl.org/info/docs/api/api_installation.html

# briefly: install the bioperl-1.2.3 code, ensembl, ensembl-variation and ensembl-compara API code (requires mysql client and DBD::mysql):

# mkdir -p ~/eutebi/ensembl_main
# cd ~/eutebi/ensembl_main
# cvs -d :pserver:cvsuser@cvs.sanger.ac.uk:/cvsroot/ensembl checkout -r branch-ensembl-62 ensembl
# cvs -d :pserver:cvsuser@cvs.sanger.ac.uk:/cvsroot/ensembl checkout -r branch-ensembl-62 ensembl-compara
# cvs -d :pserver:cvsuser@cvs.sanger.ac.uk:/cvsroot/ensembl checkout -r branch-ensembl-62 ensembl-variation
# cd ~/eutebi
# wget http://bioperl.org/DIST/old_releases/bioperl-1.2.3.tar.gz -O /tmp/bioperl-1.2.3.tar.gz
# tar xzf /tmp/bioperl-1.2.3.tar.gz

# Put it in your path:
# export PERL5LIB=$HOME/eutebi/ensembl_main/ensembl/modules:$HOME/eutebi/ensembl_main/ensembl-hive/modules:$HOME/eutebi/ensembl_main/ensembl-variation/modules:$HOME/eutebi/ensembl_main/ensembl-compara/modules

# Check that the ensembl, ensembl-variation and ensembl-compara modules are visible:
# perldoc Bio::EnsEMBL::Variation::Study
# perldoc Bio::EnsEMBL::Gene
# perldoc Bio::EnsEMBL::Compara::GenomeDB

# Find the list of gwas snps (here for human v62):

# mysql -hensembldb.ensembl.org -P5306 -uanonymous homo_sapiens_variation_62_37g -e "SELECT   v.name AS variation,   p.description AS phenotype,   va.associated_gene,   s.name AS source,   st.external_reference FROM   variation_annotation va JOIN   variation v ON (v.variation_id = va.variation_id) JOIN   phenotype p ON (p.phenotype_id = va.phenotype_id) JOIN   study st ON (st.study_id = va.study_id) JOIN   source s ON (s.source_id = st.source_id) limit 10"
# +------------+-----------+-----------------+--------------------+--------------------+
# | variation  | phenotype | associated_gene | source             | external_reference |
# +------------+-----------+-----------------+--------------------+--------------------+
# | rs11206801 | AB1-42    | Intergenic      | NHGRI_GWAS_catalog | pubmed/20932310    | 
# | rs2075650  | AB1-42    | TOMM40          | NHGRI_GWAS_catalog | pubmed/20932310    | 
# | rs239713   | AB1-42    | Intergenic      | NHGRI_GWAS_catalog | pubmed/20932310    | 
# | rs1727638  | AB1-42    | Intergenic      | NHGRI_GWAS_catalog | pubmed/20932310    | 
# | rs12534221 | AB1-42    | Intergenic      | NHGRI_GWAS_catalog | pubmed/20932310    | 
# | rs10784496 | AB1-42    | Intergenic      | NHGRI_GWAS_catalog | pubmed/20932310    | 
# | rs2899472  | AB1-42    | CYP19A1         | NHGRI_GWAS_catalog | pubmed/20932310    | 
# | rs7631605  | P-tau181p | Intergenic      | NHGRI_GWAS_catalog | pubmed/20932310    | 
# | rs7558386  | P-tau181p | Intergenic      | NHGRI_GWAS_catalog | pubmed/20932310    | 
# | rs12643654 | P-tau181p | UNC5C           | NHGRI_GWAS_catalog | pubmed/20932310    | 
# +------------+-----------+-----------------+--------------------+--------------------+

# Then save it to a file, "> input.tsv", and query the file with this script:
# perl $HOME/eutebi/ensembl_snp_gwas_epo.pl -i input.tsv -target rattus_norvegicus -flanks 2500

# It will produce an extended tabular output like this:

# rsid        phenotype  associated_gene  source              external_ref     vfid      query         query_seq_region_name  query_strand  query_start  query_end  target             target_seq_region_name  target_strand  target_start  target_end
# rs11206801  AB1-42     Intergenic       NHGRI_GWAS_catalog  pubmed/20932310  8906441   homo_sapiens  1                      1             56848186     56853186   rattus_norvegicus  5                       -1             126299643     126303992
# rs2075650   AB1-42     TOMM40           NHGRI_GWAS_catalog  pubmed/20932310  31031874  homo_sapiens  19                     1             45393119     45398119   rattus_norvegicus  1                       -1             79014495      79021635
# rs239713    AB1-42     Intergenic       NHGRI_GWAS_catalog  pubmed/20932310  1564067   homo_sapiens  21                     1             28693976     28698976   rattus_norvegicus  na                      na             na            na
# rs1727638   AB1-42     Intergenic       NHGRI_GWAS_catalog  pubmed/20932310  14142279  homo_sapiens  6                      1             72137072     72142072   rattus_norvegicus  9                       -1             22106476      22116695
# rs12534221  AB1-42     Intergenic       NHGRI_GWAS_catalog  pubmed/20932310  2523508   homo_sapiens  7                      1             131285490    131290490  rattus_norvegicus  4                       1              58692324      58695414
# rs10784496  AB1-42     Intergenic       NHGRI_GWAS_catalog  pubmed/20932310  21541316  homo_sapiens  12                     1             66158471     66163471   rattus_norvegicus  7                       -1             59810525      59814266
# rs2899472   AB1-42     CYP19A1          NHGRI_GWAS_catalog  pubmed/20932310  23872662  homo_sapiens  15                     1             51513555     51518555   rattus_norvegicus  8                       1              57661176      57664187
# rs7631605   P-tau181p  Intergenic       NHGRI_GWAS_catalog  pubmed/20932310  17769149  homo_sapiens  3                      1             37232089     37237089   rattus_norvegicus  na                      na             na            na
# rs7558386   P-tau181p  Intergenic       NHGRI_GWAS_catalog  pubmed/20932310  4614730   homo_sapiens  2                      1             227559639    227564639  rattus_norvegicus  9                       1              81552030      81559132
# rs12643654  P-tau181p  UNC5C            NHGRI_GWAS_catalog  pubmed/20932310  27696656  homo_sapiens  4                      1             96157317     96162317   rattus_norvegicus  2                       -1             239657177     239662228
