#!/usr/bin/perl
use strict;
use Getopt::Long;
use Bio::EnsEMBL::Registry;

my $self = bless {};
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
my $pop_adaptor = Bio::EnsEMBL::Registry->get_adaptor($query,"variation","population");

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
        my $target_seq_region_name = $this_rsid->{rsid}{$rsid}{vfid}{$vfid}{target}{$target}{seq_region_name};
        my $target_strand = $this_rsid->{rsid}{$rsid}{vfid}{$vfid}{target}{$target}{strand};
        my $target_start = $this_rsid->{rsid}{$rsid}{vfid}{$vfid}{target}{$target}{start};
        my $target_end = $this_rsid->{rsid}{$rsid}{vfid}{$vfid}{target}{$target}{end};
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

