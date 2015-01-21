#!/usr/bin/perl

=head1
 iTAK -- Plant Transcription factor & Protein Kinase Identifier and Classifier

 yz357@cornell.edu

 update
	[Dec-31-2014][v1.5]: combination rules for plantTFDB and plnTFDB, new classification system for future rule update
	[Jan-19-2014][v1.4]: new category system for plant protein kinase, from RKD and HMM build by Shiu et al. 2012 
		     update the hmmscan to version 3.1, 2x faster than hmm 3.0
	[Jun-22-2013][v1.3]: update some small bugs, Pfam V27, WAK and WAKL family, a switch for transponsase filter
	[Aug-26-2011][v1.2]: report unusual sequences
	[Jun-03-2011][v1.1]: remove unsignificiant domain using GA score
	[Dec-14-2010][v1.0]: first stable version 
=cut

use strict;
use warnings;
use Cwd;
use File::Basename;
use Bio::SeqIO;
use IO::File;
use Bio::SearchIO;
use Getopt::Std;
use FindBin;
use lib "$FindBin::RealBin/bin";
use itak;

my $version = 1.5;
my $debug = 0;

my %options;
getopts('a:b:c:d:e:f:g:i:j:k:l:m:n:o:p:q:r:s:t:u:v:w:x:y:z:h', \%options);
unless (defined $options{'t'} ) { usage($version); }

if	($options{'t'} eq 'identify')	{ itak_identify(\%options, \@ARGV); }
elsif	($options{'t'} eq 'database')	{ itak_database(\%options, \@ARGV); }
elsif	($options{'t'} eq 'compare')	{ itak_compare(\%options, \@ARGV); }
else	{ usage($version); }

#################################################################
# kentnf: subroutine						#
#################################################################
=head2
 itak_compare -- compare classification 
=cut
sub itak_compare
{
	my ($options, $files) = @_;

	my $usage = qq'
USAGE: $0 -t compare -f family_name ath_tf_classification.txt 

* please use ath_pep to perform identification
* the PlnTFDB and PlantTFDB is locate in database folder

';
	print $usage and exit unless $$files[0];
	die "[ERR]fdr not exist\n" unless -s $$files[0];
	my $fdr = $$files[0];
	
	my $dbs_dir = ${FindBin::RealBin}."/database";
	my $pln	    = "$dbs_dir/PlnTFDB_Ath.txt";
	my $plant   = "$dbs_dir/PlantTFDB_Ath.txt";
	my $input   = "$fdr/tf_classification.txt";
	my $align   = "$fdr/tf_alignment.txt";

	# load family information to hash
	my $family;
	$family = $$options{'f'} if defined $$options{'f'};

	my %uid;
	my %plnTFDB     = load_family($pln, 2, 3);
	my %plantTFDB   = load_family($plant, 2, 3);
	my %iTAK        = load_family($input, 1, 2);

	# load alignment informat to hash
	my %iTAK_align;
	my $fha = IO::File->new($align) || die $!;
	while(<$fha>)
	{
		chomp;
		my @a = split(/\t/, $_);
		my ($pfam_id, $pfam_name, $score) = ($a[1], $a[9], $a[11]);
		my $combine_info = "$pfam_id;$pfam_name;$score";

		if (defined $iTAK_align{$a[0]}) {
			$iTAK_align{$a[0]}.="\t".$combine_info;
		} else {
			$iTAK_align{$a[0]} = $combine_info;
		}	
	}
	$fha->close;

	# making comparison
	foreach my $ida (sort keys %plnTFDB)	{ $uid{$ida} = 1; }
	foreach my $ida (sort keys %plantTFDB)	{ $uid{$ida} = 1; }
	foreach my $ida (sort keys %iTAK)	{ $uid{$ida} = 1; }

	my @tr = qw/ARID AUX.IAA Coactivator=p15 DDT GNAT HMG IWS1 Jumonji LUG MBF1 MED6 MED7 PHD Pseudo=ARR-B RB Rcd1-like SET SNF2 SOH1 SWI.SNF-BAF60b SWI.SNF-SWI3 TRAF Alfin-like BSD CSD DBP LIM mTERF OFP PLATZ sigma70-like TAZ TIG ULT VARL zn-clus Orphans Sigma70-like TUB Tify/;
	my %tr;
	foreach my $tr (@tr) {
        	$tr =~ s/\./\//;
	        $tr =~ s/=/ /;
	        $tr{$tr} = 1;
	}
	# foreach my $t (sort keys %tr) { print $t."\n"; } exit;

	my %tbs = qw/CO-like C2C2-CO-like Dof C2C2-Dof GATA C2C2-GATA YABBY C2C2-YABBY LSD C2C2-LSD LBD LOB Nin-like RWP-RK ZF-HD zf-HD MYB_related MYB-related/;
	$tbs{'E2F/DP'} = 'E2F-DP';
	$tbs{'NZZ/SPL'} = 'NOZZLE';

	my %tas = qw/PBF-2-like Whirly/;
	$tas{'BBR/BPC'} = "BBR-BPC";

	my ($ta, $tb, $tc);
	print "ID\tPlnTFDB\tPlantTFDB\tiTAK\n";
	foreach my $id (sort keys %uid)
	{
        	defined $plnTFDB{$id} and $ta = $plnTFDB{$id} or $ta = "NA";
	        defined $plantTFDB{$id} and $tb = $plantTFDB{$id} or $tb = 'NA';
	        defined $iTAK{$id} and $tc = $iTAK{$id} or $tc = 'NA';
	        #$ta = $tas{$ta} if defined $tas{$ta};
	        #$tb = $tbs{$tb} if defined $tbs{$tb};
	        #next if ($ta eq $tb && $tb eq $tc);
	        #next if (defined $tr{$ta} && $ta eq $tc && $tb eq 'NA');

	        #if ($ta eq 'MADS' && ($tb eq 'M-type' || $tb eq 'MIKC')) {
	        #        next if $tb eq $tc;
	        #}

	        #if ($ta eq 'AP2-EREBP' && ($tb eq 'AP2' || $tb eq 'ERF' || $tb eq 'RAV')) {
	        #        next if $tb eq $tc;
	        #}

	        #if ($ta eq 'ABI3VP1' && ($tb eq 'B3' || $tb eq 'ARF')) {
	        #        next if $tb eq $tc;
	        #}

	        #if ($ta eq 'CCAAT' && ($tb eq 'NF-YA' || $tb eq 'NF-YB' || $tb eq 'NF-YC')) {
        	#        next if $tb eq $tc;
	        #}

	        #if ($ta eq 'NA' &&  $tb eq 'NF-X1' && $tc eq 'NF-X1') { next; }

	        #if ($ta eq 'NA' &&  $tb eq 'C2C2-LSD' && $tc eq 'C2C2-LSD') { next; }

	        #if ($ta eq "HB") { next; }

		my $alignment = 'No align';
		$alignment = $iTAK_align{$id} if defined $iTAK_align{$id};

		if ($family) {
			if ( $family eq $ta || $family eq $tb || $family eq $tc) { print "$id\t$ta\t$tb\t$tc\t$alignment\n"; }
		} else {
			print "$id\t$ta\t$tb\t$tc\t$alignment\n";
		}
	}
}

# load tf family
sub load_family 
{
	my ($file, $gene_col, $family_col) = @_;
	my %gene_family;
	my $fh = IO::File->new($file) || die "$file $!";
	while(<$fh>)
	{
		chomp;
		next if $_ =~ m/^#/;
		my @a = split(/\t/, $_);
		my $id = $a[$gene_col-1];
		$gene_family{$id} = $a[$family_col-1];
	}
	$fh->close;
	return %gene_family;
}

=head2
 itak_identify -- identification of TFs or PKs
=cut
sub itak_identify
{
	my ($options, $files) = @_;

	my $usage =qq'
USAGE:  perl $0 [options] input_seq 

	-m  [mode]	quick or normal (default: normal) 
        -a  [Integer]   number of CPUs used for hmmscan. (default = 1)
        -o  [String]    Name of the output directory. ( default = \'input file
                        name\' + \'_output\')

';

	# +++++ check input parameters +++++
	print $usage and exit unless $$files[0];
	foreach my $f (@$files) {
		my $output_dir = $f."_output";
		my $temp_dir = $f."_temp";
		print "[WARN]output folder exist: $output_dir\n" if -e $output_dir;
		print "[WARN]temp folder exist: $temp_dir\n" if -e $temp_dir;
		die "[ERR]input file not exist\n" unless -s $f;
	}
	
	my $cpu = '20';
	$cpu = $$options{'p'} if (defined $$options{'p'} && $$options{'p'} > 0); 

	my $mode = 'normal';
	$mode = 'quick' if (defined $$options{'m'} && $$options{'m'} eq 'quick');

	# +++++ set database and script +++++
	my $bin_dir = ${FindBin::RealBin}."/bin";
	my $dbs_dir = ${FindBin::RealBin}."/database";
	unless (-e $bin_dir) { die "[ERR]bin folder not exist.\n$bin_dir\n"; }
	unless (-e $dbs_dir) { die "[ERR]database folder not exist.\n $dbs_dir\n"; }

	my $tfam_db	= $dbs_dir."/Tfam_domain.hmm";		# database for transcription factors (subset of Pfam-A + customized). for quick mode
	my $pfam_db	= $dbs_dir."/Pfam-A.hmm";		# Pfam-A
	my $sfam_db	= $dbs_dir."/TF_selfbuild.hmm";		# self-build 
	my $plantsp_db = $dbs_dir."/PlantsPHMM3_89.hmm";	# plantsP kinase
	my $shiu_db    = $dbs_dir."/Plant_Pkinase_fam.hmm";	# shiu kinase database
	foreach my $db (($tfam_db, $pfam_db, $sfam_db, $plantsp_db, $shiu_db)) {
		die "[ERR]database file $db\n" unless (-s $db.".h3f" && -s $db.".h3i" && -s $db.".h3m" && -s $db.".h3p");
	}

	my $tf_rule = $dbs_dir."/TF_Rule.txt";              	# Rules for Transcription Factors
	my $correct_ga = $dbs_dir."/GA_table.txt";		# update GA cutoff
	my $pk_desc = $dbs_dir."/PK_class_desc.txt";		# PK family description (for PPC)
	my $hmmscan_bin = $bin_dir."/hmmscan";			# hmmscan 
	my $hmmpress_bin = $bin_dir."/hmmpress";		# hmmpress

	foreach my $f (($tf_rule, $correct_ga, $pk_desc, $hmmscan_bin, $hmmpress_bin)) {
		die "[ERR]file not exist: $f\n" unless -s $f;
	}

	my %tf_rule = load_rule($tf_rule);
	my %ga_cutoff;
	%ga_cutoff = load_ga_cutoff($pfam_db, $correct_ga, $sfam_db) if $mode eq 'normal';
	%ga_cutoff = load_ga_cutoff($tfam_db, $correct_ga, $sfam_db) if $mode eq 'quick';
	my $pkid_des = pk_to_hash($pk_desc);

	# +++++ main +++++ 
	foreach my $f (@$files)
	{
		# create folder for temp files and output files
		my $temp_dir = $f."_temp";
		my $output_dir = $f."_output";
		mkdir($temp_dir) unless -e $temp_dir;
		mkdir($output_dir) unless -e $output_dir;

		my $input_protein_f = $temp_dir."/protein_seq.fa";			# input protein sequence
		my $tmp_pfam_hmmscan = $temp_dir."/protein_seq.pfam.hmmscan.txt";	# temp hmmscan result compared to Pfam-A + customized
		my $report_info = '';

		# put seq to hash
		# key: id, alphabet, seq; value: alphabet, seq
		# proteins to temp file
		my %seq_info = seq_to_hash($f);
		$report_info = "Load ".scalar(keys(%seq_info))." sequences.\n\nNot protein sequence ID:\n";

		my $outp = IO::File->new(">".$input_protein_f) || die $!;
		foreach my $id (sort keys %seq_info) {

			if ($seq_info{$id}{'alphabet'} eq 'protein') {
				print $outp ">".$id."\n".$seq_info{$id}{'seq'}."\n";
			} else {
				$report_info.= "$id\n";
			}
		}
		$outp->close;
		print "[ERR]no input proteins\n" unless -s $input_protein_f;

		# ==== Part A TF identification ====
		# ==== A1. compare input seq with database ====
		my ($hmmscan_hit_1, $hmmscan_detail_1);

		if ($mode eq 'normal') {
			my $hmmscan_command1 = "$hmmscan_bin --acc --notextw --cpu $cpu -o $tmp_pfam_hmmscan.a $pfam_db $input_protein_f";
			my $hmmscan_command2 = "$hmmscan_bin --acc --notextw --cpu $cpu -o $tmp_pfam_hmmscan.b $sfam_db $input_protein_f";
			run_cmd($hmmscan_command1) unless -s $tmp_pfam_hmmscan.".a"; # test code
			run_cmd($hmmscan_command2) unless -s $tmp_pfam_hmmscan.".b"; # test code
			($hmmscan_hit_1, $hmmscan_detail_1) = itak::parse_hmmscan_result($tmp_pfam_hmmscan.".a");
			my ($hmmscan_hit_2, $hmmscan_detail_2) = itak::parse_hmmscan_result($tmp_pfam_hmmscan.".b");
			$hmmscan_hit_1.= $hmmscan_hit_2;
			$hmmscan_detail_1.= $hmmscan_detail_2;
		} else {
			my $hmmscan_command = "$hmmscan_bin --acc --notextw --cpu $cpu -o $tmp_pfam_hmmscan $tfam_db $input_protein_f";
			run_cmd($hmmscan_command) unless -s $tmp_pfam_hmmscan; # test code
			($hmmscan_hit_1, $hmmscan_detail_1) = itak::parse_hmmscan_result($tmp_pfam_hmmscan);
		}

		my %align_pfam_hash = aln_to_hash($hmmscan_detail_1, \%ga_cutoff);

		# ==== A2. TF identification ====
		my %qid_tid = itak_tf_identify($hmmscan_hit_1, $hmmscan_detail_1, \%ga_cutoff, \%tf_rule);

		# ==== A3. save the result ====
		my $output_sequence	  = "$output_dir/tf_sequence.txt";
		my $output_alignment	  = "$output_dir/tf_alignment.txt";
		my $output_classification = "$output_dir/tf_classification.txt";
		itak_tf_write_out(\%qid_tid, \%seq_info, $hmmscan_detail_1, \%tf_rule, $output_sequence, $output_alignment, $output_classification);
		
		# ==== Part B PK identification ====
		# ==== B1. get protein kinase sequence ====
		my %pkinase_id;
		chomp($hmmscan_hit_1); my @hit_line = split(/\n/, $hmmscan_hit_1);
		foreach my $line ( @hit_line ) {
			my @a = split(/\t/, $line);
			$a[1] =~ s/\..*//;
			if ($a[1] eq 'PF00069' ) {
				die "[ERR]no GA for PF00069" unless defined $ga_cutoff{$a[1]};
				$pkinase_id{$a[0]} = 1 if $a[2] >= $ga_cutoff{$a[1]};
			} 

			if ($a[1] eq 'PF07714') {
				die "[ERR]no GA for PF07714" unless defined $ga_cutoff{$a[1]};
				$pkinase_id{$a[0]} = 1 if $a[2] >= $ga_cutoff{$a[1]};
			}
		}

		next if scalar(keys(%pkinase_id)) == 0;
		
		my $tmp_pkinase_seq = $temp_dir."/pkinase_seq.fa"; 
		my $out1 = IO::File->new(">".$tmp_pkinase_seq) || die $!;
		foreach my $id (sort keys %pkinase_id) {
			print $out1 ">".$id."\n".$seq_info{$id}{'seq'}."\n";
		}
		$out1->close;

		# ==== B2. compare input seq with databas ====
                my $tmp_plantsp_hmmscan = "$temp_dir/protein_seq.plantsp.hmmscan.txt";
                my $tmp_shiu_hmmscan    = "$temp_dir/protein_seq.shiu.hmmscan.txt";
		#my $tmp_rkd_hmmscan     = "$temp_dir/protein_seq.rkd.hmmscan.txt";

		my $plantsp_hmmscan_cmd = "$hmmscan_bin --acc --notextw --cpu $cpu -o $tmp_plantsp_hmmscan $plantsp_db $tmp_pkinase_seq";
		my $shiu_hmmscan_cmd    = "$hmmscan_bin --acc --notextw --cpu $cpu -o $tmp_shiu_hmmscan    $shiu_db    $tmp_pkinase_seq";
		#my $rkd_hmmscan_cmd     = "$hmmscan_bin --acc --notextw --cpu $cpu -o $tmp_rkd_hmmscan     $rkd_db    $tmp_pkinase_seq";

		run_cmd($plantsp_hmmscan_cmd) unless -s $tmp_plantsp_hmmscan;
		run_cmd($shiu_hmmscan_cmd) unless -s $tmp_shiu_hmmscan;
		my ($plantsp_hit, $plantsp_detail) = itak::parse_hmmscan_result($tmp_plantsp_hmmscan);
		my ($shiu_hit, $shiu_detail)       = itak::parse_hmmscan_result($tmp_shiu_hmmscan);

		# ==== B3. PK classification ====		
		my ($plantsp_cat, $plantsp_aln) = itak_pk_classify($plantsp_detail);
                my ($shiu_cat, $shiu_aln) = itak_pk_classify($shiu_detail);

                # ==== B4 classification of sub pkinase ====
		my @wnk1 = ("$dbs_dir/Pkinase_sub_WNK1.hmm",   "30" , "PPC:4.1.5", "PPC:4.1.5.1");
		my @mak  = ("$dbs_dir/Pkinase_sub_MAK.hmm", "460.15" , "PPC:4.5.1", "PPC:4.5.1.1");
		my @sub = (\@wnk1, \@mak);

		foreach my $s ( @sub ) {
			# check array info for sub classify
			die "[ERR]sub classify info ".join(",", @$s)."\n" unless scalar(@$s) == 4;
			my ($hmm_profile, $cutoff, $cat, $sub_cat) = @$s;

			# get seq for sub classify
			my $seq_num = 0;
			my $ppc_seq = "$temp_dir/temp_ppc_seq";
			my $ppfh = IO::File->new(">".$ppc_seq) || die $!;
			foreach my $seq_id (sort keys %$plantsp_cat) {
				if ( $$plantsp_cat{$seq_id} eq $cat ) {
               				die "[ERR]seq id: $seq_id\n" unless defined $seq_info{$seq_id}{'seq'};
                                	print $ppfh ">".$seq_id."\n".$seq_info{$seq_id}{'seq'}."\n";
                                	$seq_num++;
				}
                        }
			$ppfh->close;

			# next if there is no seq
			next if $seq_num == 0;

			# hmmscan and parse hmm result
        		my $ppc_hmm_result = $temp_dir."/temp_ppc_sub_hmmscan.txt";
        		my $hmm_cmd = "$hmmscan_bin --acc --notextw --cpu $cpu -o $ppc_hmm_result $hmm_profile $ppc_seq";
			run_cmd($hmm_cmd) unless -s $ppc_hmm_result;
        		my ($ppc_hits, $ppc_detail) = itak::parse_hmmscan_result($ppc_hmm_result);
			my @hit = split(/\n/, $ppc_detail);

			foreach my $h (@hit) {
				my @a = split(/\t/, $h);
				if ( $a[9] >= $cutoff ) {
					$$plantsp_cat{$a[0]} = $sub_cat;
					$$plantsp_aln{$a[0]} = $h."\n";
				}
			}
                }

                # %pkinases_cat = get_wnk1(\%pkinases_cat, "$dbs_dir/wnk1_hmm_domain/WNK1_hmm" , "30" ,"PPC:4.1.5", "PPC:4.1.5.1");""
                # $$pkid_des{"PPC:4.1.5.1"} = "WNK like kinase - with no lysine kinase";
                # %pkinases_cat = get_wnk1(\%pkinases_cat, "$dbs_dir/mak_hmm_domain/MAK_hmm" , "460.15" ,"PPC:4.5.1", "PPC:4.5.1.1");
                # $$pkid_des{"PPC:4.5.1.1"} = "Male grem cell-associated kinase (mak)";
        
		#foreach my $pid (sort keys %pkinases_cat) {
                #        print $ca_fh1 $pid."\t".$pkinases_cat{$pid}."\t".$$pkid_des{$pkinases_cat{$pid}}."\n";
                #}
                #$ca_fh1->close;

		# ==== B5 save result =====
		# output plantsp classification
		my $ppc_cat = $output_dir."/PPC_classification.txt";
		my $ppc_aln = $output_dir."/PPC_alignment.txt";

                my $ca_fh1 = IO::File->new(">".$ppc_cat) || die $!;
                my $al_fh1 = IO::File->new(">".$ppc_aln) || die $!;
		foreach my $pid (sort keys %$plantsp_cat) { print $ca_fh1 $pid."\t".$$plantsp_cat{$pid}."\t".$$pkid_des{$$plantsp_cat{$pid}}."\n"; }

                foreach my $pid (sort keys %$plantsp_cat) {
                        if (defined $align_pfam_hash{$pid} && defined $$plantsp_cat{$pid} ) {
				print $al_fh1 $$plantsp_aln{$pid};
                                print $al_fh1 $align_pfam_hash{$pid};
                        } else {
                                die "Error! Do not have alignments in hmm3 parsed result\n";
                        }
                        # delete $pkinase_id{$pid};
                }

                #foreach my $pid (sort keys %plantsp_cat) {
                #        print $ca_fh1 $pid."\tPPC:1.Other\n";
		#
		#	if (defined $pkinase_aln{$pid}) {
                #                print $al_fh1 $pkinase_aln{$pid};
                #        } else {
                #                die "Error! Do not have alignments in hmm3 parsed result\n";
                #        }
                #}

		$ca_fh1->close;
                $al_fh1->close;
		
                # output Shiu classification
		my $shiu_cat_file = $output_dir."/shiu_classification.txt";
                my $shiu_aln_file = $output_dir."/shiu_alignment.txt";

                my $ca_fh2 = IO::File->new(">".$shiu_cat_file) || die $!;
                my $al_fh2 = IO::File->new(">".$shiu_aln_file) || die $!;
                foreach my $pid (sort keys %$shiu_cat) { print $ca_fh2 $pid."\t".$$shiu_cat{$pid}."\n"; }
		
		foreach my $pid (sort keys %$shiu_cat) { 
                        if (defined $align_pfam_hash{$pid} && defined $$shiu_cat{$pid} ) {
                                print $al_fh2 $$shiu_aln{$pid};
                                print $al_fh2 $align_pfam_hash{$pid};
                        } else {
                                die "Error! Do not have alignments in hmm3 parsed result\n";
                        }
		}

                $ca_fh2->close;
                $al_fh2->close;	

		# output pkinase sequences
		my $pkinase_seq = $output_dir."/pk_sequence.txt";
		my $out_pks = IO::File->new(">".$pkinase_seq) || die $!;
		foreach my $pid (sort keys %pkinase_id) {
			my $cat1 = 'NA'; my $cat2 = 'NA';
			$cat1 = $$plantsp_cat{$pid} if defined $$plantsp_cat{$pid};
			$cat2 = $$shiu_cat{$pid} if defined $$shiu_aln{$pid};
			print $out_pks ">$pid PlantsP:$cat1;Shiu:$cat2\n$seq_info{$pid}{'seq'}\n";
		}
		$out_pks->close;

		print $report_info."\n";
		# remove temp folder
		# run_cmd("rm -rf $temp_dir") if -s $temp_dir;
	}
}

=head2
 itak_database -- prepare itak database
=cut
sub itak_database
{
	my ($options, $files) = @_;

        my $usage =qq'
USAGE:  perl $0 -t database  ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam27.0/Pfam-A.hmm.gz 

';
	print $usage and exit unless $$files[0];

	# check file exist
	my $bin_dir = ${FindBin::RealBin}."/bin";
	my $dbs_dir = ${FindBin::RealBin}."/database";
	unless (-e $bin_dir) { die "[ERR]bin folder not exist.\n$bin_dir\n"; }
	unless (-e $dbs_dir) { die "[ERR]database folder not exist.\n $dbs_dir\n"; }

	my $TF_selfbuild = $dbs_dir."/TF_selfbuild.hmm";	# selfbuild hmm
	my $plantsp_db	 = $dbs_dir."/PlantsPHMM3_89.hmm";	# plantsP kinase
	my $shiu_db	 = $dbs_dir."/Plant_Pkinase_fam.hmm";	# shiu kinase database
	my $Psub_wnk1	 = $dbs_dir."/Pkinase_sub_WNK1.hmm";	# wnk1 hmm
	my $Psub_MAK 	 = $dbs_dir."/Pkinase_sub_MAK.hmm";	# MAK hmm
	my $hmmpress_bin = $bin_dir."/hmmpress";		# hmmspress

	foreach my $f (($TF_selfbuild,$plantsp_db,$shiu_db,$Psub_wnk1,$Psub_MAK,$hmmpress_bin)) {
		die "[ERR]file not exist: $f\n" unless -s $f;
	}

	my $pfam_db = $dbs_dir."/TFHMM_3.hmm";                  # database for transcription factors (Pfam-A + customized)
	warn "[ERR]database exist $pfam_db\n" if -s $pfam_db;

	# build database;
	my $pfam_a_gz = $$files[0];	$pfam_a_gz =~ s/.*\///;
	my $pfam_a = $pfam_a_gz;	$pfam_a =~ s/\.gz//;
	my $hmmscan_bin = $bin_dir."/hmmscan";                  # hmmscan
	
	# run_cmd("wget $$files[0]");
	# run_cmd("gunzip $pfam_a_gz");
	# run_cmd("cat $pfam_a $TF_selfbuild > $pfam_db");
	# run_cmd("$hmmpress_bin $pfam_db");
	print qq'
wget $$files[0]
gunzip $pfam_a_gz
cat $pfam_a $TF_selfbuild > $pfam_db
rm $pfam_a

# hmmpress database
$hmmpress_bin $pfam_db
$hmmpress_bin $plantsp_db
$hmmpress_bin $shiu_db
$hmmpress_bin $Psub_wnk1
$hmmpress_bin $Psub_MAK

';
	exit;
}

#################################################################
# kentnf: true subroutine					#
#################################################################
=head2
 load_rule -- load rule to hash
 # this rule is update on 20141103
 # Description of each column
 # 1. ID of rule, the rule with small order number will have high priority
 # 2. name of the rule -- subfamily
 # 3. parent name of the rule -- superfamily
 # 4. required domain
 # 5. auxiiary domain
 # 6. forbidden domain
 # 7. type
 # 8. description
=cut
sub load_rule
{
	my $rule_file = shift;
	
	my %rule_obj;
	my ($id, $name, $family, $required, $auxiiary, $forbidden, $type, $desc);
	my $fh = IO::File->new($rule_file) || die $!;
	while(<$fh>)
	{
		chomp;
		next if $_ =~ m/^#/;
		
		if ($_ =~ m/^\/\//) 	# put current rule to hash, and start a new rule
		{
			die "[ERR]undef rule member $id\n" unless($id && $name && $family && $required && $auxiiary && $forbidden && $type && $desc);
			$rule_obj{$id}{'name'} = $name;
			$rule_obj{$id}{'family'} = $family;
			$rule_obj{$id}{'required'}  = parse_domain_rule($required);
			$rule_obj{$id}{'auxiiary'}  = parse_domain_rule($auxiiary);
			$rule_obj{$id}{'forbidden'} = parse_domain_rule($forbidden);
			$rule_obj{$id}{'type'} = $type;	
			$rule_obj{$id}{'desc'} = $desc;
			($id, $name, $family, $required, $auxiiary, $forbidden, $type, $desc) = ('', ''. '', '', '', '', '', '');
		} elsif ($_ =~ m/^ID:/) {
			$id = $_; $id =~ s/^ID://;
		} elsif ($_ =~ m/^Name:/) {
			$name = $_; $name =~ s/^Name://;
		} elsif ($_ =~ m/^Family:/) {
			$family = $_; $family =~ s/^Family://;
		} elsif ($_ =~ m/^Required:/) {
			$required = $_; $required =~ s/^Required://;
		} elsif ($_ =~ m/^Auxiiary:/) {
			$auxiiary = $_; $auxiiary =~ s/^Auxiiary://;
		} elsif ($_ =~ m/^Forbidden:/) {
			$forbidden = $_; $forbidden =~ s/^Forbidden://;
		} elsif ($_ =~ m/^Type:/) {
			$type = $_; $type =~ s/^Type://;
		} elsif ($_ =~ m/^Desc:/) {
			$desc = $_; $desc =~ s/^Desc://;
		} else {
			next;
		}
	}
	$fh->close;
	return %rule_obj;
}

=head
 parse_domain_rule: parse domain rules

 Description of domain rules
 PF00001 : domain ID, without version 
 #2 : number of requred domains
 '-' : and
 ':' : or members
 ';' : or rules
 '()' : priority

 example 1: PF00001#2-PF00002#1
 mean this rule require 2 of PF00001 domain and one of PF00002

 example 2; PF00001#2;PF00002#1
 mean this rule require 2 of PF00001 domain or one of PF00002

 example 3; PF00001#2;PF00002#1-PF00003#1:PF00004#1
 requred hash: 	1. PF00001 and PF00001
		2. PF00002 and PF00003
		3. PF00002 and PF00004
=cut
sub parse_domain_rule 
{
	my $domain_rule = shift;

	return $domain_rule if $domain_rule eq 'NA';

	my %domain_combination = ();			# key, array of domains for the rule.
	my @r = split(/;/, $domain_rule);
	foreach my $r ( @r ) {

		# hash for sub domain combination
		# key, array of domains for the rule, sub of domain_combination, equal to domain_combination when @r == 1;
		my %domain_combination_sub = ();

		my @m = split(/--/, $r);
		foreach my $m ( @m ) {

			my @p = split(/:/, $m);
			die "[ERR] $m\n" if scalar(@p) < 1;

			# for the first domain combination sub
			if (scalar keys %domain_combination_sub == 0) {
				foreach my $p ( @p ) {
					my $domain_id = split_domain_num($p);
					$domain_combination_sub{$domain_id} = 1;
				}
				next;
			}

			# for the single domain in this member
			if (scalar @p == 1) {
				my $domain_id = split_domain_num($p[0]);
				foreach my $com (sort keys %domain_combination_sub) {
					delete $domain_combination_sub{$com};		# remove old record
					$com.= ",$domain_id";				# add new domain to old for new record
					$domain_combination_sub{$com} = 1;		# put new reacord to hash
				}
				next;
			}

			# for the multiply domains in this member
			if (scalar @p > 1) {
				foreach my $com (sort keys %domain_combination_sub) {
					delete $domain_combination_sub{$com};		# remove old record
					foreach my $p (@p) {
						my $domain_id = split_domain_num($p);	
						my $new_com = $com.",$domain_id";	# add new domain to old for new record
						$domain_combination_sub{$new_com} = 1;	# put new reacord to hash
					}
				}
			}
		}

		# put sub domain combination to domain combination
		foreach my $com (sort keys %domain_combination_sub) {
			$domain_combination{$com} = 1;
		}
	}
	return \%domain_combination;
}

sub split_domain_num
{
	my $domain_num = shift;
	# print "x:$domain_num\n";

	die "[ERR]domain num format 1 $domain_num\n" unless $domain_num =~ m/#/;
	my @a = split(/#/, $domain_num);
	die "[ERR]domain num format 2 $domain_num\n" unless (scalar @a == 2);
	die "[ERR]domain num format 3 $domain_num\n" unless $a[1] > 0;
	my $domain_id = '';
	for (my $i=0; $i<$a[1]; $i++) {
		$domain_id.=",".$a[0];
	}
	$domain_id =~ s/^,//;
	return $domain_id;
}

=head2
 print_rule -- print rule in hash for debug
=cut
sub print_rule
{
	my $rule_pack = shift;
	my %rule = %$rule_pack;
	foreach my $id (sort keys %rule) {
	        print $id."\n";
	        print $rule{$id}{'name'},"\n";
	        print $rule{$id}{'family'},"\n";
	        print $rule{$id}{'type'},"\n";
	        print $rule{$id}{'desc'},"\n";

	        print "Required:\n";
	        foreach my $d (sort keys %{$rule{$id}{'required'}}) {
	                print $d."\n";
	        }

	        print "Auxiiary:\n";
		foreach my $d (sort keys %{$rule{$id}{'auxiiary'}}) {
	                print $d."\n";
	        }

		print "Forbidden:\n";
        	foreach my $d (sort keys %{$rule{$id}{'forbidden'}}) {
                	print $d."\n";
        	}
	}	
	exit;
}

=head2
 load_ga_cutoff: load GA cutoff score to hash
=cut
sub load_ga_cutoff 
{
	my ($pfam_db, $correct_ga, $sfam_db) = @_;


	# put GA cutoff to hash
	# key: pfam ID 
	# value: GA score
	my %ga_cutoff;
	my ($pfam_id, $ga_score) = ('', '');

	my $fh0 = IO::File->new($sfam_db) || die $!;
	while(<$fh0>)
	{
		chomp;
		if ($_ =~ m/^ACC\s+(\S+)/) {
			$pfam_id = $1;
			$pfam_id =~ s/\..*//;
		} elsif ($_ =~ m/^GA\s+(\S+)\s+(\S+)/) {	# using score for domain : The order of the thresholds is sequence, domain
			$ga_score = $2;
		} elsif ($_ eq "//") {
			warn "[WARN]no pfam id\n" unless $pfam_id;
			warn "[WARN]no ga score $pfam_id\n" unless $ga_score;
			$ga_cutoff{$pfam_id} = $ga_score;
			$pfam_id = '';
			$ga_score = '';
		} else {
			next;
		}
	}
	$fh0->close;

	my $fh1 = IO::File->new($pfam_db) || die $!;
	while(<$fh1>)
	{
		chomp;
		if ($_ =~ m/^ACC\s+(\S+)/) {
			$pfam_id = $1;
			$pfam_id =~ s/\..*//;
		} elsif ($_ =~ m/^GA\s+(\S+)/) {
			$ga_score = $1;
		} elsif ($_ eq "//") {
			warn "[WARN]no pfam id\n" unless $pfam_id;
			warn "[WARN]no ga score $pfam_id\n" unless $ga_score;
			$ga_cutoff{$pfam_id} = $ga_score;
			$pfam_id = '';
			$ga_score = '';
		} else {
			next;
		}	
	}
	$fh1->close;

	# correct GA cutoff
	my $stop = 0;
	my $fh2 = IO::File->new($correct_ga) || die $!;
	while(<$fh2>)
	{
		chomp;
		next if $_ =~ m/^#/;
		my @a = split(/\t/, $_);
		($pfam_id, $ga_score) = @a;
		$pfam_id =~ s/\..*//;
		warn "[ERR]no correct pfam id  $_\n" and $stop = 1 unless defined $ga_cutoff{$pfam_id};
		warn "[ERR]no correct ga score $_\n" and $stop = 1 unless $ga_score;
		$ga_cutoff{$pfam_id} = $ga_score;
	}
	$fh2->close;

	exit if $stop == 1;
	# print scalar(keys(%ga_cutoff)). "record has GA score\n";
	return %ga_cutoff;
}

=head2
 itak_tf_identify -- identify transcription factors

AT1G01140.1   PF00069.20      19      274     1       260     YEMGRTLGEGSFAKVKYAKNTVTGDQAAIKILDREKVFRHKMVEQLKrEISTMKLIKHPNVVEIIEVMASKTKIYIVLEL
VNGGELFDKIAQQGRLKEDEARRYFQQLINAVDYCHSRGVYHRDLKPENLILDANGVLKVSDFGLSAFSrqVREDGLLHTACGTPNYVAPEVLSDKGYDGAAADVWSCGVILFVLMAGYLPFDEP---NLMTLYKRICKAEFSC
PPWFS----QGAKRVIKRILEPNPITRISIAELLEDEWF ye +++lG+Gsf+kV  ak+  tg++ A+Kil++e+  + k  ++l+ E++ +k ++Hpn+v+++ev+ +k+++y+vle+v+gg+lfd ++++g+l+e+e++++
++q++++++y+Hs+g++HrDLKpeN++ld +g+lk++DFGl+      ++++ l+t +gt++Y+APEvl ++ ++++++DvWs+Gvil+ l+ g lpf++    + + l+++i k + + + + s    + +k +ik++le++p
 +R++++e+l+++w+ yelleklGsGsfGkVykakekktgkkvAvKilkkeeekskkektalr.ElkilkklnHpnivkllevfeekdelylvleyveggdlfdllkkkgklseeeikkialqilegleylHsngiiHrDLKpe
NiLldekgelkiaDFGlakkl..eksseklttlvgtreYmAPEvllkakeytkkvDvWslGvilyelltgklpfsgeseedqlelirkilkkkleedepkssskseelkdlikkllekdpakRltaeeilkhpwl 241.8  6
.4e-72  Protein kinase domain   448

=cut
sub itak_tf_identify 
{
	my ($hmmscan_hit, $hmmscan_detail, $ga_cutoff, $tf_rule) = @_;

	chomp($hmmscan_hit);
	chomp($hmmscan_detail);

	# create hash for result
	# key: query id of protein
	# value: tid of TF
	my %qid_tid;

	# put hits domains to hash : query_hits
	# key: query ID, score
	# value: array of hmm hits, and score
	# * filter hits with lower score using GA cutoff
	my %query_hits;
	my ($query_id, $pfam_id, $score, $evalue);
	my @a = split(/\n/, $hmmscan_hit);
	foreach my $a (@a) 
	{
		my @b = split(/\t/, $a);
		#AT1G01140.1     PF00069.20      241.8   6.4e-72
		($query_id, $pfam_id, $score, $evalue) = @b;
		$pfam_id =~ s/\..*// if $pfam_id =~ m/^PF/;
		die "[ERR]undef GA score for $pfam_id\n" unless defined $$ga_cutoff{$pfam_id};
		next if $score < $$ga_cutoff{$pfam_id};

		if (defined $query_hits{$b[0]}) {
			$query_hits{$b[0]}{'pid'}.="\t".$pfam_id;
			$query_hits{$b[0]}{'score'}.="\t".$score;
		} else {
			$query_hits{$b[0]}{'pid'} = $pfam_id;
			$query_hits{$b[0]}{'score'} = $score;
		}
	}

	# compare query_hits with rules, 
	foreach my $qid (sort keys %query_hits)
	{
		my $hits = $query_hits{$qid}{'pid'};		# do not use hash for hit two same pfam domain
		my $score = $query_hits{$qid}{'score'};	
		# print "$qid\t$hits\n";
		my $rule_id = compare_rule($hits, $score, $tf_rule);
		$qid_tid{$qid} = $rule_id if $rule_id ne 'NA';
	}

	return %qid_tid;
}

=head2
 compare_rule : compare
 # input is filtered hmm_hit and packed rules
=cut
sub compare_rule
{
	my ($hmm_hit, $hmm_score, $rule_pack) = @_;
	
	my %rule = %$rule_pack; # unpack the rule

	# compare the hits with rules, including required, auxiiary, and forbidden domains
	# the comparison will return match status: 
	# 0, do not match
	# 1, partially match
	# 2, full match
	# the assign family to each protein according hits and ruls
	my $rule_id = 'NA';

	my @hits = split(/\t/, $hmm_hit);
	my @score = split(/\t/, $hmm_score);
	die "[ERR]hit num do not mach score num\n" unless (scalar @hits == scalar @score);

	my $total_domain = 0;   # number of required domains used in classification
	my $total_score = 0;    # total of score of required domains
	foreach my $rid (sort keys %rule) {
		my $required_h = $rule{$rid}{'required'};
		my $auxiiary_h = $rule{$rid}{'auxiiary'};
		my $forbidden_h = $rule{$rid}{'forbidden'};

		# compare forbidden with hits
		my $f_status = 0;
		if ($forbidden_h ne 'NA') {
			foreach my $forbidden (sort keys %$forbidden_h) {
				my @f = split(/,/, $forbidden);
				my ($match_status, $match_score) = compare_array(\@hits, \@score, \@f);
				$f_status = 1 if $match_status > 0;
			}
		}
		next if $f_status == 1;

		# compare required with hits
		my $r_status = 0;
		foreach my $required (sort keys %$required_h) {
			my @r = split(/,/, $required);	
			my ($match_status, $match_score) = compare_array(\@hits, \@score, \@r);
			my $domain_num = scalar(@r);
			$r_status = 1 if $match_status == 2;
			if ($match_status == 2) {

				print "$rid\t$domain_num\t$match_score\n";

				if ($domain_num > $total_domain && $match_score > $total_score) {
					$total_domain = $domain_num;
					$total_score  = $match_score;
					$rule_id = $rid;
				}
			}
		}

		# compare auxiiary with hits, may use it in the future
		# my $a_status = 0;
		# if ($auxiiary_h ne 'NA') {
		#	foreach my $auxiiary (sort keys %$auxiiary_h) {
		#		my @a = split(/,/, $auxiiary);
		#		my $match_status = compare_array(\@hits, \@a);
		#		$a_status = 1 if $match_status == 2;
		#	}
		#	push(@rule_id, $rid) if ($r_status == 1 && $a_status == 1);
		#
		# } else {
		#	push(@rule_id, $rid) if $r_status == 1;
		# }
	}
	return $rule_id;
}

# sub for compare rule
sub compare_array
{
	my ($array_A, $score, $array_B) = @_;
	my @a = @$array_A;
	my @b = @$array_B;
	my @s = @$score;

	# convert array to hash
	# %ua: key: pfam_id, value: num
	# %ub: key: required_pfam_id, value: num
	# %sa: key: pfam_id, value: score array
	my %ua; my %ub; my %sa;
	foreach my $a (@a) {
		if (defined $ua{$a}) { $ua{$a}++; } else { $ua{$a} = 1; }
		my $score = shift @s;
		if (defined $sa{$a}) { 
			$sa{$a}.="\t".$score;
		} else {
			$sa{$a} = $score;
		}
	}

	foreach my $b (@b) {
		if (defined $ub{$b}) { $ub{$b}++; } else { $ub{$b} = 1; }
	}
	
	# compare two hash
	# find the best score for match
	my $match = 0;
	my $match_score = 0;
	foreach my $d (sort keys %ub) {
		if (defined $ua{$d} && $ua{$d} >= $ub{$d}) {
			$match++;
			my @s = split(/\t/, $sa{$d});
			my $n = $ub{$d};
			foreach my $s (sort {$b<=>$a} @s) {
				$match_score = $match_score + $s;
				$n--;
				last if $n == 0;		
			}
		}
	}

	# retrun match status
	# 0, do not match
	# 1, partially match
	# 2, full match	
	my $match_status = 0;
	$match_status = 1 if $match > 0;
	if ( $match == scalar(keys(%ub)) ) {
		$match_status = 2;
	}

	return ($match_status, $match_score); 
}

=head2
 itak_tf_write_out: write out tf result to output file
=cut
sub itak_tf_write_out
{
	my ($qid_tid, $seq_info, $hmmscan_detail_1, $tf_rule, $output_sequence, $output_alignment, $output_classification) = @_;

	# put hmmscan_detail to hash
	my %q_detail;
	chomp($hmmscan_detail_1);
	my @a = split(/\n/, $hmmscan_detail_1);
	foreach my $a (@a) {
		my @b = split(/\t/, $a);
		if (defined $q_detail{$b[0]}) {
			$q_detail{$b[0]}.= $a."\n";
		} else {
			$q_detail{$b[0]} = $a."\n";
		}
	}

	my $out1 = IO::File->new(">".$output_sequence) || die $!;
	my $out2 = IO::File->new(">".$output_alignment) || die $!;
	my $out3 = IO::File->new(">".$output_classification) || die $!; 

	foreach my $qid (sort keys %$qid_tid) {
		
		my $tid     = $$qid_tid{$qid};
		my $tname   = $$tf_rule{$tid}{'name'};
		my $tfamily = $$tf_rule{$tid}{'family'};
		my $type    = $$tf_rule{$tid}{'type'};
		my $desc    = $$tf_rule{$tid}{'desc'};
		my $qseq    = $$seq_info{$qid}{'seq'};
		my $align   = $q_detail{$qid};

		print $out1 ">$qid [$type]$tname:$tfamily--$desc\n$qseq\n";
		print $out2 $align;
		print $out3 "$qid\t$tname\t$type\t$tfamily\n";
	}
	$out1->close;
	$out2->close;
	$out3->close;
}

=head2 
 itak_pk_classify -- classify pkinase
=cut
sub itak_pk_classify
{
        my $hmmscan_detail = shift;
	chomp($hmmscan_detail);
	my @hit_line = split(/\n/, $hmmscan_detail);

	# put hmmscan hit to hash
	# %hit: key: protein seq ID 
	#       value: family id
	# %score: key: protein seq ID
	# 	  value: best hit score for domain
	# %best_hit_align: key: protein seq ID
	# 		   value: best hit alignment
        my %hit; my %score; my %best_hit_align;
	foreach my $line (@hit_line) {
		my @a = split(/\t/, $line);
		if (defined $hit{$a[0]}) {
			if ($a[9] > $score{$a[0]}) {
				$hit{$a[0]} = $a[1];
				$score{$a[0]} = $a[9];
				$best_hit_align{$a[0]} = $line."\n";
			}
		} else {
			$hit{$a[0]} = $a[1];
			$score{$a[0]} = $a[9];
			$best_hit_align{$a[0]} = $line."\n";
			
		}
	}
        return (\%hit, \%best_hit_align);
}

=head2 
 aln_to_hash: put hmmscan alignment detail to hash
=cut
sub aln_to_hash
{
        my ($hmmscan_detail, $ga_cutoff) = @_;

        my %aln_hash; # key: protein seq id; value: alignment information
	chomp($hmmscan_detail);
	my @detail_line = split(/\n/, $hmmscan_detail);
	foreach my $line ( @detail_line ) {
		my @a = split(/\t/, $line);
		my ($pfam_id, $score, $evalue) = ($a[1], $a[9], $a[10]);
		$pfam_id =~ s/\..*//;
		next if (defined $$ga_cutoff{$pfam_id} && $score < $$ga_cutoff{$pfam_id});
		next if ((!defined $$ga_cutoff{$pfam_id}) && ($evalue > 1e-3));

		if ( defined $aln_hash{$a[0]} ) {
			$aln_hash{$a[0]}.= $line."\n";
		} else {
			$aln_hash{$a[0]} = $line."\n";
		}
	}
        return %aln_hash;
}

=head2
 pk_to_hash -- load plantsp family description to hash
=cut
sub pk_to_hash
{
        my $file = shift;
        my %hash;
        my $pfh = IO::File->new($file) || die "Can not open protein kinase description file: $file $!\n";
        while(<$pfh>)
        {
                chomp;
                my @pm = split(/\t/, $_, 2);
                $hash{$pm[0]} = $pm[1];
        }
        $pfh->close;
        return (\%hash);
}

=head2
 seq_to_hash: put seq info to hash
=cut
sub seq_to_hash
{
	my $input_file = shift;
	my %seq_info;
	my $in = Bio::SeqIO->new(-format=>'fasta', -file=>$input_file);
	while(my $inseq = $in->next_seq)
	{
		$seq_info{$inseq->id}{'alphabet'} = $inseq->alphabet;
		$seq_info{$inseq->id}{'seq'} = $inseq->seq;
	}
	return %seq_info	
}

=head2
 run_cmd : run command
=cut
sub run_cmd
{
	my $cmd = shift;
	print "[ERR]no command: $cmd\n" and exit unless $cmd;
	print $cmd."\n" and return(1) if $debug;
	system($cmd) && die "[ERR]cmd: $cmd\n";
}

=head2
 usage : print usage
=cut
sub usage
{
	my $version = shift;
	my $usage = qq'
VERSION: $version
USAGE:  perl $0 -t [tool] 
	
	database	prepare database files for identification
	identify	identify TFs and PKs

';
	print $usage;
	exit;
}

