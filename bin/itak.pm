#!/usr/bin/perl
package itak;

use strict;
use warnings;
use Exporter;
use IO::File;

#my ($all_hits, $all_detail) = parse_hmmscan_result($tmp_hmm_result);

=head2 parse_hmmscan_result
 
 Function: parse hmmscan result 

 Input: hmmscan result file name

 Output1: detail hits info, format is below:
	GeneID      PfamID		GA	    Evalue
	AT2G34620.1	PF02536.7	242.2	4.5e-72

 Output2: alignment of hits into, format is below:
	1. GeneID      -- AT1G01140.1
	2. PfamID      -- PF00069.18
	3. Query Start -- 19
	4. Query End   -- 274
	5. Hit Start   -- 1
	6. Hit End     -- 260
	7. Query Seq   -- YEMGRTLGEGSFAKVKYAKNTVTGDQAAIKILDREKVF....
	8. Alignment   -- ye++++lG+Gsf+kV  ak+  tg++ A+Kil++e+  ....
	9. Hit Seq     -- yelleklGsGsfGkVykakkkktgkkvAvKilkkeeek....
	10. GA Score   -- 241.4
	11. Evalue     -- 6.9e-72 
	12. Description-- Protein kinase domain
	13. Qeury Len  -- 448

 Output3: for some hits, the best hits domain may not meet the requirement of significant
          but multi-hits for one domain actually meet the GA cutoff. That why the GA score
          has two values. One for overall, another for best domain
	GeneID      PfamID      GA      Evalue
    AT2G34620.1 PF02536.7   242.2   4.5e-72
=cut
sub parse_hmmscan_result
{
	my $hmm_result = shift;

	my ($result_out1, $result_out2, $result_out3, $one_result, $align_detail, $hits);

	my $rfh = IO::File->new($hmm_result) || die "Can not open hmmscan result file : $hmm_result $! \n";
	while(<$rfh>)
	{
		unless (/^#/)
		{
			if (/^\/\//)
			{
				#################################
				# parse domain info		#
				#################################

				############################################################
				# short info for every hsp				   #
				# $hits = query_id \t domain_id \t evalue \t score desc \n #
				############################################################
	
				#$hits = parse_Hits($one_result);
				my @hits_content = split(/>>/, $one_result);

				#########################################################
				# parse hit head content like below			#
				#########################################################
				#
				# Match a line like this
				# E-value  score  bias    E-value  score  bias    exp  N   Sequence Description
				# -------  -----  -----   -------  ------ -----   ---- --  -------- -----------
				#  4e-83   285.8  10.0    5.3e-83  285.5  7.0     1.1  1   Q14SN3.1 Q14SN3_9HEPC Polyprotein (Fragment).
				#######################################################################################################
				# previous hmmer3 parse function need this part for best one domain.
				# New version just need Query name, length and no hit information in this part.
				#######################################################################################################
				my @hit_head = split(/\n/, $hits_content[0]);

				my ($query_name, $query_length);

				my $jumper = 0;
				foreach my $hit_head_line (@hit_head)
				{
					if ($hit_head_line =~ m/^Query:\s+(\S+)\s+\[M\=(\d+)\]/) 
					{
						$query_name = $1; $query_length = $2;
					}
					elsif ($hit_head_line =~ m/^Query:\s+(\S+)\s+\[L\=(\d+)\]/) 
					{
						$query_name = $1; $query_length = $2;
					}
					elsif ($hit_head_line =~ m/No hits detected that satisfy reporting thresholds/)
					{
						$jumper = 1;
						#die "$hit_head_line\n ";
					}
					else 
					{
						my @b = split(/\s+/, $hit_head_line, 11);
						next unless defined $b[2];
						next unless ($b[2] =~ m/^\d+\.\d+/);
						$result_out3.="$query_name\t$b[9]\t$b[2]\t$b[1]\n";
						#print "$query_name\t$b[9]\t$b[2]\t$b[1]\n";
					}
				}

				#########################################################
				# parse hsp part of hits				#
				#########################################################

				my $hsp_detail = ""; my $hsp_info = "";

				unless( $jumper == 1)
				{
					for(my $ih=1; $ih<@hits_content; $ih++)
					{
						my $one_hit = ">>".$hits_content[$ih];
						($hsp_info, $hsp_detail) = parse_align($one_hit, $query_name, $query_length);
						$result_out1.= $hsp_info;
						$result_out2.= $hsp_detail;
					}
				}

				#########################################################
				# init 							#
				#########################################################
			    	$one_result = "";
			}
			else 
			{
				#store all one protein hmmscan info to this char;
				$one_result.=$_;
			}
		}
	}
	$rfh->close;
	return ($result_out1, $result_out2, $result_out3);
}

sub parse_align
{
	my ($hsp_info, $query_name, $query_length) = @_;

	my $output1 = ""; my $output2 = "";
	
	#########################################################
	# get hit id, hit desc and hsp info form one hit	#
	#########################################################
	my @hsp_line = split(/\n/, $hsp_info);

	my ($hit_id, $hit_desc);
	my %info1 = ();

	for(my $i=0; $i<@hsp_line; $i++)
	{
		if ($hsp_line[$i] =~ m/^>>\s+/)
		{
			my @aa = split(/\s+/, $hsp_line[$i], 3);
			$hit_id = $aa[1];
			$hit_desc = $aa[2];
		}
		elsif ( $hsp_line[$i] =~ m/^\s+(\d+)\s+\W\s+/)
		{
			$info1{$1} = $hsp_line[$i];
		}
		else
		{
			next;	
		}
	}

	#########################################################
	# get query string, hit string, match string of HSP	#
	#########################################################
	my ($query_string, $hit_string, $match_string, $hsp_length, $align_pos, $match_start);

	my @domain = split(/== domain/, $hsp_info);

	for(my $j=1; $j<@domain; $j++)
	{
		my @domain_line = split(/\n/, $domain[$j]);

		my @info = split(/\s+/, $info1{$j});

		for(my $k=1; $k<@domain_line; $k++)
		{
			#################################################
			# get hit string of HSP				#
			#################################################
			if ($domain_line[$k] =~ m/^\s+\Q$hit_id\E\s+\d+\s+(\S+)\s+\d+/ )
			{
				my @fff = split(/\s+/, $domain_line[$k]);

				$hit_string = $fff[3];
	
				$hsp_length = length($hit_string);

				$align_pos = index($domain_line[$k], $fff[3]);

				$match_start = 1;

				if ($align_pos < 0)
				{
					die "Error! Align start position is negative: $align_pos\n$domain_line[2]\n$hit_string\n";
				}
			}

			#################################################
			# get match string of HSP			#
			#################################################
			elsif (defined $match_start && $match_start == 1 )
			{
				$match_string = substr($domain_line[$k], $align_pos, $hsp_length);
				$match_start = 0;
			}

			#################################################
			# get query string of HSP			#
			#################################################
			elsif ($domain_line[$k] =~ m/^\s+\Q$query_name\E/ ) 
			{
				my @qqq = split(/\s+/, $domain_line[$k]);
				$query_string = $qqq[3];
			}
			else
			{
				next;
			}

		}
		
		#where these code come from? 
		#$score   = $dMatch[3];
		#$evalue  = $dMatch[6];
		#$hmmfrom = $dMatch[7];
		#$hmmto	  = $dMatch[8];
		#$seqfrom = $dMatch[10];
		#$seqto   = $dMatch[11];

		$output1.="$query_name\t$hit_id\t$info[3]\t$info[6]\n";
		$output2.="$query_name\t$hit_id\t$info[10]\t$info[11]\t$info[7]\t$info[8]\t$query_string\t$match_string\t$hit_string\t$info[3]\t$info[6]\t$hit_desc\t$query_length\n";
	}

	return ($output1,$output2);
}

=head2 parse_format_result

 Function: filter hmmscan parsed result, if socre >= GA score or evalue <= 1e-3, it will be seleted.

 Input: 1. formated result of hmmscan: query name; hit name;score; evalue; 
	2. gathering score hash, PF id ; GA score; 
                                 key      vaule

 Return: hash1 -- key: seq_id,    value: domain_id1 \t domain_id2 \t ... \t domain_idn
 Return: hash2 -- key: domain_id, value: PfamID
 Retrun: hash3 -- key: domain_id, value: SeqID \t PfamID \t 
=cut
sub parse_format_result
{
	my ($in_file, $ga_score_hash) = @_;

	my %out_hash; my %hsp_hit_id; my %hsp_detail;

	my @in_file = split(/\n/, $in_file);

	my $len = length(scalar(@in_file)); 

	my $uid = 0;
	
	for(my $in=0; $in<@in_file; $in++)
	{
		#################################################
		# filter HSP detail result with GA or e-value	#
		#################################################
		chomp($in_file[$in]);
		my @fmm = split(/\t/, $in_file[$in]);

		$fmm[1] =~ s/\..*//; #this is Pfam or self-build domain ID

		my $valued = 0; my $gaScore;

		if (defined $$ga_score_hash{$fmm[1]})
		{
			$gaScore = $$ga_score_hash{$fmm[1]};

			if ( $fmm[9] >= $gaScore ) { $valued = 1; }
		}
		else
		{ 
			if ($fmm[10] <= 1e-3) { $valued = 1; }
		}

		#################################################
		# if valued means we select this domain info	#
		# then creat three hashes base on this		#
		#################################################
		if ($valued == 1)
		{
			$uid++; my $zero = "";
			my $rlen = $len-length($uid);
			for(my $l=0; $l<$rlen; $l++)
			{
				$zero.="0";
			}

			$hsp_detail{$zero.$uid} = $in_file[$in];
			$hsp_hit_id{$zero.$uid} = $fmm[1];
			
			if (defined $out_hash{$fmm[0]})
			{
				$out_hash{$fmm[0]}.= "\t".$zero.$uid;
			}
			else
			{
				$out_hash{$fmm[0]}.= $zero.$uid;
			}
		}
	}	
	return (\%out_hash, \%hsp_hit_id, \%hsp_detail);
}

1;
