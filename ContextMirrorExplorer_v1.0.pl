#!/usr/bin/perl

=head1 NAME

ContextMirrorExplorer_v1.0.pl - Retrieves co-evolutionary partners for a list of input genes from prokaryotic species

=head1 SYNOPSIS


Use:

    perl ContextMirrorExplorer_v1.0.pl [--help] [--man]
                               [--out <output file name>] [--out_type <R|N1|N2|N3|N4>]
                                --ids <list of gene_ids>

Examples:

    perl ContextMirrorExplorer_v1.0.pl --help

    perl ContextMirrorExplorer_v1.0.pl --man

    #   Load all RefSeq data:

    perl ContextMirrorExplorer_v1.0.pl --ids NP_414543.1,NP_414590.1,NP_414631.1 --out output_example --out_type N4

=head1 DESCRIPTION

This script has been created to provide an easy access to the results of the large-scale co-evolutionary analyses performed in the context of MICROME project (http://www.microme.eu/).

ContextMirror has been used for detecting co-evolutionary signals out of the evolutionary histories of the groups of orthologs detected in sets of selected species [1].

These sets have been selected in order to provide an adequate representation of the evolutionary information associated to a collection of 24 reference proteomes [2].

Although these analyses are focused on the reference species their results can be fully expanded to other species included in the analyses.

This script (ConextMirrorExplorer.pl) accepts a list of gene/protein identifiers (ENSEMBL, locus, GeneID, Genbank, EcoGene and gene names) and retrieves the list of predicted co-evolving partners for the corresponding query proteins (as pairs of locus tags).

These results can be retrieved as a report or as different combinations of network files suitable for being visualized and analysed with cytoscape software (http://www.cytoscape.org/).

This program uses orthology-based transference of the results obtained for the reference species to any query sequence from over 1,300 species.

It also combine the information provided by the different reference-centered analyses to provide an integrated set of predicted coevolving partners for each query protein.

References:

1.	Juan D, Pazos F, Valencia A (2008) High-confidence prediction of global interactomes based on genome-wide coevolutionary networks. Proc Natl Acad Sci USA 105: 934–939.
2.	Herman D, Ochoa D, Juan D, Lopez D, Valencia A, et al. (2011) Selection of organisms for the co-evolution-based study of protein interactions. BMC Bioinformatics 12: 363.

=head1 ARGUMENTS


ContextMirrorExplorer_v1.pl takes the following arguments:

=over 4

=item help

  --help

(Optional.) Displays the usage message.

=item man

  --man

(Optional.) Displays all documentation.

=item ids

  --ids

(Required.) Comma-separated list of gene identifiers (ENSEMBL, locus, GeneID, Genbank, EcoGene and gene names)

=item out

  --out

(Optional.) Prefix for all the output files (default value: output/output, it writes output files in "output" directory). Be careful it might overwrite previous results.

=item out_type

  --out_type

(Optional.) Specify the type of output (default value: R):

=over 8

=item R

Report output. Text file summarizing all the information of the results.

=item N1

Separated networks output. It generates network files for every pair of query protein and reference species.

=item N2

Separated reference species networks output. It generates network files for every reference species.

=item N3

Integrated query protein networks output. It generates network files for every query protein combining pairs from all the reference species (higher correlation values are retained).

=item N4

Integrated networks. It generates network files combining results from query proteins from the same species integrating results from all the reference species (as in N3).
				
=back

All the networks are generated as a combination of a file defining edges in .sif format (see http://wiki.cytoscape.org/Cytoscape_User_Manual/Network_Formats) 

and a file defining the correlation values as edge atributes in .eda format (see http://wiki.cytoscape.org/Cytoscape_User_Manual/Attributes).

These files can be easily parsed and directly imported in cytoscape software (http://www.cytoscape.org/) for further visualization and analysis.
				
=head1 AUTHOR

David Juan at the Spanish National Cancer Research Centre (CNIO), E<lt>dadejuan@cnio.es<gt>.

=head1 COPYRIGHT

This program is distributed under the Artistic License.

=head1 DATE

26-Apr-2013

=cut


package main;
use strict;
use warnings;
use IO::File;
use Getopt::Long ();    #   Resist name-space pollution!
use Pod::Usage ();      


MAIN:
{


# Working configuration
my $input   = "";
my $distance = 1;
my $ids_filename="data/ids.txt";
my $species_filename="data/prokaryotes2013.txt";
my $networks_dir="data/networks/";
my $out_prefix="output/output";
my $out_type="R";


#R  -> Report
#N1 -> Fully Separated Networks
#N2 -> Separated Networks by query species and reference species
#N3	-> Integrated Networks from reference species separated by query genes
#N4 -> Integrated Networks from reference species


# Internal variables
my ($help,$man,$line,$internal_id,$id,$net,$id1,$id2,$tmp_str,$sort_res,$id_res,$id_tmp,$in);
my (@input,@ids_file,@input_ids,@species_file,@tr,@tr2,@input_networks,@network_file,@output_networks,@sort_selected);
my (%ids_equiv,%fixed_ids_equiv,%species_input,%species_name,%species_taxid,%id_species,%selected_pairs,%input_output,%results,%try,%done,%integrated_results);


Getopt::Long::GetOptions ('help'          =>  \$help,
        																		'man'           =>  \$man,
        																		"out=s" 								=>  \$out_prefix,
        																		"out_type=s"    =>  \$out_type,
                        		"ids=s"   						=>  \$input);

Pod::Usage::pod2usage( -verbose => 1 ) if ( $help );
Pod::Usage::pod2usage( -exitstatus => 0, -verbose => 2 ) if ( $man  );

unless ( $input )
    {
        Pod::Usage::pod2usage( -exitstatus => 2 );
    }
$input=~s/\s+//ig;

@input_ids=split /\,/,$input;
foreach (@input_ids){if(!/[\w\_\-\.\/(\)]/){Pod::Usage::pod2usage( -exitstatus => 2 );}}


open IDS_FILE, "$ids_filename" or &error_handler('open_file',$ids_filename);
@ids_file=<IDS_FILE>;
close IDS_FILE;
foreach $line(@ids_file)
{
	foreach $id(@input_ids)
	{
		if($line=~/\b$id\b/)
		{
			chomp $line;
			@tr=split /\t/,$line;
			foreach (@tr)
			{
				push @{$fixed_ids_equiv{$_}},$tr[1];
			}
			
			$internal_id=$tr[1];
			@tr2=split /\_/,$internal_id;
			
			$id_species{$tr[1]}=$tr2[$#tr2];
			$species_input{$tr2[$#tr2]}=1;
			last;
		}
	}
}
undef @ids_file;
$in=0;
foreach $id(@input_ids)
	{
		if(!exists($fixed_ids_equiv{$id}))
		{
				print STDERR "WARNING: I coudn't find any information about $id\n";
		}else
		{
			$in=1;
		}
	}

if(!$in){Pod::Usage::pod2usage( -exitstatus => 2 );}
open SPECIES_FILE, "$species_filename" or &error_handler('open_file',$species_filename);
@species_file=<SPECIES_FILE>;
close SPECIES_FILE;
foreach $line(@species_file)
{
	chomp $line;
	@tr=split /\t/,$line;
	$species_name{$tr[3]}=$tr[0];
	$species_taxid{$tr[3]}=$tr[1];
}
undef @species_file;



@input_networks=keys(%species_input);
foreach (@input_networks)
{
	open NETWORK_FILE, "$networks_dir/$_\.r2" or &error_handler('open_file',"$networks_dir/$_\.r2");
	@network_file=<NETWORK_FILE>;
	close NETWORK_FILE;
	foreach $line(@network_file)
	{
		chomp;
		@tr=split /\t/,$line;
		if(exists($id_species{$tr[1]}))
		{
			$id1=$tr[1];
			$id1=~s/\_(\d+)$//;
			$id2=$tr[2];
			$id2=~s/\_(\d+)$//;
			$results{$tr[1]}{$tr[0]}{"$id1|$id2"}=$tr[4];
			if(!exists($integrated_results{$tr[1]}{"$tr[1]|$tr[2]"}) && !exists($integrated_results{$tr[1]}{"$tr[2]|$tr[1]"}))
			{
				$integrated_results{$tr[1]}{"$tr[1]|$tr[2]"}=$tr[4];
			}else
			{
				if(exists($integrated_results{$tr[1]}{"$tr[1]|$tr[2]"}))
				{
					if($tr[4]>$integrated_results{$tr[1]}{"$tr[1]|$tr[2]"})
					{
						$integrated_results{$tr[1]}{"$tr[1]|$tr[2]"}=$tr[4];
					}
				}else
				{
					if($tr[4]>$integrated_results{$tr[1]}{"$tr[2]|$tr[1]"})
					{
						$integrated_results{$tr[1]}{"$tr[2]|$tr[1]"}=$tr[4];
					}
				}
			}
			$input_output{$tr[1]}{"$tr[1]|$tr[2]"}=1;
		}
		if(exists($id_species{$tr[2]}))
		{
			$id1=$tr[1];
			$id1=~s/\_(\d+)$//;
			$id2=$tr[2];
			$id2=~s/\_(\d+)$//;
			$results{$tr[2]}{$tr[0]}{"$id1|$id2"}=$tr[4];
			#print"$tr[2]\t$tr[0]\n";
			if(!exists($integrated_results{$tr[2]}{"$tr[1]|$tr[2]"}) && !exists($integrated_results{$tr[2]}{"$tr[2]|$tr[1]"}))
			{
				$integrated_results{$tr[2]}{"$tr[1]|$tr[2]"}=$tr[4];
			}else
			{
				if(exists($integrated_results{$tr[2]}{"$tr[1]|$tr[2]"}))
				{
					if($tr[4]>$integrated_results{$tr[2]}{"$tr[1]|$tr[2]"})
					{
						$integrated_results{$tr[2]}{"$tr[1]|$tr[2]"}=$tr[4];
					}
				}else
				{
					if($tr[4]>$integrated_results{$tr[2]}{"$tr[2]|$tr[1]"})
					{
						$integrated_results{$tr[2]}{"$tr[2]|$tr[1]"}=$tr[4];
					}
				}
			}
			$input_output{$tr[2]}{"$tr[1]|$tr[2]"}=1;
		}
	}
}


if($out_type eq "R")
{
	$"="_";
	open OUT_FILE, ">$out_prefix\_@input_ids\_report.txt" or &error_handler('open_file',"$out_prefix\_@input_ids\_report.txt");
	$"="\t";
foreach $id(@input_ids)
{
	print OUT_FILE "\n\n# QUERY: $id\n\n";
	foreach $id_res(@{$fixed_ids_equiv{$id}})
	{
		@output_networks=keys(%{$results{$id_res}});

	#if(exists($results{$fixed_ids_equiv{$id}})){
	
	$id_tmp=$id_res;
	$id_tmp=~s/\_[0-9]+$//;
	print OUT_FILE "\n# Results for $id($id_tmp) from $species_name{$id_species{$id_res}} (TaxID: $species_taxid{$id_species{$id_res}})\n";
	foreach $net(@output_networks)
	{
		print OUT_FILE "\n# Reference Species: $species_name{$net} (TaxID: $species_taxid{$net})\n\n";
		@sort_selected=sort{$results{$id_res}{$net}{$b}<=>$results{$id_res}{$net}{$a}} (keys(%{$results{$id_res}{$net}}));
		foreach $sort_res(@sort_selected)
		{
			$tmp_str=$sort_res;
			$tmp_str=~s/\|/\t/ig;
			print OUT_FILE "\t$tmp_str\t$results{$id_res}{$net}{$sort_res}\n";
		}
	}
	}
	}
	close OUT_FILE;
}elsif($out_type =~ /^N[1-4]/)
{
	
	foreach $id(@input_ids)
	{
		foreach $id_res(@{$fixed_ids_equiv{$id}})
		{
			$id_tmp=$id_res;
			$id_tmp=~s/\_[0-9]+$//;
			if($out_type eq "N1" or $out_type eq "N2")
			{
				@output_networks=keys(%{$results{$id_res}});

				foreach $net(@output_networks)
				{
					if($out_type eq "N1")
					{
						open OUT_PAIRS_FILE, ">$out_prefix\_$id\_$id_tmp\_$species_taxid{$id_species{$id_res}}\_$species_taxid{$net}\.sif" or &error_handler('open_file',"$out_prefix\_$id\_$id_tmp\_$species_taxid{$id_species{$id_res}}\_$species_taxid{$net}\.sif");
						open OUT_ATTRIBUTES_FILE, ">$out_prefix\_$id\_$id_tmp\_$species_taxid{$id_species{$id_res}}\_$species_taxid{$net}\.eda" or &error_handler('open_file',"$out_prefix\_$id\_$id_tmp\_$species_taxid{$id_species{$id_res}}\_$species_taxid{$net}\.eda");
						print OUT_ATTRIBUTES_FILE "CoevolutionProfiles\n";
					}elsif($out_type eq "N2")
					{
						if(!exists($done{"$out_prefix\_$species_taxid{$id_species{$id_res}}\_$species_taxid{$net}"}))
						{
							$done{"$out_prefix\_$species_taxid{$id_species{$id_res}}\_$species_taxid{$net}"}=1;
							open OUT_PAIRS_FILE, ">$out_prefix\_$species_taxid{$id_species{$id_res}}\_$species_taxid{$net}\.sif" or &error_handler('open_file',"$out_prefix\_$species_taxid{$id_species{$id_res}}\_$species_taxid{$net}\.sif");
							open OUT_ATTRIBUTES_FILE, ">$out_prefix\_$species_taxid{$id_species{$id_res}}\_$species_taxid{$net}\.eda" or &error_handler('open_file',"$out_prefix\_$species_taxid{$id_species{$id_res}}\_$species_taxid{$net}\.eda");
							print OUT_ATTRIBUTES_FILE "CoevolutionProfiles\n";
						}else
						{
							open OUT_PAIRS_FILE, ">>$out_prefix\_$species_taxid{$id_species{$id_res}}\_$species_taxid{$net}\.sif" or &error_handler('open_file',"$out_prefix\_$species_taxid{$id_species{$id_res}}\_$species_taxid{$net}\.sif");
							open OUT_ATTRIBUTES_FILE, ">>$out_prefix\_$species_taxid{$id_species{$id_res}}\_$species_taxid{$net}\.eda" or &error_handler('open_file',"$out_prefix\_$species_taxid{$id_species{$id_res}}\_$species_taxid{$net}\.eda");
						}
					}
					
					@sort_selected=sort{$results{$id_res}{$net}{$b}<=>$results{$id_res}{$net}{$a}} (keys(%{$results{$id_res}{$net}}));
					foreach $sort_res(@sort_selected)
					{
						$tmp_str=$sort_res;
						$tmp_str=~s/\|/ coev /;
						print OUT_PAIRS_FILE "$tmp_str\n";
						$tmp_str=~s/ coev / (coev) /;
						print OUT_ATTRIBUTES_FILE "$tmp_str = $results{$id_res}{$net}{$sort_res}\n";
					}
					close OUT_PAIRS_FILE;
					close OUT_ATTRIBUTES_FILE;
				}
			}else
			{
				if($out_type eq "N3")
				{
						open OUT_PAIRS_FILE, ">$out_prefix\_$id\_$id_tmp\_$species_taxid{$id_species{$id_res}}\.sif" or &error_handler('open_file',"$out_prefix\_$id\_$id_tmp\_$species_taxid{$id_species{$id_res}}\.sif");
						open OUT_ATTRIBUTES_FILE, ">$out_prefix\_$id\_$id_tmp\_$species_taxid{$id_species{$id_res}}\.eda" or &error_handler('open_file',"$out_prefix\_$id\_$id_tmp\_$species_taxid{$id_species{$id_res}}\.eda");
						print OUT_ATTRIBUTES_FILE "CoevolutionProfiles\n";
				}elsif($out_type eq "N4")
				{
					if(!exists($done{"$out_prefix\_$species_taxid{$id_species{$id_res}}"}))
					{
						$done{"$out_prefix\_$species_taxid{$id_species{$id_res}}"}=1;
						open OUT_PAIRS_FILE, ">$out_prefix\_$species_taxid{$id_species{$id_res}}\.sif" or &error_handler('open_file',"$out_prefix\_$species_taxid{$id_species{$id_res}}\.sif");
						open OUT_ATTRIBUTES_FILE, ">$out_prefix\_$species_taxid{$id_species{$id_res}}\.eda" or &error_handler('open_file',"$out_prefix\_$species_taxid{$id_species{$id_res}}\.eda");
						print OUT_ATTRIBUTES_FILE "CoevolutionProfiles\n";
					}else
					{
						open OUT_PAIRS_FILE, ">>$out_prefix\_$species_taxid{$id_species{$id_res}}\.sif" or &error_handler('open_file',"$out_prefix\_$species_taxid{$id_species{$id_res}}\.sif");
						open OUT_ATTRIBUTES_FILE, ">>$out_prefix\_$species_taxid{$id_species{$id_res}}\.eda" or &error_handler('open_file',"$out_prefix\_$species_taxid{$id_species{$id_res}}\.eda");
					}
				}
				@sort_selected=sort{$integrated_results{$id_res}{$b}<=>$integrated_results{$id_res}{$a}} (keys(%{$integrated_results{$id_res}}));
				foreach $sort_res(@sort_selected)
					{
						$tmp_str=$sort_res;
						$tmp_str=~s/\|/ coev /;
						print OUT_PAIRS_FILE "$tmp_str\n";
						$tmp_str=~s/ coev / (coev) /;
						print OUT_ATTRIBUTES_FILE "$tmp_str = $integrated_results{$id_res}{$sort_res}\n";
					}
					close OUT_PAIRS_FILE;
					close OUT_ATTRIBUTES_FILE;
			}
		}
	}
	
}




sub error_handler
{
	my $error_type=$_[0];
	my $add_info=$_[1];
	if(!$error_type){die "ERROR: wrong call to error_handler()";}
	elsif($error_type eq 'open_file'){die	"ERROR: I couldn't open $add_info\n";}
}



exit 0;
}
