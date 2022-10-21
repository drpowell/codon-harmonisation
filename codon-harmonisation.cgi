#!/usr/bin/perl -Tw

use strict;

use CGI;
use CGI::Carp qw(fatalsToBrowser warningsToBrowser);

my $css=<<FOO;
th {
  background-color: yellow;
}

.f {
  color: red;
}

.t {
  color: blue;
}
FOO

my $q = new CGI;

print $q->header;
print $q->start_html(-title => "Codon Harmonisation Tool",
                     -style => {-code => $css });

print "<h1>Codon Harmonisation Tool</h1>\n";

print "Codon usage statistics from <a href=http://www.kazusa.or.jp/codon/>Codon Usage Database</a>\n";

print "<hr>\n";

my($dat, $genetic_code, $rev_genetic_code, $protein_abbrev) = read_data("data.spsum");

my $from = $q->param('from');
my $to   = $q->param('to');
my $seq  = $q->param('seq');
my $verbose = $q->param('verbose');
my $amino_norm = $q->param('amino_norm');

$verbose ||= 0;

if (defined($from) && defined($to) && defined($seq)) {
  process_input($from, $to, $seq, $verbose);
}


display_form();

print "<em>Author: <a href=mailto:david.powell\@drp.id.au>David Powell</a></em>\n";

warningsToBrowser(1);           # Any warnings go in the html as comments
print $q->end_html;


sub process_input {
  my($from, $to, $seq, $verbose) = @_;
  
  if (!exists($dat->{$from})) {
    print "<h2>Error: '$from' is not a known organism</h2>\n";
    return;
  }

  if (!exists($dat->{$to})) {
    print "<h2>Error: '$to' is not a known organism</h2>\n";
    return;
  }

  $seq =~ tr/uU/tT/;
  $seq =~ s/\s+//g;

  if ($seq !~ /^[atgc]*$/i) {
    print "<h2>Error: Bad character in sequence</h2>\n";
    return;
  }

  if (length($seq)%3 != 0) {
    print "<h2>Error: Sequence is not an integer number of codons</h2>\n";
    return;
  }



  my $mapping = codon_mapping($from, $to);
  process_sequence($from, $to, $mapping, $seq, $verbose);

  display_mapping($from, $to, $mapping) if ($verbose >= 1);

  display_codon_usage($from, $to)       if ($verbose >= 1);
}


sub process_sequence {
  my($from_species, $to_species, $mapping, $seq, $verbose) = @_;

  my $orig = '';
  my $conv = '';

  my $protein = '';

  print "<h2>Codon harmonisation for <span class=f>'$from_species'</span> to <span class=t>'$to_species'</span></h2>\n";

  my @codons;
  for(my $i=0; $i<length($seq); $i+=3) {
    my $codon = uc substr($seq, $i, 3);

    $orig .= $codon . " ";
    $conv .= sprintf "%s%s%s ",
               $codon eq $mapping->{$codon} ? "" : "<b>",
               $mapping->{$codon},
               $codon eq $mapping->{$codon} ? "" : "</b>";
    
    $protein .= $protein_abbrev->{$genetic_code->{$codon}};

#    if (($i+3) % 30 == 0) {
#      $orig .= "<br>";
#      $conv .= "<br>";
#    }
  }

  print "<h3>Original Sequence</h3><tt>$orig</tt>";
  print "<h3>Converted Sequence</h3><tt>$conv</tt>";
  print "<h4>Protein Sequence</h4><tt>$protein</tt>";

  print "<hr>\n";
}

sub codon_mapping {
  my($from_species, $to_species) = @_;
  
  my $f_dat =  $dat->{$from_species};
  my $t_dat =  $dat->{$to_species};
  my %mapping;

  for my $amino (sort amino_sort keys %$rev_genetic_code) {
    for my $from_codon ( @{$rev_genetic_code->{$amino}} ) {
      my $f_norm = ($amino_norm ? $f_dat->{AMINO_NORM}{$amino} : $f_dat->{TOTAL});
      my $t_norm = ($amino_norm ? $t_dat->{AMINO_NORM}{$amino} : $t_dat->{TOTAL});

      my $val = $f_dat->{USAGE}{$from_codon}/$f_norm;

      # Find closest in t_dat for $val
      my $best_codon;
      my $best_diff;
      for my $to_codon ( @{$rev_genetic_code->{$amino}} ) {
        my $to_val = $t_dat->{USAGE}{$to_codon}/$t_norm;
        if (!defined($best_codon) ||
            $best_diff > abs($val - $to_val)) {
          $best_codon = $to_codon;
          $best_diff = abs($val - $to_val);
        }
      }

      $mapping{$from_codon} = $best_codon;
    }
  }

  return \%mapping;
}

sub display_mapping {
  my($from_species, $to_species, $mapping) = @_;

  print "<h3>Codon harmonisation mapping from <span class=f>'$from_species'</span> to <span class=t>'$to_species'</span></h3>\n";
  print "<table border=1>";
  for my $amino (sort amino_sort keys %$rev_genetic_code) {
    my $pre_str = "<tr><th>$amino";
    for my $codon ( @{$rev_genetic_code->{$amino}} ) {
      if ($codon ne $mapping->{$codon}) {
        printf "%s<td><span class=f>%s</span> -> <span class=t>%s</span>\n",
               $pre_str, $codon, $mapping->{$codon};
        $pre_str = "";
      }
    }
  }
  print "</table>";
  print "<hr>\n";

}

sub amino_sort {
  ($a eq 'END') && return 1;
  ($b eq 'END') && return -1;
  $a cmp $b;
}

sub display_codon_usage {
  my($from_species, $to_species) = @_;

  my $f_dat =  $dat->{$from_species};
  my $t_dat =  $dat->{$to_species};

  print "<h3>Codon usage for <span class=f>'$from_species'</span> and <span class=t>'$to_species'</span></h3>\n";

  print "<table border=1>";
  for my $amino (sort amino_sort keys %$rev_genetic_code) {
    printf "<tr><th>%s (%s)",$amino,$protein_abbrev->{$amino};

    printf "<th><table><tr><td><span class=f>%.2f%%</span><tr><td><span class=t>%.2f%%</span></table>\n",
           $f_dat->{AMINO_NORM}{$amino}/$f_dat->{TOTAL} * 100,
           $t_dat->{AMINO_NORM}{$amino}/$t_dat->{TOTAL} * 100;

    my $f_norm = ($amino_norm ? $f_dat->{AMINO_NORM}{$amino} : $f_dat->{TOTAL});
    my $t_norm = ($amino_norm ? $t_dat->{AMINO_NORM}{$amino} : $t_dat->{TOTAL});
    
    for my $codon ( @{$rev_genetic_code->{$amino}} ) {
      print "<td>$codon<td>";
      print "<table>";
      printf "<tr><td><span class=f>%.2f</span>", 
              $f_dat->{USAGE}{$codon}/$f_norm * 100;
      printf "<tr><td><span class=t>%.2f</span>",
              $t_dat->{USAGE}{$codon}/$t_norm * 100;
      print "</table>\n";
    }
  }
  print "</table>\n";

  print "<hr>";
}

sub display_form {
  print $q->start_form;
  print "<TABLE>\n";
  print "<TR><TD>From Organism: ";
  print "<TD>" .$q->popup_menu(-name => 'from', -values => [sort keys %$dat],
                               -default => 'Homo sapiens');
  print "<TR><TD>Target Organism: ";
  print "<TD>" . $q->popup_menu(-name => 'to', -values => [sort keys %$dat],
                                -default => 'Escherichia coli');
  print "<TR><TD>Raw DNA sequence: ";
  print "<TD>" . $q->textarea(-name => 'seq', -rows => 10, -columns => 50);
  print "<TR><TD>Verbose: ";
  print "<TD>" . $q->popup_menu(-name => 'verbose', -values => [0,1], 
                                -labels => {0 => 'No', 1 => 'Yes'});
  print "<TR><TD>Normalise";
  print "<TD>" . $q->popup_menu(-name => 'amino_norm', -values => [1,0], 
                                -labels => {1 => 'Over amino acids', 0 => 'Over everything'});
  print "</TABLE>\n";
  print $q->submit, "\n";
  print $q->end_form, "\n";

  print "<hr>\n";
}

sub read_data {
  my($fname) = @_;

  my($code_by_num, $code_by_codon, $genetic_code, $protein_abbrev) = gen_codon_code();

  # Reverse genetic code.  ie. make amino -> list of codons
  my $rev_genetic = {};
  for my $codon (keys %$genetic_code) {
    $rev_genetic->{$genetic_code->{$codon}} ||= [];
    push(@{ $rev_genetic->{$genetic_code->{$codon}} }, $codon);
  }


  my $res = {};

  open(F, "< $fname") or die "Can't read $fname";
  while(<F>) {
    next unless (/(.*?):\s(\d+)$/);
    my($species, $numSeq) = ($1,$2);
    my $counts = <F>;
    my @counts = split /\s+/, $counts;

    if (@counts != 64) {
      print STDERR "Bad number of counts for '$species'\n";
      next;
    }

    my $sum = 0;
    map { $sum += $_ } @counts;
    map { $_ = $_/$sum } @counts;      # Convert to fractional

    my $usage= {};
    for my $i (0 .. $#counts) {
     $usage->{$code_by_num->[$i]} =  $counts[$i];
    }

    $res->{$species}{NUMSEQ} = $numSeq;
    $res->{$species}{USAGE} = $usage;
   
    # Calc totals per amino acid

    my $totalPerAmino = {};
    for my $amino (keys %$rev_genetic) {
      my $sum = 0;
      for my $codon ( @{$rev_genetic->{$amino}} ) {
        $sum += $usage->{$codon};
      }
      $totalPerAmino->{$amino} = $sum;
    }

    # Calculate relative usage of each amino acid.
    my $totalNum = 0;
    map { $totalNum += $_ } values %$totalPerAmino;
    $res->{$species}{TOTAL} = $totalNum;


    # Set the AMINO normalising value
    $res->{$species}{AMINO_NORM} = $totalPerAmino;

  }

  return ($res, $genetic_code, $rev_genetic, $protein_abbrev);
}

sub gen_codon_code {
  my $code_by_num   = [];
  my $code_by_codon = {};
  my $genetic_code  = {};

  my $abbrev = {};

  while(<DATA>) {
    next if /^\s*$/;
    if (/^(...): (.)$/) {
      $abbrev->{$1} = $2;
      next;
    }

    die "Internal error: '$_'" unless (/^\s*(\d+): (...) (...)/);
    my($num, $codon, $amino) = ($1,$2,$3);
    $genetic_code->{$codon}  = $amino;
    $code_by_num->[$num]     = $codon;
    $code_by_codon->{$codon} = $num;
  }
  return ($code_by_num, $code_by_codon, $genetic_code, $abbrev);
}

__DATA__
Arg: R
Leu: L
Ser: S
Thr: T
Pro: P
Ala: A
Gly: G
Val: V
Lys: K
Asn: N
Gln: Q
His: H
Glu: E
Asp: D
Tyr: Y
Cys: C
Phe: F
Ile: I
Met: M
Trp: W
END: *
 0: CGA Arg
 1: CGC Arg
 2: CGG Arg
 3: CGT Arg
 4: AGA Arg
 5: AGG Arg
 6: CTA Leu
 7: CTC Leu
 8: CTG Leu
 9: CTT Leu
10: TTA Leu
11: TTG Leu
12: TCA Ser
13: TCC Ser
14: TCG Ser
15: TCT Ser
16: AGC Ser
17: AGT Ser
18: ACA Thr
19: ACC Thr
20: ACG Thr
21: ACT Thr
22: CCA Pro
23: CCC Pro
24: CCG Pro
25: CCT Pro
26: GCA Ala
27: GCC Ala
28: GCG Ala
29: GCT Ala
30: GGA Gly
31: GGC Gly
32: GGG Gly
33: GGT Gly
34: GTA Val
35: GTC Val
36: GTG Val
37: GTT Val
38: AAA Lys
39: AAG Lys
40: AAC Asn
41: AAT Asn
42: CAA Gln
43: CAG Gln
44: CAC His
45: CAT His
46: GAA Glu
47: GAG Glu
48: GAC Asp
49: GAT Asp
50: TAC Tyr
51: TAT Tyr
52: TGC Cys
53: TGT Cys
54: TTC Phe
55: TTT Phe
56: ATA Ile
57: ATC Ile
58: ATT Ile
59: ATG Met
60: TGG Trp
61: TAA END
62: TAG END
63: TGA END

