##1. Bio-CUA

This is a Perl software package for analyzing codon usage bias.

##2. Purpose

The aim of this distribution is to provide comprehensive and flexible
tools to analyze codon usage bias (CUB) and relevant problems.

One amino acid can be encoded by more than one synonymous codon, and
synonymous codons are unevenly used. For example, some codons are used
more often than other synonymous ones in highly expressed genes
(I<Sharp and Li 1987>). To measure the unevenness of codon usage, multiple
metrics of codon usage bias have been developed, such as Fop
(Frequency of optimal codons), CAI (Codon Adaptation Index), tAI (tRNA
Adaptation Index), and ENC (Effective Number of Codons). Biased codon
usage is widespread, visible in all species. It is important both to
identify codons having high translational efficiency (often termed
optimal codons) and to study the distribution of codon usage among
genes (e.g., genes with more optimal codons versus genes with fewer
optimal codons).

So far, no software implements all the above CUB metrics in one place.
More importantly, parameters in existing software are often fixed,
so one can compute certain types of CUB metrics for a limited list of 
species and cannot tune parameters. For example, when one wants to
identify optimal codons in certain tissues, it may be better to use
most highly expressed genes of that tissue to calculate CAI, which is
impossible with existing software.

This package mainly solves these two problems: providing tools
computing all common CUB metrics and allowing users to tune parameters
freely. We also incorporate or extend some method variants, such as 
GC-content corrected ENC, background-data normalized CAI, etc. 
See the relevant methods' description in CUB classes for more details.


##3. INSTALLATION

To install this module, run the following commands:

a. If install after downloading the package,

	perl Makefile.PL
	make
	make test
	make install

b. install directly from CPAN as

	cpan Bio::CUA

##4. SUPPORT AND DOCUMENTATION

After installing, you can find documentation for this module with the
perldoc command.

    perldoc Bio::CUA

You can also find a tutorial on the CPAN, which is also accessible
through

    perldoc Bio::CUA::Tutorial
	
You can also look for information at:

    RT, CPAN's request tracker (report bugs here)
        http://rt.cpan.org/NoAuth/Bugs.html?Dist=Bio-CUA

    AnnoCPAN, Annotated CPAN documentation
        http://annocpan.org/dist/Bio-CUA

    CPAN Ratings
        http://cpanratings.perl.org/d/Bio-CUA

    Search CPAN
        http://search.cpan.org/dist/Bio-CUA/


You can also email me at zhangz.sci@gmail.com for help.

##5. LICENSE AND COPYRIGHT

Copyright (C) 2015 Zhenguo Zhang

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see L<http://www.gnu.org/licenses/>.

##6. CITATION

Zhenguo Zhang, CUA: a Flexible and
Comprehensive Codon Usage Analyzer (In preparation)

