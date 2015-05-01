use 5.006;
use strict;
use warnings;
use ExtUtils::MakeMaker;

WriteMakefile(
    NAME             => 'Bio::CUA',
    AUTHOR           => q{Zhenguo Zhang <zhangz.sci@gmail.com>},
    VERSION_FROM     => 'lib/Bio/CUA.pm',
    ABSTRACT_FROM    => 'lib/Bio/CUA.pm',
    LICENSE          => 'gpl_3',
	#PL_FILES         => {'other.PL' => ''}, # files run during compilation 
	EXE_FILES        => ["bin/build_cai_param.pl",
	                     "bin/build_tai_param.pl",
						 "bin/calculate_CUB.pl"
	                    ], # files in bin directory
    MIN_PERL_VERSION => 5.006,
    CONFIGURE_REQUIRES => {
        'ExtUtils::MakeMaker' => 0,
    },
    BUILD_REQUIRES => {
        'Test::More' => 0,
    },
    PREREQ_PM => {
        #'ABC'              => 1.6,
        #'Foo::Bar::Module' => 5.0401,
		'File::Sort'        => 1.01,
		'File::Temp'        => 0,
		'Fcntl'             => 0,
		'parent'            => 0,
		'Getopt::Long'      => 0,
		'Carp'              => 0
    },

	# add more information
	META_MERGE  => {
		resources => {
			homepage => 'http://molevol.altervista.org/'
		}
	},

    dist  => { COMPRESS => 'gzip -9f', SUFFIX => 'gz', },
    clean => { FILES => 'Bio-CUA-*' },
	test  => { RECURSIVE_TEST_FILES => 1, # test directories in t recursively
		# TESTS => 't/01-seqio.t'
	},
	# environment variables
	macro => {
#		'expected_file' => 't/expected.tsv' # expected results to compare
	}
);

package MY; # other content MakeMaker can't produce

sub postamble
{
	# The return will be directly write to makefile
	my $postamble = <<POST;
export expected_file = t/expected.tsv
export seq_file = t/test.fa
export tRNA_file = t/dmel_r5_apr2006.tRNA_copy
POST

	return $postamble;
}

# put subroutines here
