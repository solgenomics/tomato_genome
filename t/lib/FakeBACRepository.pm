package FakeBACRepository;
use Moose;
use FindBin;
use Path::Class;
use File::Temp;
use File::Copy ();

with 'CXGN::TomatoGenome::BACPublisher';

has '+bac_repository_root' => (
    required => 0,
    isa => 'File::Temp::Dir',
    default => sub {  File::Temp->newdir },
   );

sub cp {
    File::Copy::copy(map "$_", @_) or die "cannot copy @_";
}

sub BUILD {
    my $self = shift;
    my $tempdir = dir( $self->bac_repository_root );
    my $tf = 'C06HBa0012B01.1.v87.seq';
    my $data = dir( 't','data' );
    my $test_file = $data->file($tf);
    my $dest_file = $tempdir->subdir('chr06')->subdir('finished')->file($tf);
    $dest_file->dir->mkpath;
    cp( $test_file, $dest_file );
    cp( $test_file, $tempdir->file('bacs.v3453.seq') );
    cp( $data->file('foo2.gff3'), $tempdir->file('bacs.v3435.all.gff3') );
}


1;
