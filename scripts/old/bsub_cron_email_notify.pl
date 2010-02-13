#!/usr/bin/env perl
use strict;
use warnings;
use English;
use Carp;
use FindBin;
use Getopt::Std;

#use Data::Dumper;

use Digest::MD5 qw/md5_hex/;
use File::Basename;

use Mail::Sendmail;

use CXGN::TomatoGenome::BACPublish qw/find_submissions/;
use CXGN::TomatoGenome::BACSubmission;
use CXGN::TomatoGenome::Config;
use CXGN::IndexedLog;

####### DEFAULTS ####

my $default_logging_table = CXGN::TomatoGenome::Config->load_locked->{'bac_processing_log'};
# minimum time in seconds between successive emails about the same thing
my $notification_interval = 60*60*24*7; #< 1 week

#####################

sub usage {
  my $message = shift || '';
  $message = "Error: $message\n" if $message;

  require CXGN::TomatoGenome::Config;
  my $country_uploads_dir = CXGN::TomatoGenome::Config->load_locked->{'country_uploads_path'};

  die <<EOU;
$message
Usage:
  $FindBin::Script

  BAC pipeline cron script to notify BAC submitters of problems with
  either new submissions or current submissions.

  Options:

    -d <dir>
      (optional) set root of uploading file hierarchy
      Defaults to: $country_uploads_dir

    -l <name>
      (optional) set name of database table to use for logging.
      Defaults to: $default_logging_table

    -x
      don't actually send emails, just print what you would send to
      standard output
EOU
}
sub HELP_MESSAGE {usage()}

my %opt;
getopts('xd:l:',\%opt) or usage();

# for each thing we check, log whether we've emailed about it, and
# email again if the last time we emailed about it was more than the
# given maximum time ago
my $dbh = CXGN::DB::Connection->new;
my $log = CXGN::IndexedLog->open( DB => $dbh,
				   ($opt{l} || $default_logging_table)
				 );

my @submission_files = find_submissions($opt{d});

my %messages_to_send;
# hash of {  email address => { submitter => submitter record,
#                               errors => [ [ submission file, error text ], ...],
#                               log_strings => [ str, str, ...],
#                             },
#            ...
#         }

# go through the new submissions, check for errors
foreach my $subfile (@submission_files) {
  my ($bn,$uldir) = fileparse($subfile);

  my $sub = eval {CXGN::TomatoGenome::BACSubmission->open_stripped($subfile)};
  my $errors = $EVAL_ERROR ? 'Error untarring submission file.  Please re-upload.' 
    : $sub->validation_text;

  next unless $errors;

  print "found errors for $bn:\n$errors\n" if $opt{x};

  my @submitters = $sub->submitters;

  # check whether we've already sent this email recently, and skip this file if we have
  foreach my $submitter (@submitters) {
    my $errs_checksum = md5_hex($errors);
    my $log_string = "EMAIL_NOTIFIED_ABOUT_${errs_checksum} $subfile";
    my %log_entry = $log->lookup(content => $log_string);
    if( %log_entry && $log_entry{timestamp} >= time()-$notification_interval ) {
      # if we sent an email exactly like this up to a week ago, don't send another one
      warn "skipping email to $submitter->{name} about $bn\n";
      next;
    }

    my $messages_record = $messages_to_send{$submitter->{email}} ||= {};
    $messages_record->{submitter} = $submitter;
    push @{$messages_record->{errors}}, [ $subfile, $errors ];
    push @{$messages_record->{log_strings}}, $log_string;
  }
}


# now send all the emails we need to send
foreach my $message (values %messages_to_send) {
  my $errors_text = join "\n", map {
    my ($bn,$dir) = fileparse($_->[0]);
    $_->[1]
  } @{$message->{errors}};

  my $message_text = <<EOT;
Greetings, this is the SGN BAC submission pipeline, informing you of problems with BAC submissions for which you are registered as the contact person.

$errors_text

Please fix the problems detailed above, either by re-uploading the submissions that have intrinsic problems, or by updating external data sources (like GenBank/EMBL/DDBJ, or the SGN BAC Registry) to conform to the data in the submissions you uploaded, or by canceling the submission (use SCP to delete the uploaded file).  If you have questions, or need assistance, please feel free to email sgn-feedback\@sgn.cornell.edu.  We are here to help!

This message will repeat periodically until the problems are fixed or submissions are canceled.

Sincerely,
The SGN BAC Submission Pipeline Robot
EOT
  print (("="x50)."\n"."email to $message->{submitter}{name} <$message->{submitter}{email}>:\n$message_text");

  unless( $opt{x} ) { 
    sendmail( From => 'cxgnbacpipeline@upload.sgn.cornell.edu',
	      To => "$message->{submitter}->{name} <$message->{submitter}->{email}>",
	      Subject => 'submission problems - SGN BAC submission pipeline',
	      Body => $message_text,
	      Cc  => 'Robert Buels <rmb32@cornell.edu>',
	      'Reply-To' => 'sgn-feedback@sgn.cornell.edu',
	    )
      or die "sendmail failed: $Mail::Sendmail::error, $Mail::Sendmail::log\n";

    foreach my $log_string (@{$message->{log_strings}}) {
      $log->append($log_string);
    }
  }
}


