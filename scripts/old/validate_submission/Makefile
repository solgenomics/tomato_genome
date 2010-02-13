SOURCEFILES := ../bsub_validate_submission.pl ../../../perllib/CXGN/Tools/List.pm ../../../perllib/CXGN/Genomic/CloneIdentifiers.pm ../../../perllib/CXGN/TomatoGenome/BACSubmission.pm

script: temp_patchgen
	mv temp_patchgen validate_submission.pl
	chmod a+x validate_submission.pl
	patch -p0 < patch

temp_patchgen: $(SOURCEFILES)
	perl -pe '$$line = $$_; $$skipping = 1 if $$line =~ /BEGIN_SKIP_FOR_PORTABLE_VALIDATION/; $$_="" if $$skipping; $$skipping = 0 if $$line =~ /END_SKIP_FOR_PORTABLE_VALIDATION/;' $(SOURCEFILES) > temp_patchgen

patch:  temp_patchgen validate_submission.pl
	- diff -u temp_patchgen validate_submission.pl > patch
	rm -f temp_patchgen

install: script test
	cp validate_submission.pl /data/prod/public/tomato_genome/bacs/validate_submission.pl

test:
	prove /usr/local/lib/site_perl/t/CXGN/bacsubmission.t
	perl ./validate_submission.pl ~rob/bacs/chr01/*/C01HBa0163B20.1.*.tar.gz
clean:
	rm -f validate_submission.pl

.PHONY: test
