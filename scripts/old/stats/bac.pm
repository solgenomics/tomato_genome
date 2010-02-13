
package bac;

return 1;

sub new { 
    my $class = shift;
    my $id = shift;
    my $args = {};
    my $self = bless $args, $class;
    $self->set_id($id);
    $self->set_has_T7(0);
    $self->set_has_SP6(0);
    @{$self->{T7}} = ();
    @{$self->{SP6}} = ();
    return $self;
}

sub set_id { 
    my $self = shift;
    $self -> {name}=shift;
}

sub get_id { 
    my $self = shift;
    return $self->{name};
}

sub set_has_T7 { 
    my $self = shift;
    $self->{has_end}->{T7} = shift;
}

sub set_has_SP6 {
    my $self = shift;
    $self->{has_end}->{SP6} = shift;
}

sub set_has_end { 
    my $self = shift;
    my $end = shift;
    $self->{has_end}->{$end}=1;
}

sub has_T7 { 
    my $self = shift;
    return $self->{has_end}->{T7};
}

sub has_SP6 { 
    my $self = shift;
    return $self ->{has_end}->{SP6};
}

sub has_end {
    my $self = shift;
    my $end = shift;
    return $self->{has_end}->{$end};
}

sub set_seq { 
    my $self = shift;
    my $dir = shift;
    $self->{seq}{$dir}=shift;
}

sub get_seq { 
    my $self = shift;
    my $dir = shift;
    return $self->{seq}{$dir};
}

sub clear_T7_annotations { 
    my $self = shift;
    @{$self->{annotations}->{T7}}=();
}

sub clear_SP6_annotations { 
    my $self=shift;
    @{$self->{annotations}->{SP6}}=();
}

sub clear_annotations { 
    my $self = shift;
    my $end = shift;
    @{$self->{annotations}->{$end}}=();
}

sub add_T7_annotation { 
    my $self = shift;
    my $annotation = shift;
    push @{$self->{annotations}->{T7}}, $annotation;
}

sub add_SP6_annotation { 
    my $self = shift;
    my $annotation = shift;
    push @{$self ->{annotations}->{SP6}}, $annotation;
}

sub add_annotation { 
    my $self = shift;
    my $end = shift;
    push @{$self->{annotations}->{$end}}, $annotation;
}

sub add_unknown_annotation { 
    my $self = shift;
    my $annotation = shift;
    push @{$self->{annotations}->{unknown}}, $annotation;
}

sub get_T7_annotations { 
    my $self = shift;
    return @{$self->{annotations}->{T7}};
}

sub get_SP6_annotations { 
    my $self = shift;
    return @{$self->{annotations}->{SP6}};
}

# sub get_annotations { 
#     my $self = shift;
#     my $end = shift;
#     if (!$end) { print STDERR "get_annotations: which end?\n"; exit(0); }
#     return @{$self->{annotations}->{$end}}; 
# }

sub set_T7_class { 
    my $self = shift;
    $self->{end_class}->{T7} = shift;
}

sub set_SP6_class { 
    my $self = shift;
    $self->{end_class}->{SP6} = shift;
}

sub set_end_class { 
    my $self = shift;
    my $end = shift;
    $self->{end_class}->{$end} = shift;
}

sub get_SP6_class { 
    my $self = shift;
    return $self->{end_class}->{SP6};
}

sub get_T7_class { 
    my $self = shift;
    return $self->{end_class}->{T7};
}

sub get_end_class {
    my $self = shift;
    my $end = shift;
    if (!$end) { print STDERR "get_end_class: which end?\n"; exit(); }
    return $self->{end_class}->{$end};
}

=head2 function get_T7_blast_coords

Synopsis:	
Arguments:	
Returns:	
Side effects:	
Description:	

=cut

sub get_T7_blast_coords { 
    my $self=shift;
    return join ",", @{$self->{blast_coords}->{T7}};
}

=head2 function set_T7_blast

Synopsis:	
Arguments:	
Returns:	
Side effects:	
Description:	

=cut

sub add_T7_blast_coords { 
    my $self=shift;
    my $s_start = shift;
    my $s_end = shift;
    push @{$self->{blast_coords}->{T7}}, "$s_start..$s_end";
}

=head2 function get_SP6_blast

Synopsis:	
Arguments:	
Returns:	
Side effects:	
Description:	

=cut

sub get_SP6_blast_coords { 
    my $self=shift;
    return join ",", @{$self->{blast_coords}->{SP6}};
}

sub get_blast_coords { 
    my $self = shift;
    my $end = shift;
    return join ",", @{$self->{blast_coords}->{$end}};
}

=head2 function set_SP6_blast_coords

Synopsis:	
Arguments:	
Returns:	
Side effects:	
Description:	

=cut
    
    sub add_SP6_blast_coords { 
	my $self=shift;
	my $s_start = shift;
	my $s_end = shift;
	push @{$self->{blast_coords}->{SP6}}, "$s_start..$s_end";
    }

sub add_blast_coords { 
    my $self = shift;
    my $direction = shift;
    my $s_start = shift;
    my $s_end = shift;
    
    push @{$self->{blast_coords}->{$direction}}, "$s_start..$s_end"; 
}


sub clear_blast_coords { 
    my $self =shift;
    my $direction = shift;
    if ($direction) { 
	@{$self->{blast_coords}->{$direction}}=(); 
    }
    else { print STDERR "clear_blast_coords: which end?\n"; exit(0); }
}


=head2 function get_bac_class

Synopsis:	
Arguments:	
Returns:	
Side effects:	
Description:	

=cut

sub get_bac_class { 
    $self=shift;
    return $self->{bac_class};
}

=head2 function set_bac_class

Synopsis:	
Arguments:	
Returns:	
Side effects:	
Description:	

=cut

sub set_bac_class { 
    $self=shift;
    $self->{bac_class}=shift;
}



sub get_end_seq_count { 
    my $self = shift;
    if ($self->has_T7() && $self->has_SP6()) { 
	return 2;
    }
    elsif ($self->has_SP6() || $self->has_T7()) { 
	return 1;
    }
    else { return 0; }
}
    
sub add_annotation { 
    my $self = shift;
    my $direction = shift;
    my $annotation = shift;
    if (!$direction =~ /T7|SP6/i) { $self->add_unknown_anntotation($annotation); }
    if ($direction=~/T7/i) { $self->add_T7_annotation($annotation); }
    if ($direction=~/SP6/i) { $self->add_SP6_annotation($annotation); }
    
}

sub set_annotation { 
    my $self = shift;
    my $direction = shift;
    if (!$direction) { print STDERR "Need direction.\n"; exit(0); }
    $self->{annotation}{$direction} = shift;
    $self->{evalue}{$direction} = shift;
    $self->{bitscore}{$direction} = shift;
}

sub get_annotation { 
    my $self = shift;
    my $direction = shift;
    if (!$direction) { print STDERR "Need direction.\n"; exit(0); }
    return ($self->{annotation}{$direction}, $self->{evalue}{$direction},$self->{bitscore}{$direction});
}
#
=head2 function get_end_class_path

Synopsis:	
Arguments:	
Returns:	
Side effects:	
Description:	

=cut
#
sub get_end_class_path { 
    my $self=shift;
    my $direction = shift;
    return $self->{end_class_path}{$direction};
}

#
=head2 function set_end_class_path

Synopsis:	
Arguments:	
Returns:	
Side effects:	
Description:	

=cut
#   

sub set_end_class_path { 
    my $self=shift;
    my $direction = shift;
    $self->{end_class_path}{$direction}=shift;
}

   

sub set_chromosome { 
    my $self = shift;
    $self->{chromosome} = shift;
}

sub get_chromosome{ 
    my $self  = shift;
    return $self->{chromosome};
}

sub set_position { 
    my $self = shift;
    $self->{position} = shift;
}

sub get_position { 
    my $self = shift;
    return $self ->{position};
}

sub get_GC_content { 
    my $self = shift;
    my $end = shift;
    if (!$end) { print STDERR "get_GC_content: which end?\n"; exit(0); }
    return $self->{GC}->{$end};
    
}

sub set_GC_content { 
    my $self =shift;
    my $end = shift;
    if (!$end) { print STDERR  "set_GC_content: which end?\n"; exit(0); }
    $self->{GC}->{$end}=shift;
}

sub set_intrinsic_repeat { 
    my $self = shift;
    my $direction = shift;
    $self->{intrinsic_repeat}{$direction} = shift;
}

sub get_intrinsic_repeat { 
    my $self = shift;
    my $direction = shift;
    return $self->{intrinsic_repeat}{$direction};
}
