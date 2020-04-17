'''

Changelog:
    -0.1 initial version
    -0.2 remove k-mer connections that are too far away from the mean
    -0.3 bugfix with erroneous connections
    -0.31 fix typo
'''
from argparse import ArgumentParser
import sys
import logging
import glob
from Bio import SeqIO, Seq
from contig import Contig
from kmer import Kmer

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(levelname)-8s %(message)s',
                    datefmt='%d.%m %H:%M:%S')
logger = logging.getLogger()

SKIPSIZE = 400            # how many k-mers we skip
UNPLACE_CUTOFF = 5        # if a contig has only this many hits we don't place it
OVERLAP_CUTOFF = 80       # if a contig pair has less than this shared unique k-mers then the pair is ignored
VERSION = 0.31


def get_contigs_recursively(contig, set_of_contigs, contigs_to_ignore):
    ''' Takes a contig and keeps on iterating over the connections of that contil until nothing left

    Takes: a contig Object
           a set of contig names (String)
           a set of contig names to ignore
    returns: an updated set of contig names (String)
    '''

    this_contig = dict_of_contigs[contig.name]
    this_contig_connections = this_contig.connections
    if not this_contig_connections:
        # no connections left, return
        return set_of_contigs
    for c in this_contig_connections:
        for other_c in this_contig_connections[c]:
            this_hit, this_contig = other_c
            # c = 60231 
            # other_c = (4584, Object:tig00000001XXXuniq_kmer_id_chr1A.soap.74.contigs.fasta)
            if tuple(sorted([contig.name, other_c[1].name])) in contigs_to_ignore:
                continue
            if this_contig.name not in set_of_contigs:
                set_of_contigs.add(this_contig.name)
                # we found a new contig, get the connections of that contig
                return get_contigs_recursively(this_contig, set_of_contigs, contigs_to_ignore)
    # we have connections but we have seen them all before
    return set_of_contigs


parser = ArgumentParser(description='This script takes unique k-mers in tabular formats and a bunch of fastas and stitches them based on shared k-mers.')
parser.add_argument('file_of_kmers', help='File detailing the paths to all k-mer files, in order')
parser.add_argument('file_of_fastas', help='File detailing paths to all fasta files, in same order')
parser.add_argument('txt_out', help='Output file - txt file which has the order of scaffolds and notes')
parser.add_argument('fasta_out', help='Output file - fasta file which has merged and reversed contigs')

# optional arguments here
parser.add_argument('--skipsize', default=SKIPSIZE,
                    help='Which n-th k-mer is checked. Set to 0 to skip nothing. Default: %s'%SKIPSIZE)
parser.add_argument('--min_placement', default=UNPLACE_CUTOFF,
                    help='How many k-mers have to hit a contig for that contig to be placed. Default: %s'%UNPLACE_CUTOFF)
parser.add_argument('--min_pair_overlap', default=OVERLAP_CUTOFF,
                    help='Contig pairs with less than this many shared k-mers are ignored. Default: %s'%(OVERLAP_CUTOFF))
parser.add_argument("-v", "--verbose", help="Increase output verbosity",
                    action="store_true")
parser.add_argument('-V','--version', help='Print version number and exit.',
                    action='version', version='Version: %s'%VERSION)

args = parser.parse_args()

if args.verbose:
    logger.setLevel(logging.DEBUG)

fasta_out = open(args.fasta_out, 'w')
txt_out = open(args.txt_out, 'w')
logging.info('Writing to %s and %s'%(fasta_out.name, txt_out.name))

skipsize = int(args.skipsize)
min_placement = int(args.min_placement)
overlap_cutoff = int(args.min_pair_overlap)

kmers = open(args.file_of_kmers)
# Parse the name of all fasta files into a list
fastas = open(args.file_of_fastas)
index_of_fastas = {}
list_of_fastas = []
for fastacounter, f in enumerate(fastas):
    f = f.rstrip()
    list_of_fastas.append(f)

fastas.close()

#  iterate over all files of k-mers - then iterate over the three contigs surrounding this

first_kmer = '' # we store the kmers as a linked list of kmers
previous_kmer = ''

dict_of_contigs = {} # key: contig name, value: object of class Contig

logging.info('Going over all kmer files now...')
for counter, k in enumerate(kmers):
    # k = uniq_kmer_id_chr1A.soap.88
    k = k.rstrip()
    logging.info('Looking at k-mers in %s and getting positions in surrounding fastas'%k)
    fh = open(k)
    inside_counter = 0
    
    # now we get the surrounding fastas and see where this k-mer hits
    fastas = []
    try:
        fastas.append(list_of_fastas[counter-1])
    except IndexError:
        # first fasta
        pass

    fastas.append(list_of_fastas[counter])

    try:
        fastas.append(list_of_fastas[counter+1])
    except IndexError:
        # end of fastas
        pass

    # for NOW, I'm only getting the 3 fastas around.
 
    for linecounter, line in enumerate(fh):
        if skipsize > 0 and inside_counter / skipsize != 1:
            logging.debug('skipping counter %s with skipsize %s'%(linecounter, skipsize))
            inside_counter += 1
            continue
        inside_counter = 1

        ll = line.split()
        # ['896171689', 'ACTCCAACATGTCTCACAAGCCTGCGTCAGCCCCAAAGACC', '+', 'chr1A', '239482626']
        name, sequence, chrom = ll[0], ll[1], ll[3]
        this_kmer = Kmer(name, sequence, chrom)
        rev_this_kmer = str(Seq.Seq(sequence).reverse_complement())
        if not first_kmer:
            # initialise the linked list
            first_kmer = this_kmer

        if previous_kmer:
            # link them up!
            previous_kmer.set_next_kmer(this_kmer)
            this_kmer.set_previous_kmer(previous_kmer)
       
        # THE FOLLOWING WOULD BE WAY FASTER IF WE'D STORE THE CONTIGS IN AN FM-INDEX
        # BUT THERE IS NO NICE PYTHON MODULE AND I'M LAZY

        for fasta in fastas:
            sequences_seen = 0 # SeqIO does not complain if an input file is not a fasta
            for seq in SeqIO.parse(fasta, 'fasta'):
                # we have to have a unique name, combo of contig name + filename should work
                long_name = '%sXXX%s'%(seq.id, fasta) 
                sequences_seen += 1
                if long_name in dict_of_contigs:
                    this_contig = dict_of_contigs[long_name]
                else:
                    # make the new Contig object
                    name = seq.id
                    filename = fasta
                    sequence = str(seq.seq)
                    this_contig = Contig(name, filename, sequence)
                    dict_of_contigs[long_name] = this_contig

                sequence = this_contig.sequence
                # search for our k-mer in this sequence
                try:
                    this_hit = sequence.index(this_kmer.seq)
                except ValueError:
                    this_hit = ''
                # try the reverse complement too!
                try:
                    rev_hit = sequence.index(rev_this_kmer)
                except ValueError:
                    rev_hit = ''

                # this_hit is now the position of this_kmer on the current contig
                # store the hit in the contig's datastructure
                # and store the hit in the k-mers data structure
                if rev_hit:
                    this_contig.revs_fors['rev'] += 1
                    this_kmer.set_of_hits.add( (rev_hit, this_contig) )
                if this_hit:
                    this_contig.revs_fors['for'] += 1
                    this_kmer.set_of_hits.add( (this_hit, this_contig) )
                this_contig.kmer_hits[this_kmer.name] = this_kmer
            if not sequences_seen:
                logging.error('You do not have any entries in the fasta file %s. Exiting.'%(fasta))
                sys.exit(1)

        previous_kmer = this_kmer

    fh.close()
    
list_of_contigs = [] # stores the order of contigs, before joining/merging step
index_of_contigs = {} # key: contig, value: index in the order of contigs
index_counter = 1

logging.info('Iterating over the k-mers now, creating connections between contigs')
# now iterate over all k-mers
current_kmer = first_kmer
while current_kmer.next_kmer != 'END':
    this_set_of_hits = current_kmer.set_of_hits
    if not this_set_of_hits:
        logging.debug('Kmer %s has no hits'%current_kmer.name)
        current_kmer = current_kmer.next_kmer
        continue
    # we have to sort by the first element explicitly, if two connections share the same position it will try and sort by Contig-objects and crash
    hits = sorted(this_set_of_hits, key=lambda x:x[0])
    # [(59669, tig00000001XXXuniq_kmer_id_chr1A.soap.533.contigs.fasta), (272044, tig00000061XXXuniq_kmer_id_chr1A.soap.532.contigs.fasta)]
    # note: the contigs are not strings but of type Contig

    for h in hits:
        contig = h[1]
        if contig.name not in index_of_contigs:
            list_of_contigs.append(contig)
            index_of_contigs[contig.name] = index_counter
            index_counter += 1
        # make the link here!
        for other_h in hits:
            if h == other_h: continue
            this_pos, this_contig = h # (59669, tig00000001XXXuniq_kmer_id_chr1A.soap.533.contigs.fasta)
            other_pos, other_contig = other_h # (272044, tig00000061XXXuniq_kmer_id_chr1A.soap.532.contigs.fasta)
            # we now link from 59669 to 272044 in this_contig
            if this_pos in this_contig.connections:
                this_contig.connections[this_pos].append( (other_pos, other_contig) )
            else:
                this_contig.connections[this_pos] = [ (other_pos, other_contig) ]
            if other_contig in this_contig.connections_stats:
                this_contig.connections_stats[other_contig] += 1
            else:
                this_contig.connections_stats[other_contig] = 1

    current_kmer = current_kmer.next_kmer

logging.info('A last stroll over the contigs, with merging and writing out')
# c is of type Contig
finished_contigs = set() # we need to store which contigs we've merged
supercontig_counter = 1 # we need a new name for the merged contigs :)
for c in list_of_contigs:
    if c.name in finished_contigs:
        # we have merged this contig into another one
        logging.info('Contig %s has been written out or merged somewhere else, skipping'%c.name)
        continue
    revs_fors = c.revs_fors
    revs_fors_total = sum(revs_fors.values())
    if revs_fors_total < min_placement:
        logging.info('Contig %s has too few hits in total (%s), skipping'%(c.name, revs_fors_total))
        txt_out.write('%s\tSkipped, few hits (%s)\n'%(c.name, revs_fors_total))
        continue

    # is the contig chimeric?
    if (revs_fors['rev'] > 0) and (revs_fors['for'] > 0):
        logging.info('Contig %s has %s forward and %s reverse hits, could be chimeric'%(c.name, revs_fors['for'], revs_fors['rev']))
        c.chimeric = 'Maybe'
        txt_out.write('%s\tChimeric (%s forward and %s reverse hits)\n'%(c.name, revs_fors['for'], revs_fors['rev']))

    if revs_fors['rev'] > revs_fors['for']:
        logging.info('Reverse complementing contig %s'%(c.name))
        c.rev_complement()
        txt_out.write('%s\tReverse complemented\n'%c.name)
    c.checked_reverse = True

    # are there links for this contig?
    if not c.connections:
        # no links, just write out
        logging.info('Contig %s has no hits, writing out'%c.name)
        txt_out.write('%s\tNo merging\n'%(c.name))
        fasta_out.write(c.fasta())
        continue

    # there are links, we have to merge as much as possible
    this_connections = c.connections
    this_connections_stats = c.connections_stats

    # first, we have to find out which connections with contigs we have to ignore - example:
    # {tig00000005XXXuniq_kmer_id_chr1A.soap.532.contigs.fasta: 2, 
    # tig00000007XXXuniq_kmer_id_chr1A.soap.532.contigs.fasta: 5, 
    # tig00000001XXXuniq_kmer_id_chr1A.soap.533.contigs.fasta: 648}
    # in this case, only the last connection is meaningful

    contigs_to_ignore = set()
    for other_contig in this_connections_stats:
        contig_connection_count = this_connections_stats[other_contig]
        if contig_connection_count < overlap_cutoff:
            logging.info('Ignoring connection between contig %s and contig %s, not enough hits (%s)'%(c.name, other_contig.name, contig_connection_count))
            # we have to store this as a pair - there can still be meaningful connections with other contigs
            contigs_to_ignore.add( tuple(sorted([c.name, other_contig.name])) )

            # we also have to clean those connections between this contig and the other contig that we want deleted
            new_c_connections = {}
            for connection in c.connections:
                # connection : 348954
                for position in c.connections[connection]:
                    # position: (7156, Object:tig00000023XXXuniq_kmer_id_chr3D.soap.590.contigs.fasta)
                    position_pos, position_contig = position
                    if other_contig != position_contig:
                        if connection not in new_c_connections:
                            new_c_connections[connection] = [position]
                        else:
                            new_c_connections[connection].append(position)

            c.connections = dict(new_c_connections)
            # now clean the other too
            new_c_connections = {}
            for connection in other_contig.connections:
                for position in other_contig.connections[connection]:
                    position_pos, position_contig = position
                    if c != position_contig:
                        if connection not in new_c_connections:
                            new_c_connections[connection] = [position]
                        else:
                            new_c_connections[connection].append(position)
            other_contig.connections = dict(new_c_connections)

    # do we have connections left??
    if len(contigs_to_ignore) == len(this_connections_stats):
        logging.info('After cleaning, contig %s has no overlaps left, writing out'%(c.name))
        fasta_out.write(c.fasta())
        txt_out.write('%s\tNo merging after filtering of connections\n'%(c.name))
        continue

    # now do the actual merging
    logging.info('Contig %s has overlaps, merging'%(c.name))
    # we need to store that we merged this contig
    finished_contigs.add(c.name)

    # we need to store, for each connected contig, the max and the min position of the connection

    # first, traverse the connections of this contig to the next contig, 
    # to the next contig etc. until we have all connected contigs
    set_of_contigs_connected_here = set([c.name])
    for connection in this_connections:
        other_connections = this_connections[connection]

        # connection: 462
        # other_connection = [(134691, tig00000005XXXuniq_kmer_id_chr1A.soap.532.contigs.fasta)]
        for other_connection in other_connections:
            other_position, other_contig = other_connection
            if tuple(sorted([c.name, other_contig.name])) in contigs_to_ignore:
                continue
            set_of_contigs_connected_here.add(other_contig.name)
            # RECURSION MAGIC WOOOOO
            set_of_contigs_connected_here = get_contigs_recursively(other_contig, set_of_contigs_connected_here, contigs_to_ignore)

    # set_of_contigs_connected_here:
    # {'tig00000005XXXuniq_kmer_id_chr1A.soap.532.contigs.fasta', 
    # 'tig00000001XXXuniq_kmer_id_chr1A.soap.533.contigs.fasta', 
    # 'tig00000007XXXuniq_kmer_id_chr1A.soap.532.contigs.fasta',
    # 'tig00000061XXXuniq_kmer_id_chr1A.soap.532.contigs.fasta'}
    # now we have to check which of those contigs should be reverse complemented
    for contig_name in set_of_contigs_connected_here:
        contig = dict_of_contigs[contig_name]
        if not contig.checked_reverse:
            if contig.revs_fors['rev'] > contig.revs_fors['for']:
                logging.info('Reverse complementing contig %s'%(contig.name))
                contig.rev_complement()
                txt_out.write('%s\tReverse complemented\n'%contig.name)
            this_contig.checked_reverse = True

    # now iterate over those contigs
    min_max_dict = {} # key: contig, value: another dictionary, key: other_contig, value: min, max of connections

    for contig_name in set_of_contigs_connected_here:
        # c is of type String
        this_dict = {}
        
        contig = dict_of_contigs[contig_name]
        # now iterate over all hits of this contig, store the min and max for each partner - 
        # usually there is only one partner but we have to make sure...
        thisseen = set()
        for connection in contig.connections:
            # c, connection, contig.connections[connection]
            # tig00000007XXXuniq_kmer_id_chr1A.soap.532.contigs.fasta 1997 [(32958, Object:tig00000061XXXuniq_kmer_id_chr1A.soap.532.contigs.fasta)]
            # 1997 is the position of the connection on tig00000007XXXuniq_kmer_id_chr1A.soap.532.contigs.fasta
            # 32958 is the position of the connection on tig00000061XXXuniq_kmer_id_chr1A.soap.532.contigs.fasta
            if contig.connections[connection] and contig.connections[connection][0][1] not in thisseen:
                thisseen.add(contig.connections[connection][0][1])
            for other_connection in contig.connections[connection]:
                # other_connection is now (32958, Object:tig00000061XXXuniq_kmer_id_chr1A.soap.532.contigs.fasta)
                other_position, other_contig = other_connection

                # above we learned to ignore the connections we have to skip, same here!
                if tuple(sorted([contig.name, other_contig.name])) in contigs_to_ignore:
                    continue

                if other_contig.name not in set_of_contigs_connected_here:
                    continue

                if other_contig not in this_dict:
                    this_dict[other_contig] = [connection]
                else:
                    this_dict[other_contig].append(connection)

        new_dict = {}
        for k in this_dict:
            values = sorted(this_dict[k])
            if len(values) < overlap_cutoff:
                continue
            new_dict[k] = (min(values), max(values), len(values))
        this_dict = new_dict
        # do we even have connections after skipping everything?
        if this_dict:
            min_max_dict[contig_name] = this_dict


    logging.info('Merging %s contigs (%s)'%(len(min_max_dict), ','.join(min_max_dict.keys())))
    # we need the order of the contigs first, then we can just iterate over them and merge up
    order = {} # key: index of contig, value: contig name

    for contig in min_max_dict:
        this_index = index_of_contigs[contig]
        order[this_index] = contig

    '''
    order = {2: 'tig00000061XXXuniq_kmer_id_chr1A.soap.532.contigs.fasta', 4: 'tig00000001XXXuniq_kmer_id_chr1A.soap.533.contigs.fasta'}
    # yes, that's 2 and 4 - since we ignore some overlaps we lost contigs 1 and contigs 3

    min_max_dict = 
    {'tig00000061XXXuniq_kmer_id_chr1A.soap.532.contigs.fasta': 
            {Object:tig00000001XXXuniq_kmer_id_chr1A.soap.533.contigs.fasta: [212717, 280254]}, 
    'tig00000001XXXuniq_kmer_id_chr1A.soap.533.contigs.fasta': 
            {Object:tig00000061XXXuniq_kmer_id_chr1A.soap.532.contigs.fasta: [350, 67876]}}

    '''
    sorted_indexes = sorted(order.keys())

    # annoying thing - sometimes, my 'trick' of sortig the contigs like this does not work
    # for example, contig 1 may have a link with contig 3, but not with contig 2
    # so we have to make sure that all contigs in the order are actually connected, and if necessary, reshuffle.
    # ASSUMPTION: contig 1 is actually contig 1 and not contig 4

    new_order = {} # Tell me how do I feel 

    index_counter = 1
    this_contig = order[sorted_indexes[0]]

    # TODO
    # PERHAPS DELETE THE FOLLOWING
    # TODO
    # FIND OUT WHY THE LENGTH OF ORDER IS 4, SHOULD BE 2
    new_order[index_counter] = this_contig
    # now walk through the min/max for each partner
    seen = set([this_contig])
    while True:
        longest_length = 0
        longest_partner = ''
        for partner in min_max_dict[this_contig]:
            # partner is Object: Contig
            if partner.name in seen:
                # if we don't skip partners we've seen we'll 'look backwards' from the second contig to the first
                continue
            seen.add(partner.name)
            this_min, this_max, this_length = min_max_dict[this_contig][partner]
            length = abs(this_min - this_max)
            if length > longest_length:
                longest_length = length
                longest_partner = partner.name
        index_counter += 1
        if not longest_partner:
            # stop searching
            break
        index_counter += 1
        new_order[index_counter]= longest_partner
        finished_contigs.add(longest_partner)
        this_contig = longest_partner
            
    #for i in sorted(new_order):
    #print(i, new_order[i])
    # TODO: CHECK IF THIS WORKS
    # TODO: SOME CONTIGS WHO HAVE BEEN PREVIOUSLY IN order ARE NOT BEING WRITTEN OUT, CHECK
    order = new_order
    sorted_indexes = sorted(order.keys())

    # we have to store which contigs we actually ended up merging
    finished_contigs.add(c.name)

    # now merge
    supercontig_name = 'Supercontig_%s'%(supercontig_counter)
    supercontig_counter += 1
    supercontig_sequence = []

    all_names = []
    for position, index in enumerate(sorted_indexes):
        this_contig_name = order[index]
        all_names.append(this_contig_name) # need this for the stats file
        this_contig = dict_of_contigs[this_contig_name]
        this_contig_sequence = this_contig.sequence
        if position == 0:
            # |this_contig|
            #         |other_contig|
            next_contig_name = order[sorted_indexes[position+1]]
            next_contig = dict_of_contigs[next_contig_name]
            this_min, this_max, this_length = min_max_dict[this_contig_name][next_contig]
            logging.info('Pulling out until %s in contig %s (total length: %s)'%(this_min, this_contig_name, len(this_contig.sequence)))
            this_piece = this_contig_sequence[:this_min]
        elif position == len(sorted_indexes) - 1:
            # this is the last contig, cut towards the end
            # |other_contig|
            #       |this_contig|
            previous_contig_name = order[sorted_indexes[position-1]]
            previous_contig = dict_of_contigs[previous_contig_name]
            this_min, this_max, this_length = min_max_dict[this_contig_name][previous_contig]
            this_piece = this_contig_sequence[this_max:]
            logging.info('Pulling out from %s in contig %s (total length: %s)'%(this_max, this_contig_name, len(this_contig.sequence)))
        else:
            # |other_contig_a|    |other_contig_b|
            #            |this_contig|
            previous_contig_name = order[sorted_indexes[position-1]]
            previous_contig = dict_of_contigs[previous_contig_name]
            
            next_contig_name = order[sorted_indexes[position+1]]
            next_contig = dict_of_contigs[next_contig_name]

            previous_min, previous_max, previous_length = min_max_dict[this_contig_name][previous_contig]
            next_min, next_max, next_length = min_max_dict[this_contig_name][next_contig]
            # this_contig is cut from 'max' to the previous contig
            # to 'min' to the next contig
            this_piece = this_contig_sequence[previous_max:next_min]
            logging.info('Pulling out from %s to %s in contig %s'%(previous_max, next_min, this_contig_name))
            if not this_piece:
                logging.info('After pulling out, contig %s was reduced to nothing, i.e., it may be completely covered by other contigs'%(this_contig_name))

        supercontig_sequence.append(this_piece)

    supercontig_sequence = ''.join(supercontig_sequence)
    supercontig = Contig(supercontig_name, 'Merged', supercontig_sequence)
    logging.info('New supercontig: %s'%(supercontig.name))
    txt_out.write('%s\tSupercontig: connects %s\n'%(supercontig.name, ','.join(all_names)))
    fasta_out.write(supercontig.fasta())

logging.info('Done! Have a nice day :)')
sys.exit(0)
