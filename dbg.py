# IMPORT GRAPH BASIC STRUCTURES
from node import Node
from edge import Edge
from contig import Contig
# IMPORT BUBBLES AND TIPS FOR TOPOLOGICALLY BASED PRUNING
from bubble import Bubble
from tip import Tip
# IMPORT HELPER FUNCTIONS
from kmer import kmer_generate
from parser import parser
from utils import rc, median_qual

from statistics import median

import warnings

class Dbg:
    # COSTRUZIONE
    ### Dizionario nodi e archi
    ### Arrivano le sequenze parsate
    ### Genero k-mers
    ### Calcolo prefisso e suffisso
    ### Se k-mero nel dizionario degli archi aumento molteplicità
    ### Altrimenti, controllo se ho un oggetto tra nodo prefix e suffix
    ### Altrimenti creo oggetto nodo prefix e/o oggetto nodo suffix 
    ### Se creo nodi alla fine creo arco
    ### In generale aggiungo archi alla lista di archi in entrata e uscita dei nodi

    # PRUNING
    ### Tolgo edges a troppo bassa molteplicity

    # CAMMINO
    ### Parto casualmente, arrivo fino a un hub (outbound edges > 2)
    ### Se riesco risolvo, altrimenti mi fermo
    """

    """

    def __init__(self, k, file_name, file_type, paired = False, file_nameR2 = None, reverse_r2 = False):
        if not isinstance(k, int):
            raise TypeError('k parameter must be integer')
        if not isinstance(file_name, str):
            raise TypeError('file_name must be string')
        if not isinstance(file_type, str):
            raise TypeError('file_type must be string')
        if not isinstance(paired, bool):
            raise TypeError('paired must be boolean')
        if not isinstance(reverse_r2, bool):
            raise TypeError("reverse_r2 must be boolean")
        if file_type not in ['fa', 'fq', 'fasta', 'fastq']:
            raise IOError("file_type must be one of the following: ['fa', 'fq', 'fasta', 'fastq']")
        if not (3 <= k <= 101):
            raise ValueError(f"k must be between 3 and 101. Received: {k}")
        if k % 2 == 0:
            warnings.warn("k is even. Using an odd k is highly recommended to avoid topological complications")
        self.k = k #lunghezza kmeri
        self.file_nameR1 = file_name
        if file_type in ['fa', 'fasta']:
            self.file_type = 'fasta'
        elif file_type in ['fq', 'fastq']:
            self.file_type = 'fastq'
        if file_type in ['fa', 'fasta'] and paired:
            print("Data in fasta format cannot be paired")
            paired = False
            print("paired set to False")
        self.paired = paired
        self.reverse_r2 = False
        if self.paired == True:
            if not isinstance(file_nameR2, str):
                raise TypeError("file_nameR2 must be string")
            else:
                self.file_nameR2 = file_nameR2
                self.reverse_r2 = reverse_r2
        else:
            if file_nameR2 is not None:
                print("You inputed file_nameR2 but paired is false or file_type is fasta. It will not be used")
            if reverse_r2:
                print("You set reverse_r2 to True but paired is false. It will not be used")

        self.graph_built = False
        self.pruned_graph = False
        self.bubbles_processed = False
        self.tips_processed = False
        self.isolated_contigs_removed = False
        self.contigs_assembled = False
        self.contigs_forked_in_hub = False
        self.contigs_written = False
        self.contigs_output_fasta = ''
        self.nodes = {}  #dizionario nodi
        self.edges = {}  #dizionario archi
        self.contigs = []
        self.bubbles = {}
        self.tips = {}

    def build_graph(self, quality_check = False, quality_threshold = 30, use_both_strands = True):
        def put_reads_graph(reads, quality_check, use_both_strands, reverse_complement = False):
            reads_discarded = 0
            for read_tup in reads:
                read = read_tup[0]
                if reverse_complement:
                    read = rc(read)
                if quality_check:
                    read_qual = read_tup[1]
                    if median_qual(read_qual) > quality_threshold:
                        self.add_read(read, use_both_strands = use_both_strands)
                    else:
                        reads_discarded += 1
                else:
                    self.add_read(read, use_both_strands = use_both_strands)
            if quality_check:
                print(f"Discarded {reads_discarded} reads due to low median quality")
            return reads_discarded

        if not isinstance(quality_check, bool):
            raise TypeError("quality_check must be boolean")
        if not isinstance(quality_threshold, int):
            raise TypeError("quality_threshold must be integer")
        if not isinstance(use_both_strands, bool):
            raise TypeError("use_both_strands must be boolean")
        if self.file_type == 'fasta' and quality_check:
            print("You can't perform a quality check if your reads are in fasta format.")
            quality_check = False
            print("quality check disabled")
        if self.reverse_r2 and use_both_strands:
            print("You are using both strands after setting reverse_r2 to True. This wil override reverse_r2")

        reads1 = parser(file_name=self.file_nameR1, file_type=self.file_type, quality= quality_check) #generator
        put_reads_graph(reads1, quality_check, reverse_complement = False, use_both_strands = use_both_strands)
        if self.paired:
            reads2 = parser(file_name=self.file_nameR2, file_type=self.file_type, quality= quality_check) #generator
            put_reads_graph(reads2, quality_check, reverse_complement = self.reverse_r2, use_both_strands = use_both_strands)
        self.graph_built = True
        
    def add_read(self, read: str, use_both_strands: bool):
        """
            Decomposes a single DNA read into overlapping k-mers and updates the graph.

            For every k-mer in the read:
            1. It identifies the Prefix (Source Node) and Suffix (Destination Node).
            2. It creates or updates the Edge connecting them.
            3. It increases the Edge weight (coverage).

            :param read: The DNA sequence string.
            :type read: str
            :param use_both_strands: Whether to generate k-mers for the reverse complement too.
            :type use_both_strands: bool
        """
        if not isinstance(read, str):
            raise TypeError("read must be str")
        if not isinstance(use_both_strands, bool):
            raise TypeError("use_both_strands must be boolean")

        for kmer in kmer_generate(read, self.k, reverse= use_both_strands):
            prefix = kmer[:-1] #calcolo prefisso e suffisso
            suffix = kmer[1:]

            if prefix not in self.nodes: #se non ci sono nel dizionario creo i nodi e li aggiungo al dict
                self.nodes[prefix] = Node(prefix)
            if suffix not in self.nodes:
                self.nodes[suffix] = Node(suffix)

            node_prefix = self.nodes[prefix]
            node_suffix = self.nodes[suffix]

            if kmer in self.edges: #se l'arco è nel dizionario aumento la molteplicità
                self.edges[kmer].weight += 1
            else:
                new_edge = Edge(kmer, node_prefix, node_suffix) #altrimenti creo l'arco, lo aggiungo al dizionario e...
                self.edges[kmer] = new_edge
                node_prefix.outgoing.append(new_edge) #aggiorno gli archi in entrata e in uscita dai nodi
                node_suffix.incoming.append(new_edge)

    
    def prune_graph(self, threshold: int=2):
        """
            Removes edges with low coverage to clean up sequencing errors.

            Low-frequency k-mers are usually  caused by random sequencing errors.
            Removing them drastically simplifies the graph topology.

            :param threshold: Edges with weight < threshold will be deleted. Defaults to 2.
            :type threshold: int
            :raises RuntimeError: If called before the graph is built.
        """
        if not self.graph_built:
            raise RuntimeError("You cannot prune a graph if it does not exist")
        if not isinstance(threshold, int):
            raise TypeError("threshold must be integer")
        kmers_to_remove = []
        for kmer, edge in self.edges.items():
            if edge.weight < threshold:
                kmers_to_remove.append(kmer) #aggiungo i kmeri da rimuovere alla lista

        count = 0
        for kmer in kmers_to_remove:
            edge_to_delete = self.edges[kmer]
            
            if edge_to_delete in edge_to_delete.source.outgoing:
                edge_to_delete.source.outgoing.remove(edge_to_delete) #rimuovo dalla lista outgoing l'arco

            if edge_to_delete in edge_to_delete.destination.incoming:
                edge_to_delete.destination.incoming.remove(edge_to_delete) #rimuovo dalla lista incoming l'arco

        
            del self.edges[kmer] #rimuovo dal dict
            count += 1
        
        print(f"{count} edges removed")
        self._rm_isolated_nodes()
        self.pruned_graph = True

    def get_graph_stats(self):
        """
            Prints a summary of the current graph status.

            It reports:
            - Total nodes and edges.
            - Number of hubs (complex junctions).
            - Average coverage.
            - Status of cleaning steps (Bubbles/Tips).
        """
        if not self.graph_built:
            raise RuntimeError("You cannot get stats of a graph if it does not exist")
        num_nodes = len(self.nodes) #numero nodi
        num_edges = len(self.edges)
        total_weight = 0
        
        if num_edges > 0:
            total_weight = sum(edge.weight for edge in self.edges.values())
            
            avg_coverage = total_weight / num_edges
        else:
            avg_coverage = 0
        
        num_hubs = len(list(node for node in self.nodes.values() if node.is_hub()))
        num_rec_hubs = len(list(node for node in self.nodes.values() if node.is_receiving_hub()))
        print("Graph Statistics: ")
        print(f"Nodes: {num_nodes}")
        print(f"  Of which {num_hubs} hubs and {num_rec_hubs} receiving hubs")
        print(f"Edges: {num_edges}")
        print(f"Average coverage: {avg_coverage:.2f}X")
        if self.bubbles_processed:
            print(f"Bubbles remaining after processing: {len(self.bubbles)}")
        if self.tips_processed:
            print(f"Tips remaining after processing: {len(self.tips)}")

    def write_contigs(self, output_filename: str ="contigs.fasta", min_len: int =0, verbose: bool = False) -> list:
        """
            Saves the assembled contigs to a FASTA file.

            :param output_filename: Name of the output file. Defaults to "contigs.fasta".
            :type output_filename: str
            :param min_len: Minimum length (bp) for a contig to be saved.
            :type min_len: int
            :param verbose: If True, prints details about every saved contig.
            :type verbose: bool
            :return: A list of lengths of the saved contigs (useful for N50 calculation).
            :rtype: list[int]
        """
        if not isinstance(output_filename, str):
            raise TypeError("output_filename must be string")
        if not isinstance(min_len, int):
            raise TypeError("min_len must be integer")
        if not isinstance(verbose, bool):
            raise TypeError("verbose must be boolean")
        if not self.contigs_assembled or not self.contigs:
            warnings.warn("Warning! No contigs in memory")
            return []

        saved_lengths = []
        count = 0

        try:
            with open(output_filename, "w") as f:
                i = 0
                for contig in self.contigs:
                    if len(contig) < min_len:
                        continue
                    median_cov = median(contig.walk_weights) if contig.walk_weights else 0
                    f.write(f">contig_{i + 1} len={len(contig)} tag={contig.tag} median_cov={median_cov}\n")
                    if verbose:
                        print(f"Written contig_{i + 1} len={len(contig)} tag={contig.tag} median_cov={median_cov}")
                    f.write(f"{contig.seq}\n")
                    saved_lengths.append(len(contig))
                    count += 1
                    i += 1

            print(f"File saved: {output_filename} ({count} contigs written, min_len={min_len})")
            self.contigs_written = True
            self.contigs_output_fasta = output_filename
            return saved_lengths

        except IOError as e:
            print(f"Error writing file: {e}")
            return []

    
    @staticmethod
    def n50(lengths: list[int]) -> int:
        """
            Calculates the N50 metric for assembly quality.

            N50 is a standard bioinformatic metric. It is the length of the shortest contig
            such that the sum of contigs of this length or longer is at least 50%
            of the total assembly length.

            :param lengths: A list of contig lengths.
            :type lengths: list[int]
            :return: The N50 value.
            :rtype: int
        """
        #lista di lunghezze
        if not isinstance(lengths, list):
            raise TypeError("lengths must be list")
        if not lengths:
            return 0
            
        sorted_lengths = sorted(lengths, reverse=True) #ordinamento decres
        total_len = sum(sorted_lengths) #somma lengths
        threshold = total_len / 2 #soglia
        
        current_sum = 0
        for length in sorted_lengths:
            current_sum += length
            if current_sum >= threshold:
                return length #n50


    def _rm_isolated_nodes(self):
        """
            Helper method to remove nodes that have no connections.
        """
        nodes_to_rm = []
        for kmer, node in self.nodes.items():
            if len(node.incoming) == 0 and len(node.outgoing) == 0:
                nodes_to_rm.append(kmer)
        
        for kmer in nodes_to_rm:
            del self.nodes[kmer]
        
        print(f"removed {len(nodes_to_rm)} isolated nodes")


    def _rm_path_edges(self, path: list["Edge"]):
        """
             Helper method to physically delete a list of edges (a path) from the graph.
             It disconnects the edges from their source/destination nodes and removes
             them from the main edge dictionary.

            :param path: The sequence of edges to remove
            :type path: list[Edge]
        """

        for edge in path:
            if edge in edge.source.outgoing:
                edge.source.outgoing.remove(edge)
            if edge in edge.destination.incoming:
                edge.destination.incoming.remove(edge)
            if edge.seq in self.edges:
                del self.edges[edge.seq]

    def _detect_bubbles(self, max_len: int):
        """
            Generator that searches for bubble structures in the graph.

            It looks for "forks" (nodes with multiple outgoing edges) and tries to find
            paths that reconverge at a specific "sink" node within a certain distance (`max_len`).

            :param max_len: The maximum length (in edges) of a bubble path.
            :type max_len: int
            :yield: Bubble objects found in the graph.
        """
        #generatore

        forks = [n for n in self.nodes.values() if n.is_hub()] #lista nodi hub

        for fork_node in forks:
            candidates = {}  #chiave: sequenza nodo arrivo, valore: lista path(lista di di edges)


            for edge in fork_node.outgoing:
                path = [edge]
                curr = edge.destination
                valid = True #diventa False quando il nodo è un hub o tip

                while len(path) < max_len: #se è lungo piu di 2k+1, lo tengo
                    if curr.is_receiving_hub():
                        break #chiusura bolla

                    if len(curr.outgoing) != 1:
                        valid = False #tanto itero per la gestione delle bolle annidate
                        break

                    next_edge = curr.outgoing[0]
                    path.append(next_edge) #aggiungo edge al path
                    curr = next_edge.destination

                #se il percorso è valido e chiuso
                if valid and curr.is_receiving_hub():
                    dest_seq = curr.seq
                    if dest_seq not in candidates: 
                        candidates[dest_seq] = []
                    candidates[dest_seq].append(path) #aggiorno dict


            for dest_seq, paths in candidates.items():
                if len(paths) > 1: #cioè, se quel nodo è davvero un hub, ha più path in entrata:
                    dest_node = self.nodes[dest_seq]
                    yield Bubble(fork_node, dest_node, paths) #restituisco l'oggetto bolla

    def process_bubbles(self, cov_threshold: int=5) -> int:
        """
            Iteratively finds and resolves bubbles in the graph.
            This method keeps calling  `_detect_bubbles` and solving them until no more
            bubbles can be resolved. It handles nested bubbles by running multiple iterations.

            :param cov_threshold: Coverage threshold used to distinguish errors from real variants.
            :type cov_threshold: int
            :return: Total number of bubbles resolved (popped).
            :rtype: int
        """
        if not self.graph_built:
            raise RuntimeError("You cannot process bubbles in a graph if it does not exist")
        if not isinstance(cov_threshold, int):
            raise TypeError("cov_threshold must be integer")

        max_len = (self.k * 2) + 1
        iteration = 1
        removed_bubbles = 0
        while True: #itero, vengono risolte prima le bolle interne
            popped_count = 0
            found_bubbles = False
            try:
                for bubble in self._detect_bubbles(max_len):
                    found_bubbles = True
                    if bubble._resolve(self, cov_threshold):
                        popped_count += 1
            except ValueError:
                print("There have been problem in tips detection. You can't perform the operation")
                break
            if not found_bubbles: #sono tutte risolte
                break

            print(f"Iteration {iteration}: {popped_count} bubbles popped")
            if popped_count == 0: #magari ci sono ancora bolle, ma non vengono risolte
                for bubble in self._detect_bubbles(max_len):
                    self.bubbles[bubble.source.seq] = bubble
                break
            removed_bubbles += popped_count
            iteration += 1

        self._rm_isolated_nodes()
        self.bubbles_processed = True
        return removed_bubbles

    def _detect_tips(self, max_len: int):
        """
            Generator that searches for dead ends (Tips) in the graph.

            It scans for:
            1. Type I Tips: Paths starting at a hub and dying out.
            2. Type II Tips: Paths starting from nowhere and merging into a hub.

            :param max_len: Maximum length of a tip to be considered. Longer dead ends are kept as valid contigs.
            :type max_len: int
            :yield: Tip objects found in the graph.
        """
        # Detection tips di tipo I
        potential_inital_nodes_I = [node for node in self.nodes.values() if node.is_hub()]
        for fork_node in potential_inital_nodes_I:
            for edge in fork_node.outgoing:
                path = [edge]
                curr = edge.destination
                valid = True #diventa False quando il nodo è un hub o tip

                while len(path) < max_len: #se è lungo piu di 2k+1, lo tengo
                    if len(curr.outgoing) == 0:
                        break #chiusura tip

                    if curr.is_hub():
                        valid = False # Tip annidati
                        break

                    next_edge = curr.outgoing[0]
                    path.append(next_edge) #aggiungo edge al path
                    curr = next_edge.destination

                #se il percorso è valido e viene troncato prima di max_len
                if valid and len(curr.outgoing) == 0:
                    yield Tip(fork_node, curr, path, 1) #restituisco l'oggetto Tip

       
        # Detection tips di tipo II
        potential_initial_nodes_II = [node for node in self.nodes.values() if len(node.incoming) == 0 and len(node.outgoing) == 1]
        for initial_node in potential_initial_nodes_II:
            #if initial_node.is_hub():
            #    continue # sarebbe trattato come un potenziale tips di tipo I
            # ^ Inutile se metto outgoing come 1
            edge = initial_node.outgoing[0]
            path = [edge]
            curr = edge.destination
            valid = True

            while len(path) < max_len:
                if curr.is_receiving_hub(): # tip rilevata
                    break
                if len(curr.outgoing) == 0:
                    valid = False # Praticamente un mini contig inferiore a max_len
                    break
                if curr.is_hub():
                    valid = False
                    break
                next_edge = curr.outgoing[0]
                path.append(next_edge)
                curr = next_edge.destination
            
            if valid and curr.is_receiving_hub():
                yield Tip(initial_node, curr, path, 2)


    def process_tips(self, cov_threshold: int=5) -> int:
        """
            Iteratively removes dead ends (Tips) from the graph.

            This effectively "shaves" the hairy parts of the graph caused by
            sequencing errors at the ends of reads.

            :param cov_threshold: Coverage threshold to confirm a tip is an error.
            :type cov_threshold: int
            :return: Total number of tips removed.
            :rtype: int
        """
        if not self.graph_built:
            raise RuntimeError("You cannot process tips in a graph if it does not exist")
        if not isinstance(cov_threshold, int):
            raise TypeError("cov_threshold must be integer")
        max_len = (self.k * 2) + 1
        iteration = 1
        removed_tips = 0
        while True: #itero, vengono risolte prima le bolle interne
            removed_count = 0
            found_tips = False

            for tip in self._detect_tips(max_len):
                found_tips = True
                if tip._resolve(self, cov_threshold):
                    removed_count += 1

            if not found_tips:
                break

            print(f"Iteration {iteration}: {removed_count} tips popped")
            if removed_count == 0: #magari ci sono ancora tips ma non vanno eliminate
                for tip in self._detect_tips(max_len):
                    self.tips[tip.initial_node.seq] = tip
                break
            removed_tips += removed_count
            iteration += 1

        self._rm_isolated_nodes()
        self.tips_processed = True
        return removed_tips

    def process_short_islands(self):
        """
            Iteratively removes dead ends (Tips) from the graph.

            This effectively "shaves" the hairy parts of the graph caused by
            sequencing errors at the ends of reads.

            :param cov_threshold: Coverage threshold to confirm a tip is an error.
            :type cov_threshold: int
            :return: Total number of tips removed.
            :rtype: int
        """
        if not self.graph_built:
            raise RuntimeError("You cannot process islands in a graph if it does not exist")
        min_len = (self.k * 2) + 1

        visited_nodes = set() #set nodi visti
        nodes_to_delete = set() #set nodi da cancellare
        deleted_islands = 0

        #itero su tutti i nodi del grafo, potrebbero esserci isole circolari con len(self.incoming)!=0
        for start_node_seq in list(self.nodes.keys()):
            if start_node_seq in visited_nodes:
                continue

            if start_node_seq in nodes_to_delete:
                continue

            island_nodes = set()
            to_visit = [self.nodes[start_node_seq]] #lista con oggetto nodo
            island_nodes.add(start_node_seq)
            visited_nodes.add(start_node_seq)  #nodo visitato

            while to_visit:
                curr_node = to_visit.pop(0) #estragggo nodo

                #guardo vicini sia in entrata che in uscita
                all_neighbors = []
                for edge in curr_node.outgoing:
                    all_neighbors.append(edge.destination)
                for edge in curr_node.incoming:
                    all_neighbors.append(edge.source)

                for neighbor in all_neighbors:
                    if neighbor.seq not in island_nodes:
                        island_nodes.add(neighbor.seq)
                        visited_nodes.add(neighbor.seq)
                        to_visit.append(neighbor)


            #Formula: L = k + (N - 1)
            num_nodes = len(island_nodes)
            seq_len_bp = self.k + (num_nodes - 1)


            if seq_len_bp < min_len:
                nodes_to_delete.update(island_nodes)
                deleted_islands += 1

        #cancellazione
        if nodes_to_delete:
            #rimozione archi
            edges_to_remove = [k for k, e in self.edges.items()
                               if e.source.seq in nodes_to_delete or e.destination.seq in nodes_to_delete]
            for k in edges_to_remove:
                del self.edges[k]

            #rimozione nodi
            for k in nodes_to_delete:
                if k in self.nodes:
                    del self.nodes[k]

            print(
                f"Island Removal: Deleted {deleted_islands} short islands (< {min_len}bp). Removed {len(nodes_to_delete)} nodes.")
        else:
            print("Island Removal: No short islands found.")
        self.isolated_contigs_removed = True
        return deleted_islands


    def _walk_linear(self, start_edge):
        start_edge.visited = True
        seq = start_edge.seq
        curr = start_edge.destination
        walk_weights = [start_edge.weight]
        receiving_hubs_visited = set()
        if curr.is_receiving_hub():
            receiving_hubs_visited.add(curr.seq)
        # Continua solo se c'è un'unica uscita (niente hub)
        while len(curr.outgoing) == 1:
            edge = curr.outgoing[0]
            if edge.visited:
                break
            edge.visited = True
            seq += edge.seq[-1]
            walk_weights.append(edge.weight)
            curr = edge.destination
            if curr.is_receiving_hub():
                receiving_hubs_visited.add(curr.seq)
        return seq, curr, walk_weights, receiving_hubs_visited


    def assemble_contigs(self, fork_in_hub = False):
        if not self.graph_built:
            raise RuntimeError("You cannot assemble contigs in a graph if it does not exist")
        if not self.bubbles_processed:
            warnings.warn("You are assembling contigs before processing bubbles. This may cause topological complications and inaccuracies")
            if fork_in_hub:
                warnings.warn("Using fork_in_hub without processing bubbles first is strongly discouraged")
        if not self.tips_processed:
            warnings.warn("You are assembling contigs before processing tips. This may cause topological complications and inaccuracies")
        if not self.isolated_contigs_removed:
            warnings.warn("You are assembling contigs before removing isolated contigs. This may create small unrealistic contigs")
        
        all_starts = [node for node in self.nodes.values() if len(node.incoming) == 0 and len(node.outgoing) > 0]
        if not fork_in_hub:
            internal_hubs = [node for node in self.nodes.values() if len(node.incoming) > 0 and node.is_hub()] #considera anche gli hub come punti di inzio
            all_starts += internal_hubs
        contigs = []
        contigs_forked_in_hub = False
        for contig_start in all_starts:
            for start_edge in contig_start.outgoing:
                if start_edge.visited:
                    continue
                contig = Contig.from_walk_linear(self._walk_linear(start_edge), self.k, contig_start)
                forked = False
                if contig.final_node.is_hub() and fork_in_hub and not (contig.is_lariat() or contig.is_circular()):
                    contigs_forked_in_hub = True
                    forked = True
                    common_contig = contig.copy()
                    if contig.final_node.seq in self.bubbles and not self.bubbles[contig.final_node.seq].visited: #è l'inizio di una bolla non già visitata
                        bubble_to_assemble = self.bubbles[contig.final_node.seq] #prendo la bolla
                        for path in bubble_to_assemble.paths: #per ciasun path della bolla
                            for edge in path: #per ciascun edge nel path
                                edge.visited = True
                        alternative_paths = list(bubble_to_assemble._get_contigs_from_path(self.k)) #lista di oggetti contig con sequenza e lista pesi
                        new_start_node = bubble_to_assemble.dest_node
                        bubble_to_assemble.visited = True
                        if len(new_start_node.outgoing) == 0 or all([edge.visited for edge in new_start_node.outgoing]): #cioè, se non continua dopo la chiusura della bolla
                            for alt_contig in alternative_paths:
                                whole_contig = common_contig + alt_contig 
                                if whole_contig.is_circular():
                                    whole_contig.tag = "circular_bubble_resolved_ending" 
                                elif whole_contig.is_lariat():
                                    whole_contig.tag = "lariat_bubble_resolved_ending"
                                else:
                                    whole_contig.tag = "linear_bubble_resolved_ending"
                                contigs.append(whole_contig)
                        else:
                            for new_start_edge in new_start_node.outgoing:
                                if new_start_edge.visited:
                                    continue
                                new_contig = Contig.from_walk_linear(self._walk_linear(new_start_edge), self.k, new_start_node) #continua a camminare dopo la chiusura della bolla
                                for alt_contig in alternative_paths:
                                    whole_contig = common_contig + alt_contig + new_contig
                                    if whole_contig.is_lariat(): #########################
                                        whole_contig.tag = "lariat_with_bubble"
                                    elif whole_contig.is_circular():
                                        whole_contig.tag = "circular_with_bubble"
                                    else:
                                        whole_contig.tag = "linear_bubble_resolved_extended"
                                    contigs.append(whole_contig)
                        for remaining_edge in [edge for edge in contig.final_node.outgoing if not edge.visited]: # Edges che partono da next_node ma non sono della bolla
                            remaining_contig = Contig.from_walk_linear(self._walk_linear(remaining_edge), self.k, contig.final_node)
                            whole_contig = common_contig + remaining_contig
                            if whole_contig.is_lariat():
                                whole_contig.tag = "lariat_through_hub"
                            elif whole_contig.is_circular():
                                whole_contig.tag = "circular_through_hub"
                            else:
                                whole_contig.tag = "linear_through_hub"
                            contigs.append(whole_contig)
                    elif contig.final_node.seq in self.tips and not self.tips[contig.final_node.seq].visited:
                        tip_to_assemble = self.tips[contig.final_node.seq]
                        for edge in tip_to_assemble.path:
                            edge.visited = True
                        tip_contig = tip_to_assemble._get_contig_from_path(self.k)
                        tip_to_assemble.visited = True
                        whole_contig = common_contig + tip_contig
                        whole_contig.tag = "linear_with_tip"
                        contigs.append(whole_contig)
                    else: #non è una bolla
                        new_start_node = contig.final_node
                        for new_start_edge in new_start_node.outgoing:
                            if new_start_edge.visited:
                                continue
                            new_contig = Contig.from_walk_linear(self._walk_linear(new_start_edge), self.k, new_start_node)
                            whole_contig = common_contig + new_contig
                            if whole_contig.is_lariat():
                                whole_contig.tag = "lariat_through_hub"
                            elif whole_contig.is_circular():
                                whole_contig.tag = "circular_through_hub"
                            else:
                                whole_contig.tag = "linear_through_hub"
                            contigs.append(whole_contig)
                if not forked: #non ho trovato nessun bivio, oppure fork_in_hub spento
                    if contig.is_circular():
                        if contig.final_node.is_hub():
                            contig.tag = "circular_ending_hub"
                        else:
                            contig.tag = "circular_simple" #potrebbe essere un circolare che parte da un hub interno
                    elif contig.is_lariat():
                        if contig.final_node.is_hub():
                            contig.tag = "lariat_hub_branching"
                        else:
                            contig.tag = "lariat_simple"
                    else:
                        contig.tag = "linear_simple" 
                    contigs.append(contig)

        if fork_in_hub:
            # hubs restanti poiché fork_in_hub si ferma dopo un livello
            remaining_hubs_seqs = {node.seq for node in self.nodes.values() if len(node.incoming) > 0 and node.is_hub() and any(not edge.visited for edge in node.outgoing)}
            while remaining_hubs_seqs:
                for contig_start_seq in list(remaining_hubs_seqs):
                    contig_start = self.nodes[contig_start_seq]
                    for start_edge in contig_start.outgoing:
                        if start_edge.visited:
                            continue
                        remaining_contig = Contig.from_walk_linear(self._walk_linear(start_edge), self.k, contig_start)
                        if remaining_contig.is_circular():
                            remaining_contig.tag = "circular_fork_hub_residual"
                        elif remaining_contig.is_lariat():
                            remaining_contig.tag = "lariat_fork_hub_residual"
                        else:
                            remaining_contig.tag = "linear_fork_hub_residual"
                        contigs.append(remaining_contig)
                    if all(edge.visited for edge in contig_start.outgoing):
                        remaining_hubs_seqs.remove(contig_start_seq)

        not_visited_edges = [edge for edge in self.edges.values() if not edge.visited]
        if len(not_visited_edges) == 0:
            print("All edges have been visited. Eulerian path(s) found")
        else: #se non sono stati visitati tutti gli edge, vuol dire che c'è un path circolare
            print(f"There are still {len(not_visited_edges)} edges not visited. Starting circular contigs search.")
            if fork_in_hub:
                print("Fork in hub is not available for circular search")
            circular_contigs = 0
            unknown_contigs = 0
            while not_visited_edges:
                start_edge = not_visited_edges[0]
                start_node = start_edge.source
                circ_contig = Contig.from_walk_linear(self._walk_linear(start_edge), self.k, start_node)
                if circ_contig.is_circular():
                    circ_contig.tag = "circular_isolated"
                elif circ_contig.is_lariat():
                    circ_contig.tag = "unknown_contig_type_lariat_like"
                else:
                    circ_contig.tag = "unknown_contig_type"
                contigs.append(circ_contig)
                if circ_contig.tag == "circular_isolated":
                    circular_contigs += 1
                else:
                    print(f"Contig found with unknown topology {len(circ_contig)}")
                    unknown_contigs += 1
                not_visited_edges = [edge for edge in self.edges.values() if not edge.visited]
            print(f"All edges have been now visited. Saved {circular_contigs} circular contigs")
            print(f"Saved {unknown_contigs} contigs of unknown topology")
        # Reset visited per poter riusare l'assembly
        for edge in self.edges.values():
            edge.visited = False
        for bubble in self.bubbles.values():
            bubble.visited = False
        for tip in self.tips.values():
            tip.visited = False
        # Deduplicazione contigs
        deduplicated_contigs = []
        checked = set()
        for contig in contigs:
            if contig.seq not in checked:
                deduplicated_contigs.append(contig)
            checked.add(contig.seq)
            checked.add(rc(contig.seq))
        self.contigs_forked_in_hub = contigs_forked_in_hub
        self.contigs = deduplicated_contigs
        self.contigs_assembled = True
        if contigs_forked_in_hub:
            print("Some contigs have been forked in hub. N50 may not be totally accurate")
        print(f"Assembled a total of {len(deduplicated_contigs)} contigs")
        return deduplicated_contigs
