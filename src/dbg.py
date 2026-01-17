# IMPORT GRAPH BASIC STRUCTURES
from src.node import Node
from src.edge import Edge
from src.contig import Contig
# IMPORT BUBBLES AND TIPS FOR TOPOLOGICALLY BASED PRUNING
from src.bubble import Bubble
from src.tip import Tip
# IMPORT HELPER FUNCTIONS
from src.kmer import kmer_generate
from src.parser import parser
from src.utils import rc, median_qual

from collections.abc import Iterator, Generator
from statistics import median

import warnings

class Dbg:
    """
    De Bruijn Graph (DBG) Assembly Engine.
    
    This class manages the full genome assembly process: from read parsing 
    and graph construction to topological pruning (bubbles/tips/short islands) and 
    final contig extraction and writing to output files .

    :ivar k: K-mer length used for graph construction.
    :vartype k: int
    :ivar file_nameR1: Path to the primary sequence file.
    :vartype file_nameR1: str
    :ivar file_type: File format (must be either 'fasta' or 'fastq').
    :vartype file_type: str
    :ivar paired: Whether the assembly must proceed on paired-end reads. Only fastq can be paired
    :vartype paired: bool
    :ivar file_nameR2: Path to the secondary sequence file if paired.
    :vartype file_nameR2: str or None
    :ivar reverse_r2: Whether R2 reads are already reverted.
    :vartype reverse_r2: bool
    :ivar graph_built: Flag indicating if the graph has been constructed.
    :vartype graph_built: bool
    :ivar pruned_graph: Flag indicating if low-coverage edges were removed.
    :vartype pruned_graph: bool
    :ivar bubbles_processed: Flag indicating if bubbles were resolved.
    :vartype bubbles_processed: bool
    :ivar tips_processed: Flag indicating if dead-ends were removed.
    :vartype tips_processed: bool
    :ivar isolated_contigs_removed: Flag indicating if short isolated islands were cleaned.
    :vartype isolated_contigs_removed: bool
    :ivar contigs_assembled: Flag indicating if the path traversal has been completed.
    :vartype contigs_assembled: bool
    :ivar contigs_forked_in_hub: Flag indicating if deeper exploration beyond hubs was enabled.
    :vartype contigs_forked_in_hub: bool
    :ivar contigs_written: Flag indicating if the results were exported to a FASTA file.
    :vartype contigs_written: bool
    :ivar contigs_output_fasta: Name/path of the generated output file.
    :vartype contigs_output_fasta: str
    :ivar nodes: Storage for (k-1)-mer objects {sequence: Node}.
    :vartype nodes: dict
    :ivar edges: Storage for k-mer objects {sequence: Edge}.
    :vartype edges: dict
    :ivar contigs: List of assembled Contig objects.
    :vartype contigs: list
    :ivar bubbles: Collection of detected Bubble structures.
    :vartype bubbles: dict
    :ivar tips: Collection of detected Tip structures.
    :vartype tips: dict

    """

    def __init__(self, k: int, file_name: str, file_type: str, paired: bool = False, 
                 file_nameR2: str = None, reverse_r2: bool = False, offset = 33):
        """
        Initializes the DBG with assembly parameters and input file metadata.
        DOES NOT build the graph

        :param k: K-mer length (must be odd for better graph structure).
        :type k: int
        :param file_name: Path to the R1 sequence file.
        :type file_name: str
        :param file_type: Format that must be one of ('fasta', 'fastq', 'fa', 'fq').
        :type file_type: str
        :param offset: Represents the fastq offset for ASCII quality encoding. Can be 33 (default) or 64.
        :type offset: int
        :param paired: Boolean for paired-end data. If true, a second file needs to be provided and both will be used 
                       to build graph. Default is false.
        :type paired: bool
        :param file_nameR2: Path to the R2 sequence file (optional).
        :type file_nameR2: str or None
        :param reverse_r2: Whether to reverse complement R2 reads. Useful if you don't want to use both strands and you
                           know that your sequencer already reverse complemented the reads
        :type reverse_r2: bool
        """
        # Valido parametri
        if not isinstance(k, int):
            raise TypeError('k parameter must be integer')
        if not isinstance(file_name, str):
            raise TypeError('file_name must be string')
        if not isinstance(file_type, str):
            raise TypeError('file_type must be string')
        if not isinstance(offset, int):
            raise TypeError("offset must be integer")
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
        if offset not in [33,64]:
            raise ValueError("offset is not a standard offset value")
        # Setup grafo per costruzione

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
        self.offset = offset
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

        # Inizializzo stato del grafo
        self.graph_built = False
        self.pruned_graph = False
        self.bubbles_processed = False
        self.tips_processed = False
        self.isolated_contigs_removed = False
        self.contigs_assembled = False
        self.contigs_forked_in_hub = False
        self.contigs_written = False
        self.contigs_output_fasta = ''

        # Inizializzo strutture dati del grafo
        self.nodes = {} #dizionario nodi
        self.edges = {} #dizionario archi
        self.contigs = []
        self.bubbles = {}
        self.tips = {}

    def build_graph(self, quality_check: bool = False, quality_threshold: int = 30, use_both_strands: bool = True, verbose: bool = True):
        """
        Drives the construction of the De Bruijn Graph from input FASTA/FASTQ files.

        This method coordinates the parsing of sequencing reads (R1 and optionally R2),
        applies quality filters (brutal whole read median trimming) if requested (quality_check),
        and populates the graph by adding k-mers as edges and k-1 mers as nodes. 
        It handles both single-end and paired-end data.
        Can use k-mers from both strands of the reads to generate edges.

        :param quality_check: Enable filtering reads based on their median Phred score, defaults to False.
        :type quality_check: bool
        :param quality_threshold: Minimum median Phred score required to keep a read, defaults to 30.
        :type quality_threshold: int
        :param use_both_strands: If True, adds both the forward and the reverse complement of each k-mer, defaults to True.
        :type use_both_strands: bool
        :param verbose: If True, prints status messages and final graph statistics, defaults to True.
        :type verbose: bool
        
        :raises TypeError: If quality_check, quality_threshold, or use_both_strands are not of the expected type.
        :return: None
        """
        def put_reads_graph(reads: Iterator[tuple[str, str | None]], quality_check: bool, use_both_strands: bool, reverse_complement: bool = False, verbose: bool =True) -> int:
            """
            Internal helper to iterate over reads and add them to the graph.

            :param reads: An iterable (generator) of read tuples (sequence, quality|None) from the parser.
            :type reads: iterator
            :param quality_check: Whether to perform quality filtering.
            :type quality_check: bool
            :param use_both_strands: Whether to add reverse complements to the graph.
            :type use_both_strands: bool
            :param reverse_complement: If True, reverse complements the entire read before processing, defaults to False.
                                       Totally useless if use_both_strands is True
            :type reverse_complement: bool
            :param verbose: If True, prints the number of discarded reads, defaults to True.
            :type verbose: bool

            :return: The number of reads discarded due to low quality or N presence.
            :rtype: int
            """
            reads_discarded = 0
            reads_with_N = 0
            for read_tup in reads: #reads, generatore di tuple [0]seq, [1]qual
                read = read_tup[0]
                if 'N' in read:
                    reads_with_N += 1
                    continue
                if reverse_complement:
                    read = rc(read)
                if quality_check:
                    read_qual = read_tup[1]
                    if median_qual(read_qual, self.offset) > quality_threshold: # La mediana della quality dell'intera riga deve essere sopra soglia
                        self.add_read(read, use_both_strands = use_both_strands)
                    else:
                        reads_discarded += 1
                else:
                    self.add_read(read, use_both_strands = use_both_strands)
            if quality_check and verbose:
                print(f"Discarded {reads_discarded} reads due to low median quality and {reads_with_N} due to N presence")
            elif verbose and not quality_check:
                print(f"Discarded {reads_with_N} reads due to N presence")
            return reads_discarded+reads_with_N

        if not isinstance(quality_check, bool):
            raise TypeError("quality_check must be boolean")
        if not isinstance(quality_threshold, int):
            raise TypeError("quality_threshold must be integer")
        if not isinstance(use_both_strands, bool):
            raise TypeError("use_both_strands must be boolean")
        if not isinstance(verbose, bool):
            raise TypeError("verbose must be boolean")
        if self.file_type == 'fasta' and quality_check:
            print("You can't perform a quality check if your reads are in fasta format.")
            quality_check = False
            print("quality check disabled")
        if self.reverse_r2 and use_both_strands:
            print("You are using both strands after setting reverse_r2 to True. This wil override reverse_r2")

        reads1 = parser(file_name=self.file_nameR1, file_type=self.file_type, quality= quality_check, verbose = verbose) #generator
        put_reads_graph(reads1, quality_check, reverse_complement = False, use_both_strands = use_both_strands, verbose = verbose)
        if self.paired:
            reads2 = parser(file_name=self.file_nameR2, file_type=self.file_type, quality= quality_check, verbose = verbose) #generator
            put_reads_graph(reads2, quality_check, reverse_complement = self.reverse_r2, use_both_strands = use_both_strands, verbose = verbose)
        self.graph_built = True
        if not self.nodes:
            warnings.warn("No valid k-mers generated from reads. Check k value, read quality, or input files.")
        elif verbose:
            print(f"Graph built: {len(self.nodes)} nodes, {len(self.edges)} edges")
        
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

    
    def prune_graph(self, threshold: int=2, verbose: bool=True):
        """
            Removes edges with low coverage to clean up sequencing errors.

            Low-frequency k-mers are usually  caused by random sequencing errors.
            Removing them drastically simplifies the graph topology.

            :param threshold: Edges with weight < threshold will be deleted. Defaults to 2.
            :type threshold: int
            :param verbose: If True, it will print details about the pruned edges and nodes
            :type verbose: bool
            :raises RuntimeError: If called before the graph is built.
        """
        if not self.graph_built:
            raise RuntimeError("You cannot prune a graph if it does not exist")
        if not isinstance(threshold, int):
            raise TypeError("threshold must be integer")
        if not isinstance(verbose, bool):
            raise TypeError("verbose must be boolean")
        kmers_to_remove = []
        for kmer, edge in self.edges.items():
            if edge.weight < threshold:
                kmers_to_remove.append(kmer) #aggiungo i kmeri da rimuovere alla lista

        count = 0
        for kmer in kmers_to_remove:
            edge_to_delete = self.edges[kmer] #oggetto edge da eliminare
            
            if edge_to_delete in edge_to_delete.source.outgoing:
                edge_to_delete.source.outgoing.remove(edge_to_delete) #rimuovo dalla lista outgoing del nodo l'arco

            if edge_to_delete in edge_to_delete.destination.incoming:
                edge_to_delete.destination.incoming.remove(edge_to_delete) #rimuovo dalla lista incoming del nodo l'arco

            del self.edges[kmer] #rimuovo dal dict
            count += 1
        
        if verbose:
            print(f"{count} edges removed")
        self._rm_isolated_nodes(verbose)
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
                for i, contig in enumerate(self.contigs):
                    if len(contig) < min_len:
                        continue
                    
                    # Chiamata al metodo che genera l'header complesso
                    header = contig.get_fasta_header(i + 1)
                    
                    f.write(f"{header}\n")
                    f.write(f"{contig.seq}\n")
                    
                    saved_lengths.append(len(contig))
                    count += 1
                    
                    if verbose:
                        print(f"Written {header}")
            if verbose:
                print(f"File saved: {output_filename} ({count} contigs, min_len={min_len})")
            self.contigs_written = True
            self.contigs_output_fasta = output_filename
            return saved_lengths

        except IOError as e:
            print(f"Error writing file: {e}")
            return []

    
    @staticmethod
    def n50(lengths: list[int]) -> int:
        """
            Static method that calculates the N50 metric for assembly quality.

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


    def _rm_isolated_nodes(self, verbose: bool=True):
        """
            Helper method to remove nodes that have no connections.
        """
        if not isinstance(verbose, bool):
            raise TypeError("verbose must be boolean")
        nodes_to_rm = []
        for kmer, node in self.nodes.items():
            if len(node.incoming) == 0 and len(node.outgoing) == 0:
                nodes_to_rm.append(kmer)
        
        for kmer in nodes_to_rm:
            del self.nodes[kmer]
        if verbose:
            print(f"Removed {len(nodes_to_rm)} isolated nodes")


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
                edge.source.outgoing.remove(edge) #rimozione dalla lista di archi nel nodo da cui parte
            if edge in edge.destination.incoming:
                edge.destination.incoming.remove(edge) #rimozione dalla lista di archi nel nodo a cui arriva
            if edge.seq in self.edges:
                del self.edges[edge.seq]

    def _detect_bubbles(self, max_len: int) -> Iterator["Bubble"]:
        """
            Generator that searches for bubble structures in the graph, if shorter than a max_len parameter.

            It looks for "forks" (nodes with multiple outgoing edges) and tries to find
            paths that reconverge at a specific "sink" node within a certain distance (`max_len`).

            :param max_len: The maximum length (in edges) of a bubble path.
            :type max_len: int
            :yield: Bubble objects found in the graph.
            :rtype: Iterator[Bubble]
        """
        #generatore
        forks = [n for n in self.nodes.values() if n.is_hub()] #lista nodi hub

        for fork_node in forks:
            candidates = {}  #chiave: sequenza nodo arrivo, valore: lista path(lista di edges)

            for edge in fork_node.outgoing:
                path = [edge] #lista di edge
                curr = edge.destination #nodo arrivo edge
                valid = True #diventa False quando il nodo è un hub o tip

                while len(path) < max_len: #fin quando il path è corto
                    if curr.is_receiving_hub():
                        break #chiusura bolla

                    if len(curr.outgoing) != 1:
                        valid = False #tanto itero per la gestione delle bolle annidate
                        break

                    next_edge = curr.outgoing[0]
                    path.append(next_edge) #aggiungo edge al path
                    curr = next_edge.destination

                #se il percorso è valido e chiuso e non sono tornato al nodo iniziale
                if valid and curr.is_receiving_hub() and curr.seq != fork_node.seq:
                    dest_seq = curr.seq
                    if dest_seq not in candidates: 
                        candidates[dest_seq] = []
                    candidates[dest_seq].append(path) #aggiorno dict


            for dest_seq, paths in candidates.items():
                if len(paths) > 1: #cioè, se quel nodo è davvero un hub, ha più path in entrata:
                    dest_node = self.nodes[dest_seq] #prendo il nodo
                    yield Bubble(fork_node, dest_node, paths) #restituisco l'oggetto bolla

    def process_bubbles(self, cov_threshold: int=5, ratio_threshold: float=0.1, verbose: bool=True) -> int:
        """
            Iteratively finds and resolves bubbles in the graph.
            This method keeps calling  `_detect_bubbles` with maximum len 2k+1 and resolves 
            them until no more bubbles can be resolved.
            It handles nested bubbles by running multiple iterations.

            :param cov_threshold: Max coverage threshold used to distinguish errors from real variants. Defaults to 5
            :type cov_threshold: int
            :param ratio_threshold: Max coverage ratio between paths and the best coverage one to collapse them. Defaults to 0.1.
            :param ratio_threshold: float
            :param verbose: If True, prints details about every iteration of processing bubbles.
            :type verbose: bool
            :return: Total number of bubbles resolved (popped).
            :rtype: int
        """
        if not self.graph_built:
            raise RuntimeError("You cannot process bubbles in a graph if it does not exist")
        if not isinstance(cov_threshold, int):
            raise TypeError("cov_threshold must be integer")
        elif cov_threshold < 0:
            raise ValueError("cov_threshold must be non-negative")
        if not isinstance(ratio_threshold, float):
            raise TypeError("ratio_threshold must be float")
        elif ratio_threshold <= 0 or ratio_threshold > 1:
            raise ValueError("ratio_threshold must be between 0 and 1")
        if not isinstance(verbose, bool):
            raise TypeError("verbose must be boolean")
        max_len = (self.k * 2) + 1
        iteration = 1
        removed_bubbles = 0
        while True: #itero, vengono risolte prima le bolle interne
            popped_count = 0
            found_bubbles = False
            try:
                for bubble in self._detect_bubbles(max_len):
                    found_bubbles = True
                    if bubble._resolve(self, cov_threshold, ratio_threshold):
                        popped_count += 1
            except ValueError:
                print("There have been problems in bubbles detection. You can't perform the operation")
                break
            if not found_bubbles: #sono tutte risolte
                break
            if verbose:
                print(f"Iteration {iteration}: {popped_count} bubbles popped")
            if popped_count == 0: #magari ci sono ancora bolle, ma non vengono risolte
                for bubble in self._detect_bubbles(max_len):
                    self.bubbles[bubble.source_node.seq] = bubble
                break
            removed_bubbles += popped_count
            iteration += 1

        self._rm_isolated_nodes(verbose)
        self.bubbles_processed = True
        return removed_bubbles

    def _detect_tips(self, max_len: int) -> Iterator["Tip"]:
        """
        Generator that searches for dead ends (Tips) in the graph, if shorter than a max_len parameter.

        It scans for:
        1. Type I Tips: Paths starting at a hub and stopping in a dead end.
        2. Type II Tips: Paths starting from nowhere and merging into a hub.

        :param max_len: Maximum length of a tip to be considered. Longer dead ends are kept as valid contigs.
        :type max_len: int
        :yield: Tip objects found in the graph.
        :rtype: Iterator[Tip]
        """
        # Detection tips di tipo I
        potential_inital_nodes_I = [node for node in self.nodes.values() if node.is_hub()] #lista nodi hub
        for fork_node in potential_inital_nodes_I:
            for edge in fork_node.outgoing:
                path = [edge]
                curr = edge.destination
                valid = True #diventa False quando il nodo è un hub

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


    def process_tips(self, cov_threshold: int=5, verbose: bool=True) -> int:
        """
        Iteratively removes dead ends (Tips) from the graph.

        This effectively "shaves" the hairy parts of the graph caused by
        sequencing errors at the ends of reads.

        :param cov_threshold: Coverage threshold to confirm a tip is an error.
        :type cov_threshold: int
        :param verbose: If True, prints details about every iteration of processing tips.
        :type verbose: bool
        :return: Total number of tips removed.
        :rtype: int
        """
        if not self.graph_built:
            raise RuntimeError("You cannot process tips in a graph if it does not exist")
        if not isinstance(cov_threshold, int):
            raise TypeError("cov_threshold must be integer")
        elif cov_threshold < 0:
            raise ValueError("cov_threshold must be non-negative")
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
            if verbose:
                print(f"Iteration {iteration}: {removed_count} tips popped")
            if removed_count == 0: #magari ci sono ancora tips ma non vanno eliminate
                for tip in self._detect_tips(max_len):
                    self.tips[tip.initial_node.seq] = tip
                break
            removed_tips += removed_count
            iteration += 1

        self._rm_isolated_nodes(verbose)
        self.tips_processed = True
        return removed_tips

    def process_short_islands(self, verbose: bool = True) -> int:
        """
        Identifies and removes isolated connected components ("islands") that are
        shorter than a safety threshold.
        This method performs a global cleanup of the graph by finding all disconnected
        subgraphs.

        :param verbose: If True, prints statistics about the removal process to stdout.
        :type verbose: bool
        :return: The number of independent islands removed.
        :rtype: int
        :raises RuntimeError: If the graph has not been built yet.
        """
        if not self.graph_built:
            raise RuntimeError("You cannot process islands in a graph if it does not exist")
        if not isinstance(verbose, bool):
            raise TypeError("verbose must be boolean")

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

            while to_visit: #fin quando to_visit ha un nodo
                curr_node = to_visit.pop(0) #estraggo nodo

                #guardo vicini sia in entrata che in uscita
                all_neighbors = [] #lista di nodi vicini al nodo curr
                for edge in curr_node.outgoing:
                    all_neighbors.append(edge.destination)
                for edge in curr_node.incoming:
                    all_neighbors.append(edge.source)

                for neighbor in all_neighbors:
                    if neighbor.seq not in island_nodes:
                        island_nodes.add(neighbor.seq)
                        visited_nodes.add(neighbor.seq)
                        to_visit.append(neighbor)


            #formula: L = k + (N - 1)
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
            if verbose:
                print(
                    f"Island Removal: Deleted {deleted_islands} short islands (< {min_len}bp). Removed {len(nodes_to_delete)} nodes.")
        else:
            if verbose:
                print("Island Removal: No short islands found.")
        self.isolated_contigs_removed = True
        return deleted_islands

    def is_eulerian_linear(self) -> bool:
        """
        Checks if the current graph satisfies the conditions for a linear Eulerian path.
        
        A linear Eulerian path exists if there is exactly one source (in-degree 0),
        one exit (out-degree 0), and no hubs (branching points).

        :return: True if the graph is a single linear path, False otherwise.
        :rtype: bool
        """
        number_0_in = [node for node in self.nodes.values() if len(node.incoming) == 0]
        number_0_out = [node for node in self.nodes.values() if len(node.outgoing) == 0]
        number_hub = [node for node in self.nodes.values() if node.is_hub()]
        return len(number_0_in) == len(number_0_out) and len(number_0_in) == 1 and len(number_hub) == 0

    def is_eulerian_circular(self) -> bool:
        """
        Checks if the current graph satisfies the conditions for a circular Eulerian circuit.
        
        A circular Eulerian circuit exists if every node has in-degree == out-degree == 1
        and there are no sources, exits, or hubs.

        :return: True if the graph is a single perfect cycle, False otherwise.
        :rtype: bool
        """
        number_0_in = [node for node in self.nodes.values() if len(node.incoming) == 0]
        number_0_out = [node for node in self.nodes.values() if len(node.outgoing) == 0]
        number_hub = [node for node in self.nodes.values() if node.is_hub()]
        return len(number_0_in) == len(number_0_out) and len(number_0_in) == 0 and len(number_hub) == 0

    def _walk_linear(self, start_edge: "Edge") -> tuple[str, "Node", list[int], set[str]]:
        """
        Performs a linear traversal of the graph starting from a specific edge.
        
        The traversal continues as long as nodes have an out-degree of 1 and 
        the next edge has not been visited yet. 
        It assembles linear segments of DNA and stops as soon as there is no 
        more linearity or the path merges with an already visited one.

        :param start_edge: The Edge object where the traversal begins.
        :type start_edge: Edge
        :return: A tuple containing:
            - seq (str): The assembled DNA sequence.
            - curr (Node): The final node reached (where the path forks or ends).
            - walk_weights (list[int]): List of coverage weights for each edge walked.
            - receiving_hubs_visited (set[str]): K-mers of convergence (receiving hubs) nodes encountered.
        :rtype: tuple
        """
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

    def _checks_and_create_starts(self, fork_in_hub: bool) -> list["Node"]:
        """
        Validates graph state and identifies starting points for assembly.
        
        Gives warning if bubbles, tips, and isolated contigs have not been processed 
        before starting. If fork_in_hub is True, only nodes with in-degree 0 
        are considered starts. Otherwise, internal hubs are also included because walking 
        linearly will stop at each hub.

        :param fork_in_hub: If True, restricts starting points to source nodes (In-degree 0).
        :type fork_in_hub: bool
        :raises RuntimeError: If the graph has not been built.
        :return: A list of Node objects from which assembly will initiate.
        :rtype: list[Node]
        """
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
        return all_starts

    def _assemble_bubble(self, bubble_to_assemble: "Bubble", common_contig: "Contig") -> list["Contig"]:
        """
        Resolves a bubble structure by integrating all its paths with the current 
        preceeding contig.
        
        For each path in the bubble, it creates a new Contig instance that 
        combines the prefix (common_contig), the bubble branch, and any 
        sequence following the bubble's destination node using the solid addition
        logic of the contig class.
        It creates up to n_bubble_paths * n_final_node_paths contigs 

        :param bubble_to_assemble: The Bubble object to resolve.
        :type bubble_to_assemble: Bubble
        :param common_contig: The Contig instance preceding the bubble.
        :type common_contig: Contig
        :return: A list of fully assembled and tagged Contig objects (one for each of the possible contigs).
        :rtype: list[Contig]
        """
        contigs = []
        for path in bubble_to_assemble.paths: #per ciascun path della bolla
            for edge in path: #per ciascun edge nel path
                edge.visited = True
        alternative_paths = list(bubble_to_assemble._get_contigs_from_path(self.k)) #lista di oggetti contig con sequenza e lista pesi
        new_start_node = bubble_to_assemble.dest_node
        bubble_to_assemble.visited = True
        if len(new_start_node.outgoing) == 0 or all([edge.visited for edge in new_start_node.outgoing]): #cioè, se non continua dopo la chiusura della bolla
            for alt_contig in alternative_paths:
                whole_contig = common_contig + alt_contig #magicmtethod _add_ di Contig
                whole_contig.get_tag()
                contigs.append(whole_contig)
        else: #nodo fine bolla non è terminante... può essere un hub
            if new_start_node.is_hub():
                new_hub_seq = new_start_node.seq
            else:
                new_hub_seq = None
            for new_start_edge in new_start_node.outgoing:
                if new_start_edge.visited:
                    continue
                new_contig = Contig.from_walk_linear(self._walk_linear(new_start_edge), self.k, new_start_node, hub = new_start_node.is_hub(), hub_seq = new_hub_seq) #continua a camminare dopo la chiusura della bolla
                for alt_contig in alternative_paths:
                    whole_contig = common_contig + alt_contig + new_contig
                    whole_contig.get_tag()
                    contigs.append(whole_contig)
        for remaining_edge in [edge for edge in common_contig.final_node.outgoing if not edge.visited]: # Edges che partono da next_node ma non sono della bolla
            remaining_contig = Contig.from_walk_linear(self._walk_linear(remaining_edge), self.k, common_contig.final_node, hub = True, hub_seq = bubble_to_assemble.source_node.seq)
            whole_contig = common_contig + remaining_contig
            whole_contig.get_tag()
            contigs.append(whole_contig)
        return contigs

    def assemble_contigs(self, fork_in_hub: bool = False) -> list["Contig"]:
        """
        The full assembly engine. Orchestrates the global traversal of the graph.
        
        The process follows four main stages:
        
        1. Primary traversal from identified start nodes, stopping as soon as a linear contig ends if
           fork_in_hub is false, otherwise trying to solve it.
        2. Structural resolution (Bubbles/Tips/simple forks) if an intersection is met and fork_in_hub is True.
        3. Residual hub processing (if fork_in_hub is True) to capture segments starting in secondary hubs.
        4. "Catch-all" circular search for isolated loops (circular DNA elements). It keeps going until all edges are visited

        All generated contigs are automatically classified (tagged), 
        deduplicated, and sorted by length for easier stats calculation and visualization.

        :param fork_in_hub: If False, breaks assembly at hub nodes to increase accuracy.
        :type fork_in_hub: bool
        :return: A list of unique, topologically classified Contig objects.
        :rtype: list[Contig]
        """
        all_starts = self._checks_and_create_starts(fork_in_hub) #lista nodi di inizio
        if not self.nodes or not self.edges:
            warnings.warn("Empty graph: no nodes or edges to assemble")
            self.contigs_assembled = True
            return []

        contigs = []
        contigs_forked_in_hub = False
        if self.is_eulerian_linear() or self.is_eulerian_circular():
            expected = "circular" if self.is_eulerian_circular() else "linear"
            print(f"Only a single eulerian path found. Expected a single {expected} contig")
        for contig_start in all_starts:
            for start_edge in contig_start.outgoing:
                if start_edge.visited: #strada già percorsa
                    continue
                #creazione oggetto conting con il classmethod, uso _walk.linear che cammina fin quando trova bivio o finisce la strada
                contig = Contig.from_walk_linear(self._walk_linear(start_edge), self.k, contig_start)
                forked = False
                if contig.final_node.is_hub() and fork_in_hub and not (contig.is_lariat() or contig.is_circular()):
                    contigs_forked_in_hub = True
                    forked = True
                    common_contig = contig.copy()
                    if contig.final_node.seq in self.bubbles and not self.bubbles[contig.final_node.seq].visited: #è l'inizio di una bolla non già visitata
                        bubble_to_assemble = self.bubbles[contig.final_node.seq] #prendo la bolla
                        bubble_contigs = self._assemble_bubble(bubble_to_assemble, common_contig) #lista di contigs
                        contigs.extend(bubble_contigs) #extend altrimenti append inserirebbe una lista
                    elif contig.final_node.seq in self.tips and not self.tips[contig.final_node.seq].visited:
                        tip_to_assemble = self.tips[contig.final_node.seq]  #prendo la tip
                        for edge in tip_to_assemble.path:
                            edge.visited = True
                        tip_contig = tip_to_assemble._get_contig_from_path(self.k) #creo l'oggetto contig
                        tip_to_assemble.visited = True
                        whole_contig = common_contig + tip_contig #grazie al magicmethod
                        whole_contig.get_tag()
                        contigs.append(whole_contig)
                    else: #non è una bolla o una tip
                        new_start_node = contig.final_node
                        for new_start_edge in new_start_node.outgoing:
                            if new_start_edge.visited:
                                continue
                            new_contig = Contig.from_walk_linear(self._walk_linear(new_start_edge), self.k, new_start_node, hub = True, hub_seq = common_contig.final_node.seq) #creo contig
                            whole_contig = common_contig + new_contig
                            whole_contig.get_tag()
                            contigs.append(whole_contig)
                if not forked: #non ho trovato nessun bivio, oppure fork_in_hub spento
                    # Marco che ho trovato un contig che inizia in una bolla o un tip
                    if contig.start_node.seq in self.bubbles:
                        if self.bubbles[contig.start_node.seq].dest_node.seq in contig.seq: #controllo che abbia percorso tutta la bolla
                            contig.bubble = True
                            contig.bubble_ends = (contig.start_node.seq, self.bubbles[contig.start_node.seq].dest_node.seq)
                    elif contig.start_node.seq in self.tips:
                        if self.tips[contig.start_node.seq].final_node.seq in contig.seq:
                            contig.tip = True
                            contig.tip_hub = contig.start_node.seq
                    contig.get_tag()
                    contigs.append(contig)

        if fork_in_hub:
            #hubs restanti poiché fork_in_hub si ferma dopo un livello
            #set di hub con almeno un arco non esplorato
            remaining_hubs_seqs = {node.seq for node in self.nodes.values() if len(node.incoming) > 0 and node.is_hub() and any(not edge.visited for edge in node.outgoing)}
            while remaining_hubs_seqs:
                for contig_start_seq in list(remaining_hubs_seqs):
                    contig_start = self.nodes[contig_start_seq] #nodo inizio
                    for start_edge in contig_start.outgoing:
                        if start_edge.visited:
                            continue
                        remaining_contig = Contig.from_walk_linear(self._walk_linear(start_edge), self.k, contig_start, remaining = True) #creo contig
                        if remaining_contig.start_node.seq in self.bubbles:
                            if self.bubbles[remaining_contig.start_node.seq].dest_node.seq in remaining_contig.seq:
                                remaining_contig.bubble = True
                                remaining_contig.bubble_ends = (remaining_contig.start_node.seq, self.bubbles[remaining_contig.start_node.seq].dest_node.seq)
                        elif remaining_contig.start_node.seq in self.tips:
                            if self.tips[remaining_contig.start_node.seq].final_node.seq in remaining_contig.seq:
                                remaining_contig.tip = True
                                remaining_contig.tip_hub = remaining_contig.start_node.seq
                        remaining_contig.get_tag()
                        contigs.append(remaining_contig)
                    if all(edge.visited for edge in contig_start.outgoing): #controllo che abbiamo visitato tutte le uscite dall'hub
                        remaining_hubs_seqs.remove(contig_start_seq) #se si, rimuovo

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
                start_edge = not_visited_edges[0] #un arco a caso
                start_node = start_edge.source
                circ_contig = Contig.from_walk_linear(self._walk_linear(start_edge), self.k, start_node, from_catch_all = True) #creo contig
                circ_contig.get_tag()
                contigs.append(circ_contig)
                if circ_contig.tag == "circular_isolated":
                    circular_contigs += 1
                else:
                    print(f"Contig found with unknown topology {len(circ_contig)}")
                    unknown_contigs += 1
                not_visited_edges = [edge for edge in self.edges.values() if not edge.visited] #aggiorno la lista
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
        deduplicated_contigs.sort(key=lambda x: len(x.seq), reverse=True) #ordina dal più lungo al più corto
        self.contigs_forked_in_hub = contigs_forked_in_hub
        self.contigs = deduplicated_contigs
        self.contigs_assembled = True
        if contigs_forked_in_hub:
            print("Some contigs have been forked in hub. N50 may not be totally accurate")
        print(f"Assembled a total of {len(deduplicated_contigs)} contigs")
        return deduplicated_contigs
