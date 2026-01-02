from kmer import kmer_generate
from node import Node
from edge import Edge
from parser import parser
from utils import rc

class Dbg:
    # DUBBI: FILTERING PER BASSA QUALITY?


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
        self.k = k #lunghezza kmeri
        self.nodes = {}  #dizionario nodi
        self.edges = {}  #dizionario archi
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


    def build_graph(self, quality_check = False, quality_threshold = 30, use_both_strands = True):
        
        def median_qual(qual_string):
            from statistics import median
            qual_list = [ord(x) - 33 for x in qual_string]
            return median(qual_list)


        def put_reads_graph(reads, quality_check, use_both_strands, reverse_complement = False):
            for read_tup in reads:
                read = read_tup[0]
                if reverse_complement:
                    read = rc(read)
                if quality_check:
                    read_qual = read_tup[1]
                    if median_qual(read_qual) > quality_threshold:
                        self.add_read(read, use_both_strands = use_both_strands)
                else:
                    self.add_read(read, use_both_strands = use_both_strands)
    

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

        
    def add_read(self, read, use_both_strands):
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

    
    def prune_graph(self, threshold=2):
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
        

    def get_graph_stats(self):
    
        num_nodes = len(self.nodes) #numero nodi
        num_edges = len(self.edges)
        total_weight = 0
        
        if num_edges > 0:
            total_weight = sum(edge.weight for edge in self.edges.values())
            
            avg_coverage = total_weight / num_edges
        else:
            avg_coverage = 0
            
        print("Graph Statistics: ")
        print(f"Nodes: {num_nodes}")
        print(f"Edges: {num_edges}")
        print(f"Average coverage: {avg_coverage:.2f}X")
        
        
    def write_contigs(self, output_filename="contigs.fasta"):
        if not isinstance(output_filename, str):
            raise TypeError("output_filename must be string")
        contigs = self.get_contigs()  #get_contigs da fare
        try:
            with open(output_filename, "w") as f:
                i=0
                for contig in contigs:
                    f.write(f">contig_{i+1} len={len(contig)}\n")
                    f.write(f"{contig}\n")
                    i+=1
            print(f"saved file: {output_filename}")
            return [len(c) for c in contigs] #lista lunghezze per n50
            
        except IOError as e:
            print(f"Error writing to file {output_filename}: {e}")
            return[]

    
    @staticmethod
    def n50(lengths): #lista di lunghezze
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
                
    def rm_isolated_nodes(self):
        nodes_to_rm = []
        for kmer, node in self.nodes.items():
            if len(node.incoming) == 0 and len(node.outgoing) == 0:
                nodes_to_rm.append(kmer)
        
        for kmer in nodes_to_rm:
            del self.nodes[kmer]
        
        print(f"removed {len(nodes_to_rm)} nodes")
        
