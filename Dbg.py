from kmer import kmer_generate
from node import Node
from edge import Edge

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
    def __init__(self, k):
		self.k = k #lunghezza kmeri
        self.nodes = {}  #dizionario nodi
        self.edges = {}  #dizionario archi

	def add_read(self, read):
		list_of_kmers = kmer_generate(read, self.k) #generazione kmeri
		
		for kmer in list_of_kmers:
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
		
