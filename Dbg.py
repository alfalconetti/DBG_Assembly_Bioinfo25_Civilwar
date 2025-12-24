from kmer import kmer_generate

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
    def __init__(self):
        self.nodes = {}  #dizionario nodi
        self.edges = {}  #dizionario archi
		
