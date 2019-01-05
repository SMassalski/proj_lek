import objects


def read_drugs(file):

    drugs = {}
    o = open(file, 'r')
    o.readline()  # header
    for line in o:
        line = line.strip().split('\t')
        dbid = line[0]
        name = line[1]
        groups = line[3].split('|')
        categories = line[5].split('|')
        drugs[dbid] = objects.Drug(dbid, name, groups, categories)
    return drugs


def read_proteins(file, drugs):

    proteins = {}
    edges = []
    o = open(file, 'r')
    o.readline()  # header
    for line in o:
        line = line.strip().split('\t')
        upid = line[2]
        organism = line[4]
        if organism.strip()[-1] == ')':
            org = organism.split(' ')
            organism = ''
            for el in org:
                if not el.startswith('('):
                    organism += el
                else:
                    break
        category = line[1]
        drug = drugs[line[0]]
        protein = objects.Protein(upid, organism)
        proteins[upid] = protein
        edge = objects.Edge(drug, protein, category)
        edges.append(edge)
        drug.add_edge(edge)
    return proteins, edges


directory = '/home/marni/Dokumenty/projektowanie_lekow/git_drugbank/data/'
drugs = read_drugs('%sdrugbank.tsv' % directory)
proteins, edges = read_proteins('%sprotein.tsv' % directory, drugs)
