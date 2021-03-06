import xml.etree.ElementTree as ET


class Drug:
    def __init__(self, dbid, name, groups, classification, patent_date, smiles):
        self.dbid = dbid
        self.smiles = smiles
        self.name = name
        self.edges = []
        self.classification = classification
        self.groups = groups
        self.patent_date = patent_date


class Protein:
    def __init__(self, upid, name, location, processes):
        self.upid = upid
        self.name = name
        self.subcell_location = location
        self.bio_processes = processes


def parse(xml_path):
    xml_file = open(xml_path, 'r')
    tree = ET.iterparse(xml_file, events=('start', 'end'))

    ns = '{http://www.drugbank.ca}'

    rows = []
    drug = None
    drugs = {}
    proteins = {}

    for event, element in tree:
        if event == 'start' and element.tag == ns + 'drug' and not drug:
            drug = element
        elif event == 'end' and element.tag == ns + 'drug' and element is drug:
            dbid = drug.findtext(ns + "drugbank-id[@primary='true']")
            name = drug.findtext(ns + "name")
            groups = [group.text for group in drug.findall("{ns}groups/{ns}group".format(ns=ns))]
            levels = ['kingdom', 'superclass', 'class', 'subclass', 'direct-parent']
            classification = [drug.findtext("{ns}classification/{ns}{level}".format(ns=ns, level=level)) for level in
                              levels]
            smiles = drug.findtext("{ns}calculated-properties/{ns}property[{ns}kind='SMILES']/{ns}value".format(ns=ns))
            patents = [date.text for date in drug.findall("{ns}patents/{ns}patent/{ns}approved".format(ns=ns))]
            if len(patents) == 0:
                patent_date = None
            else:
                patent_date = min(patents, key=lambda date: int(''.join(date.split('-'))))
            d = Drug(dbid, name, groups, classification, patent_date, smiles)
            drugs[dbid] = d

            # Proteins
            for category in ['target', 'enzyme', 'carrier', 'transporter']:

                targets = drug.findall('{ns}{cat}s/{ns}{cat}'.format(ns=ns, cat=category))
                for t in targets:
                    polypeptides = t.findall('{ns}polypeptide'.format(ns=ns))
                    for polyp in polypeptides:
                        upid = polyp.findtext(
                            "{ns}external-identifiers/{ns}external-identifier[{ns}resource='UniProtKB']/{ns}identifier".format(
                                ns=ns))

                        if upid not in proteins:
                            name = polyp.findtext(ns + "name")
                            location = polyp.findtext(ns + 'cellular-location')
                            bio_processes = set([process.text for process in polyp.findall(
                                "{ns}go-classifiers/{ns}go-classifier[{ns}category='process']/{ns}description".format(
                                    ns=ns))])

                            proteins[upid] = Protein(upid, name, location, bio_processes)

                        d.edges.append((category, proteins[upid]))
                    if len(polypeptides) == 0:
                        name = t.findtext(ns + "name")
                        d.edges.append((category, name))

            drug.clear()
            drug = None
    xml_file.close()
    return drugs, proteins


def write_tsv(drug_list, filename):

    f = open(filename,'w')
    f.write('\t'.join(['drug', 'type', 'target'])+'\n')
    for drug in drug_list:
        for edge in drug.edges:
            target = edge[1] if isinstance(edge[1], str) else edge[1].upid
            f.write('\t'.join([drug.dbid, edge[0], target])+'\n')
    f.close()
