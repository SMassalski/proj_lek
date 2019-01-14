import matplotlib.pylab as plt
import util


def plot_date_location(drugs, proteins, normalised=False, linear=False):
    loc = set()
    for p in proteins.values():
        loc.add(p.subcell_location)
    names = {'Unknown/other': ['Microsome', 'Glycosome', 'Peroxisome', 'Cell junction', 'Golgi apparatus', 'reticulum',
                               'Early endosome', 'Plastid', 'Melanosome', 'Vacuole', 'Cell projection', 'Lysosome',
                               'Vacuole', 'Fimbrium', 'Virion', 'Midbody', '', 'Cytoplasmic vesicle membrane',
                               'Golgi apparatus membrane', 'Virion membrane', 'Lysosome membrane',
                               'Early endosome membrane', 'Parasitophorous vacuole membrane',
                               'Sarcoplasmic reticulum membrane', 'Mitochondrion outer membrane', 'Melanosome membrane',
                               'Endoplasmic reticulum membrane', 'Late endosome membrane', 'Golgi apparatus lumen',
                               'Endoplasmic reticulum', 'Apical cell membrane', 'Rough endoplasmic reticulum',
                               'Sarcoplasmic reticulum lumen', 'Microsome membrane', 'Cellular thylakoid membrane',
                               'Host endoplasmic reticulum membrane', 'Endosome membrane',
                               'Endoplasmic reticulum-Golgi intermediate compartment membrane',
                               'Endoplasmic reticulum-Golgi intermediate compartment',
                               'Rough endoplasmic reticulum membrane', 'Endoplasmic reticulum lumen',
                               'Recycling endosome membrane', 'Endomembrane system', 'Peroxisome membrane',
                               'Vacuole membrane'],
             'Plasma membrane': ['Nucleus membrane', 'Membrane', 'Periplasm', 'Host membrane',
                                 'Cellular chromatophore membrane', 'Cell outer membrane', 'Cell membrane',
                                 'Cell inner membrane; multi-pass membrane protein', 'Cell inner membrane',
                                 'Cell inner membrane; peripheral membrane protein', 'Host cell membrane',
                                 'Basolateral cell membrane', 'Lateral cell membrane', 'Basal cell membrane'],
             'Cytosol': ['Cytoplasmic', 'cytoplasmic', 'Cytoplasm', 'cytoplasm', 'Host cytoplasm',
                         'Cytoplasmic granule', 'Cytoplasm (By similarity)'],
             'Extracellular space': ['Cell surface', 'Secreted'],
             'Mitochondrion': ['Mitochondrion', 'Mitochondrion matrix', 'Mitochondrion membrane',
                               'Mitochondrion inner membrane', 'Mitochondrion intermembrane space'],
             'Nucleus': ['Nucleus', 'nucleus', 'Chromosome', 'Nucleus membrane', 'Nucleus speckle', 'Nucleus envelope',
                         'Nucleus inner membrane', 'Nucleus matrix', 'Host nucleus', 'Nucleus outer membrane']
             }

    groups = {'Unknown/other': 0, 'Plasma membrane': 0, 'Cytosol': 0, 'Nucleus': 0, 'Extracellular space': 0,
              'Mitochondrion': 0}

    locations = {}
    for i, drug in enumerate(drugs.values()):
        if drug.patent_date is None:
            continue
        year = drug.patent_date[:4]
        if year not in locations:
            locations[year] = groups.copy()
        for e in [e[1] for e in drug.edges]:
            if isinstance(e, str):
                cellloc = ''
            else:
                cellloc = e.subcell_location
            for l in names.keys():
                if cellloc in names[l]:
                    locations[year][l] += 1

    years = list(locations.keys())
    for y in years:
        if sum(locations[y].values()) == 0:
            locations.pop(y)

    years_sort = list(map(int, list(locations.keys())))
    years_sort.sort()

    normloc = {}
    for l in locations.keys():
        s = sum(locations[l].values())
        normloc[l] = {}
        for ll in locations[l].keys():
            normloc[l][ll] = locations[l][ll] / s

    if linear:

        for g in groups.keys():
            xx = []
            yy = []
            for year in years_sort:
                xx.append(year)
                yy.append(normloc[str(year)][g])
            plt.plot(xx, yy, label=g)

        plt.legend(loc='upper right', fontsize=16)
        plt.ylabel('Normalised number of targets', fontsize=12)

    else:

        plots = []
        leg = []
        yy = [0 for i in range(len(locations))]

        for g in groups.keys():
            h = []
            for year in years_sort:
                if not normalised:
                    h.append(locations[str(year)][g])
                else:
                    h.append(normloc[str(year)][g])
            p = plt.bar(years_sort, h, bottom=yy)
            yy = [y + v for y, v in zip(yy, h)]
            plots.append(p)
            leg.append(g)

        plt.legend(plots, leg, fontsize=16)

        if not normalised:
            plt.ylabel('Number of targets', fontsize=12)
        else:
            plt.ylabel('Normalised number of targets', fontsize=12)

    plt.title("Cellular location of the target depending on the drug's patent year", fontsize=20)
    plt.yticks(fontsize=12)
    plt.xticks(years_sort, rotation='vertical', fontsize=12)
    plt.show()


drugs, proteins = util.parse('/home/marni/Dokumenty/projektowanie_lekow/full_database.xml')

# Plot bar chart
plot_date_location(drugs, proteins)

# Plot normalised bar chart
plot_date_location(drugs, proteins, normalised=True)

# Plot standard linear plot
plot_date_location(drugs, proteins, linear=True)
