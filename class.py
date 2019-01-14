import util
import matplotlib.pylab as plt
import numpy as np


def superclasses(drugs):

    camount = {}
    cdrugs = {}
    ctargets = {}
    for drug in drugs.values():
        if not drug.classification:
            cc = 'None'
        else:
            cc = drug.classification[1]
        if cc not in camount:
            camount[cc] = 1
        else:
            camount[cc] += 1
        if cc not in cdrugs:
            cdrugs[cc] = [drug]
        else:
            cdrugs[cc].append(drug)
        if cc not in ctargets:
            ctargets[cc] = {}
        for edge in [e[0] for e in drug.edges]:
            if edge not in ctargets[cc]:
                ctargets[cc][edge] = 1
            else:
                ctargets[cc][edge] += 1
    try:
        del (ctargets[''])
        del (cdrugs[''])
    except KeyError:
        pass
    return camount, cdrugs, ctargets


def plot_classes(ctargets):

    #plt.figure(figsize=(20, 15))
    font = {'size': 20}
    plt.rc('font', **font)

    names = []
    types = {'target': [], 'enzyme': [], 'transporter': [], 'carrier': []}
    for group in ctargets.keys():
        if group is None:
            names.append('None')
        else:
            names.append(group)
        for t in types.keys():
            if t in ctargets[group]:
                types[t].append(ctargets[group][t])
            else:
                types[t].append(0)
    ind = np.arange(len(names))

    bot = [0 for i in range(len(names))]
    plots = []
    leg = []
    for t in types.keys():
        p = plt.bar(ind, types[t], bottom=bot)
        bot = [b + v for b, v in zip(bot, types[t])]
        plots.append(p)
        leg.append(t)

    plt.title("Superclasses of drugs depending on type of drugs' targets", fontsize=20)
    plt.legend(plots, leg, fontsize=16)

    plt.ylabel('Number of edges', fontsize=12)
    plt.xticks(ind, names, rotation='vertical', fontsize=12)
    # plt.yticks(np.arange(0, 81, 10))
    plt.yticks(fontsize=12)
    plt.show()


def plot_classdeg(cdrugs, average):
    cand = {}
    for c in cdrugs.keys():
        deg = 0
        for i, d in enumerate(cdrugs[c]):
            deg += len(d.edges)
        cand[c] = deg/(i+1)

    xx = np.arange(len(cand))
    plt.scatter(xx, cand.values(), marker='o')
    plt.plot([0, len(cand)], [average]*2, '--', c='g')
    plt.xticks(xx, cand.keys(), rotation='vertical', fontsize=12)
    plt.ylabel('Average degree of nodes', fontsize=12)
    plt.yticks(fontsize=12)
    plt.title('Average degree of nodes in each superclass', fontsize=20)
    plt.show()


drugs, proteins = util.parse('/home/marni/Dokumenty/projektowanie_lekow/full_database.xml')
camount, cdrugs, ctargets = superclasses(drugs)
plot_classes(ctargets)
plot_classdeg(cdrugs, 3.44)
