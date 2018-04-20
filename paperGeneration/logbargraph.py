#!/usr/bin/env python

def logplota():
    import matplotlib.pyplot as plt
    import numpy as np
    import seaborn as sns

    # plot 1
    #data = [0.5056852381, 0.008218755132, 0.001361672189, 0.0003120390789]
    #115227896
    donor = [0.506773,0.498897,0.508709]
    # all points
    a = [0.002194,0.002193,0.002223,0.007349,0.002556,0.002735,0.002592,0.006216]
    # without high freq probes
    a = [0.002194,0.002193,0.002223,0.002556,0.002735,0.002592]
    b = [0.000204,0.000291]
    c = [0.0001029646784]
    labels = ['Donor', '1:500', '1:10000','Recipient']
    sns.set_context("paper", font_scale=2)

    plt.yscale('log')
    plt.xlabel("Sample", {'size':'20'})
    plt.ylabel("VAF log10", {'size':'20'})
    plt.title("VAF Dilutions", {'size':'20'})

    y = [np.mean(donor), np.mean(a), np.mean(b),np.mean(c)]
    err = [np.std(donor),np.std(a),np.std(b),np.std(c)]
    x = np.arange(len(y))
    #sns.barplot(x, y, palette='deep')
    sns.barplot(x, y, color='grey')
    plt.errorbar(x, y, err, linestyle='None', color='#575859')
    plt.xticks(range(len(y)),labels)
    plt.ylim(0.00001,1)

    #sns.despine()
    sns.despine(offset=10, trim=True)
    plt.tight_layout()
    plt.savefig('1d.png', bbox_inches='tight',dpi=1200)
    plt.show()

def logplotb():
    import matplotlib.pyplot as plt
    import numpy as np
    import seaborn as sns
    # plot 2
    # single dilutions
    a = [0.5056852381,0.51115195,0.50577619]
    b = [0.00850932521,0.007668253692,0.008478686493]
    c = [0.0014687054,0.001470865448,0.001145445721]
    d = [0.0003120390789,0.0001721584978,0.0002848341876]
    e = [0.00021234]
    labels = ['Donor', '1:50', '1:500', '1:5000','Recipient']

    plt.yscale('log')
    plt.xlabel("Sample", {'size':'20'})
    plt.ylabel("VAF log10", {'size':'20'})
    plt.title("VAF Dilutions", {'size':'20'})

    y = [np.mean(a), np.mean(b), np.mean(c), np.mean(d),np.mean(e)]
    err = [np.std(a),np.std(b),np.std(c), np.std(d),np.std(e)]
    x = np.arange(len(y))
    sns.barplot(x, y, color='grey')
    plt.errorbar(x, y, err, linestyle='None', color='#575859')
    plt.xticks(range(len(y)),labels)
    plt.ylim(0.0001,1)

    sns.set_context("paper", font_scale=2)
    sns.despine(offset=10, trim=True)
    #sns.despine()
    plt.tight_layout()
    plt.savefig('1c.png', bbox_inches='tight',dpi=1200)
    plt.show()
