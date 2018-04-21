#!/usr/bin/env python

# this is a python replacement and hopefully improvement
# for plotvafRepeatability.R
def builddf(commonVars, avgData, principleData):
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns

    # combine all date into a pandas df
    df = pd.DataFrame({'x':principleData, 'y':avgData})

    # deal with unique variants
    '''
    if commonVars:
        print 'common'
    else:
        df = df.dropna()
    '''
    df = df[df['x'] < 0.1]
    df = df[df['y'] < 0.1]

    # plot data
    plt.scatter(df['x'], df['y'], color = 'gray', s=5) # s controls point size
    plt.plot([0.14,0.86],[1,1], linewidth=22)
    plt.xlim(0,0.02)
    plt.ylim(0,0.02)

    sns.set_context("paper", font_scale=3)
    plt.xlabel('Sample x', {'size':'20'})
    plt.ylabel('Sample y', {'size':'20'})
    plt.title('VAF Comparison')

    sns.despine(offset=10, trim=True)
    plt.tight_layout()
    plt.show()
    #plt.bar(range(len(means)), means, color=colors)






if __name__ == '__main__':
    builddf(commonVars, avgData, principleData)

