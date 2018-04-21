#!/usr/bin/env python

# this is a python replacement and hopefully improvement
# for plotvafRepeatability.R
def builddf(commonVars, avgData, principleData):
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns

    # combine all date into a pandas df
    df = pd.DataFrame({'x':principleData, 'y':avgData})

    # this will retain the uniques by converting NaN to 0
    if not commonVars:
        df = df.fillna(value=0)

    # this should be defined at runtime
    df = df[df['x'] < 0.1]
    df = df[df['y'] < 0.1]

    # plot data
    plt.scatter(df['x'], df['y'], color = 'gray', s=5) # s controls point size
    plt.xlim(-0.0002,0.02)
    plt.ylim(-0.0002,0.02)

    # add y=x line
    plt.plot([0,1],[0,1], lw=2, color='black')

    sns.set_context("paper", font_scale=2)
    plt.xlabel('Sample x', {'size':'20'})
    plt.ylabel('Sample y', {'size':'20'})
    plt.title('VAF Comparison')

    sns.despine(offset=10, trim=True)
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    builddf(commonVars, avgData, principleData)

