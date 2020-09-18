import pandas as pd
import argparse
import matplotlib
import matplotlib.pyplot as plt



Metrics = ["fname",  # Metrics we used for Profiling
           "ColsPerGroup",
           "n",
           "nnz",
           "nnzPerRow",
           "ngroup",
           "nnz_access",
           "nnz_reuse",
           "Flops",
           "nlevels",
           "n_sys",
           "aver_parallelism",
           "nnzPreLevels",
           "Aver_MaxDiff",
           "Var_MaxDiff",
           "Sum_MaxDiff",
           "Num_threads",
           "t_serial",
           "t_group_level",
           "t_level",
           ""]

def ProcessCSV(csvfile):
    data = pd.read_csv(csvfile, sep=',', header=None)
    print(len(Metrics))
    data.columns = Metrics
    return data

def GetData(csvdata, fname, xlabel, ylabels, yname):
     sdata=csvdata.loc[csvdata['fname']==fname]
     print(sdata)

     xdata = sdata[xlabel]
     print(xdata)


     ydata = [sdata[y] for y in ylabels ]

     for y in ydata:
        plt.plot(xdata, y, '+-')

     plt.xlabel(xlabel)
     # plt.xticks(xdata)
     plt.ylabel(yname)
     plt.yscale('log')
     plt.title("Results on " + fname)
     plt.legend()
     plt.show()

def PlotOneMatrix(fname, data):
    xlabel_="ColsPerGroup"
    ylabel_="Execution Time (seconds)"

    metrics=["t_serial","t_group_level", "t_level"]

    GetData(data, fname, xlabel_, metrics, ylabel_)

    metrics=["nlevels"]
    ylabel_= "Number of Levels"

    GetData(data, fname, xlabel_, metrics, ylabel_)

    metrics=["Sum_MaxDiff"]
    ylabel_= "Sum_MaxDiff"

    GetData(data, fname, xlabel_, metrics, ylabel_)
def PlotAllMatrix(data):
    fdata=data
    metrics=["fname", "t_serial","t_group_level", "t_level"]
    data = data[metrics]

    data['t_group_level'] = data['t_group_level']/data['t_serial']
    data['t_level'] = data['t_level']/data['t_serial']
    data['t_serial'] = data['t_serial']/data['t_serial']

    data=data.set_index('fname')

    print(data)

    ax = data.plot.bar(rot=75)
    # plt.yscale('log')
    plt.ylabel("Normalized execution time")
    plt.show()

    metrics=["fname", "nnz_access", "nnz_reuse"]
    data = fdata[metrics]
    data=data.set_index('fname')


    data['reuse_ratio']=data['nnz_reuse']/data['nnz_access']
    data = data['reuse_ratio']

    print(data)
    ax = data.plot.bar(rot=75)
    plt.ylabel("Reuse_nnz/ Access_nnz")
    plt.show()





if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--f", help="directory to ")

    args = parser.parse_args()

    data = ProcessCSV(args.f)
    data.to_csv("splocality.csv", sep=',')

    # fname="nd24k.mtx"
    #
    # PlotOneMatrix(fname, data)

    # PlotAllMatrix(data.loc[0:33])
    # PlotAllMatrix(data.loc[34:67])