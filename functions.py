import sys
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.fftpack import fft
from scipy import signal


#def load(fname):
#    f = open(fname,'r')
#    data = []
#    for line in f.readlines()[:-1]:
#        data.append(line.replace('\n','').split(' '))
#    f.close()
#    return zip(*data)

##def column(matrix,i):
##    arr = np.array(matrix[i]).transpose()
##    return arr
##    #return [row[i] for row in matrix]

def load_initial():
    f = open("./0/include/initialConditions", 'r').readlines()
    f_l = len(f)
    i = 0
    for i in range(f_l):
        if f[i].strip().startswith('omega'):
            f[i] = f[i].replace("\t"," ").replace("(","").replace(")","").replace(";","")
            omega = " ".join(f[i].split()).split(" ")
        if f[i].strip().startswith('current'):
            f[i] = f[i].replace("\t"," ").replace("(","").replace(")","").replace(";","")
            current = " ".join(f[i].split()).split(" ")
        if f[i].strip().startswith('waveNumber'):
            f[i] = f[i].replace("\t"," ").replace("(","").replace(")","").replace(";","")
            k = " ".join(f[i].split()).split(" ")
    i +=1
    return current[1],omega[1],k[1]

def query_yes_no(question, default="yes"):
    valid = {"yes":"yes",   "y":"yes",  "ye":"yes",
             "no":"no",     "n":"no"}
    if default == None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)

    while 1:
        sys.stdout.write(question + prompt)
        choice = raw_input().lower()
        if default is not None and choice == '':
            return default
        elif choice in valid.keys():
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' "\
                             "(or 'y' or 'n').\n")

def plotgraph(x,y,y2,pp,title,xlabel,ylabel,rp):
    if all(v==0 for v in y):
        print('There is nothing in "%s"' % title)
    else:
        plt.title(title)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.plot(x[rp[0]:rp[1]],y[rp[0]:rp[1]])
        if y2 != 0: ### <>
            a = np.empty(np.size(x[rp[0]:rp[1]]))
            a.fill(y2)
            plt.plot(x[rp[0]:rp[1]],a)
        if pp == 1:
            peakFind = signal.find_peaks_cwt(y,np.arange(0.01,1))
            plt.plot(x[peakFind],y[peakFind],'ro')
        plt.show()

def plotxy(x,y,title,xlabel,ylabel):
    if all(v==0 for v in y):
        print('There is nothing in "%s"' % title)
    else:
        plt.title(title)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.plot(x,y)
        plt.show()

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def file_len(inputFile):
    with open(inputFile) as f:
        for i,l in enumerate(f):
            pass
    return i

def num_iter(inputFile,f_l):
    f = open(inputFile, 'r').readlines()
    i = k = 0
    for i in range(f_l):
        if f[i].startswith('Time = '):
            k += 1
    return k

def calc_mean(x,y,rp):
    x = x[rp[0]:rp[1]]
    y = y[rp[0]:rp[1]]
    i = a0 = 0
    N = len(y)-1
    for i in range(N):
        dx = float(x[i+1])-float(x[i])
        dy = 0.5*(float(y[i+1])+float(y[i]))
        X = float(x[N])-float(x[0])
        a0  = a0 + 1.0/X*dx*dy
        i += 1
    return a0

def calc_mean_ampl(x,y,rp):
    y_mean = calc_mean(x,y,rp)
    y = y[rp[0]:rp[1]]-y_mean
    y_sq = y**3
    y_peaks = signal.find_peaks_cwt(y_sq,np.arange(0.1,1))
    print (y_peaks)
    if len(y_peaks) > 0:
        mean = np.mean(np.abs(y[y_peaks]))
        min = np.min(np.abs(y[y_peaks]))
        max = np.max(np.abs(y[y_peaks]))
    else:
        min = max = 0
        mean = y_mean
    mid = y_mean
    return mean, min, max, mid

def calc_mean_ampl_new(x,y,rp):
    initial = load_initial()
    U_X = float(initial[0])
    omega = float(initial[1])
    k = float(initial[2])
    print ("U_x = %s\nomega = %s\nk = %s" % (U_X,omega,k))
    x = x[rp:]
    y = y[rp:]
    i = a = b = X = 0
    N = len(y)-1
    for i in range(N):
        dx = float(x[i+1])-float(x[i])
        dy = 0.5*(float(y[i+1])+float(y[i]))
        X = float(x[N])-float(x[0])
        a  = a + 2.0/X*dx*dy*np.cos((U_X*k+omega)*float(x[i]))
        b  = b + 2.0/X*dx*dy*np.sin((U_X*k+omega)*float(x[i]))
        i += 1
    amplX = np.sqrt(a**2+b**2)
    phaseX = np.arctan(b/a)
    return amplX,phaseX

def find_nearest(array, value):
    idx = (np.abs(array-value)).argmin()
    return idx

def calc_phase(x,y,rp):
    initial = load_initial()
    U_X = float(initial[0])
    omega = float(initial[1])
    k = float(initial[2])
    print ("U_x = %s\nomega = %s\nk = %s" % (U_X,omega,k))
    x = x[rp:]
    y = y[rp:]
    i = a = b = X = 0
    N = len(y)-1
    for i in range(N):
        dx = float(x[i+1])-float(x[i])
        dy = 0.5*(float(y[i+1])+float(y[i]))
        X = float(x[N])-float(x[0])
        a  = a + 2.0/X*dx*dy*np.cos((U_X*k+omega)*float(x[i]))
        b  = b + 2.0/X*dx*dy*np.sin((U_X*k+omega)*float(x[i]))
        i += 1
    amplX = np.sqrt(a**2+b**2)
    phaseX = np.arctan(b/a)
    return amplX,phaseX

def fft_welch(x,y):
    dx = (x[-1]-x[0])/np.size(x)
    x_rad,y_fft=signal.welch(y,fs=1.0/dx,window='hanning',nperseg=128,scaling='density')
    x_rad=x_rad/2.0/math.pi
    return x_rad,y_fft

def fft_my(x,y):
    N = np.size(x)
    dx = (x[-1]-x[0])/N
    y = y*signal.hann(N)
    x_rad = np.linspace(0.0,1.0/(2.0*dx),N/2)
    y_fft = np.abs(fft(y,n=N))*2.0/N
    x_rad=x_rad*2.0*math.pi
    return x_rad,y_fft[:N/2]

def printValues(inputFile,i,k):
    data = np.loadtxt(inputFile, skiprows=4).transpose()
    peakFind = signal.find_peaks_cwt(data[i],np.arange(0,10))
    if k == 0:
        return data[k][peakFind]
    else:
        return data[i][peakFind]
