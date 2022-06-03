import numpy as np
import scipy.special as sps
import matplotlib.pyplot as plt
import matplotlib
matplotlib.style.use('classic')
#plt.rcParams['text.usetex'] = True

# Particles size disribution of NMC particles
#Gamma distribution with alpha=4.4152, beta=1.9944, gives particle diameter in m
if __name__ == '__main__':
    shape, scale = 4.4152, 1.9944  # mean=4, std=2*sqrt(2)
    s = np.random.gamma(shape, scale, 50000)
#    s = s[ (s>=2*2) & (s<=12*2)]
    print(len(s))
    count, bins, ignored = plt.hist(s, 100, density=True)

#    xarr = np.array(np.linspace(0,50,100))
#    print(type(xarr))
#    y = xarr ** (shape - 1) * (np.exp(-xarr / scale) / (sps.gamma(shape) * scale ** shape))

    y = bins ** (shape - 1) * (np.exp(-bins / scale) /(sps.gamma(shape) * scale ** shape))
    plt.plot(bins, y, linewidth=2, color='r')
    plt.xlabel('Particle diameter [µm]', fontsize=16, weight='bold')
#    plt.ylabel('Number of particles [-]')
    plt.ylabel('Probability density function [-]', fontsize=16, weight='bold')
    plt.xlim([4, 24])
   # plt.legend(loc='best')
    plt.text(20,.1,s=r'$k=4.42 $'+'\n'+r'$\theta=1.99$' , fontsize=16)



    PDF_fig,ax = plt.subplots()
    ax.plot(bins, y, linewidth=3, color='b')
    ax.text(19,.09,s='Gamma\n'+r'$k=4.42 $'+'\n'+r'$\theta=1.99$' , fontsize=16)
    #ax.xlim([4, 24])
    ax.set_xlim(4, 24)
    ax.set_xlabel('Particle diameter [µm]', fontsize=16, weight='bold')
    #    plt.ylabel('Number of particles [-]')
    ax.set_ylabel('Probability density function [-]', fontsize=16, weight='bold')
    ax.tick_params(axis='both', which='major', labelsize=16)
    plt.show()