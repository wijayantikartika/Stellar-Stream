import numpy as np
import matplotlib.pyplot as plt
import time
'''
veer, vefi, vezi = np.loadtxt("semua.csv", delimiter=",",
                              usecols=(30, 31, 32), skiprows=1, unpack=True)
veerr, veteta, vefii = np.loadtxt("semua.csv", delimiter=",",
                             usecols=(33, 34, 35), skiprows=1, unpack=True)
'''
u_in, v_in, w_in = np.loadtxt("thick.csv", delimiter=",",
                              usecols=(27, 28, 29), skiprows=1, unpack=True)
E_in, lz_in, ltot_in, lperp_in, vaz_in, vdele_in = np.loadtxt(
        "thick.csv", delimiter=",",usecols=(38,39,40,41,42,43), 
        skiprows=1, unpack=True)

uuvv_in = np.power((np.power(u_in,2) + 2*np.power(v_in,2)),0.5)
'''
vaz1_in = vaz_in[vaz_in > 0]
vdele1_in = vdele_in[vaz_in > 0]
vaz2_in = vaz_in[vaz_in <= 0]
vdele2_in = vdele_in[vaz_in <= 0]

u_halo, v_halo, w_halo = np.loadtxt("halo1.csv", delimiter=",",
                              usecols=(27, 28, 29), skiprows=1, unpack=True)
E_halo, lz_halo, ltot_halo, lperp_halo, vaz_halo, vdele_halo = np.loadtxt(
        "halo1.csv", delimiter=",",usecols=(38,39,40,41,42,43), 
        skiprows=1, unpack=True)
uuvv_halo = np.power((np.power(u_halo,2) + 2*np.power(v_halo,2)),0.5)

u_thick, v_thick, w_thick = np.loadtxt("thick.csv", delimiter=",",
                              usecols=(27, 28, 29), skiprows=1, unpack=True)
E_thick, lz_thick, ltot_thick, lperp_thick, vaz_thick, vdele_thick = np.loadtxt(
        "thick.csv", delimiter=",",usecols=(38,39,40,41,42,43), 
        skiprows=1, unpack=True)
uuvv_thick = np.power((np.power(u_thick,2) + 2*np.power(v_thick,2)),0.5)
'''
#===========================================================================#
#                             Velocity Distribution                         #
#===========================================================================#
'''
#distribution of cylindrical velocities
fig, ax = plt.subplots(1, 3, dpi=200, figsize=(30,10))

ax[0].hist(veer, bins=100)
ax[0].set_xlabel(r'${V_R}$ (km/s)', fontsize=12)
ax[0].set_ylabel('N')

ax[1].hist(vefi, bins=100)
ax[1].set_xlabel(r'${V_\phi}$ (km/s)', fontsize=12)

ax[2].hist(vezi, bins=100)
#distribution of spherical velocities
ax[2].set_xlabel(r'${V_z}$ (km/s)', fontsize=12)

plt.savefig('vdist.png')
plt.show()

fig, ax = plt.subplots(1, 3, dpi=200, figsize=(30,10))
ax[0].hist(veerr, bins=100)
ax[0].set_xlabel(r'${V_r}$ (km/s)', fontsize=12)
ax[0].set_ylabel('N')

ax[1].hist(veteta, bins=100)
ax[1].set_xlabel(r'${V_\theta}$ (km/s)', fontsize=12)

ax[2].hist(vefii, bins=100)
ax[2].set_xlabel(r'${V_\phi}$ (km/s)', fontsize=12)

plt.savefig('vdistsph.png')
plt.show()
'''
#===========================================================================#
#                                Toomre Diagram                             #
#===========================================================================#
#u_thin, v_thin, w_thin = np.loadtxt("thin.csv", delimiter=",",
#                              usecols=(27, 28, 29), skiprows=1, unpack=True)
#u_thick, v_thick, w_thick = np.loadtxt("thick.csv", delimiter=",",
#                              usecols=(27, 28, 29), skiprows=1, unpack=True)
#u_halo, v_halo, w_halo = np.loadtxt("halo.csv", delimiter=",",
#                              usecols=(27, 28, 29), skiprows=1, unpack=True)
#u_halo1, v_halo1, w_halo1 = np.loadtxt("halo1.csv", delimiter=",",
#                              usecols=(27, 28, 29), skiprows=1, unpack=True)
'''
y_toomre_thin = np.power((np.power(u_thin,2) + np.power(w_thin,2)),0.5)
y_toomre_thick= np.power((np.power(u_thick,2) + np.power(w_thick,2)),0.5)
y_toomre_halo = np.power((np.power(u_halo,2) + np.power(w_halo,2)),0.5)
y_toomre_halo1= np.power((np.power(u_halo1,2) + np.power(w_halo1,2)),0.5)

fig, ax = plt.subplots()
for i in range(1,8):
    circle1 = plt.Circle((0,0), 100*i, color='black', linestyle='--', fill=False)
    ax.add_artist(circle1)
plt.scatter(v_thin, y_toomre_thin, c='blue', s=3, alpha=0.3, label='thin disk')
plt.scatter(v_thick, y_toomre_thick, c='orange', s=3, alpha=0.3, label='thick disk')
plt.scatter(v_halo, y_toomre_halo, c='green', s=3, alpha=0.3, label='halo (TD/H < 1)')
plt.ylim(0,600)
plt.xlabel('V (km/s)')
plt.ylabel(r'$\sqrt{U^2+W^2}$ (km/s)')
plt.legend()
plt.savefig('toomre1.png', dpi=300)
plt.show()

fig, ax = plt.subplots()
for i in range(1,8):
    circle1 = plt.Circle((0,0), 100*i, color='black', linestyle='--', fill=False)
    ax.add_artist(circle1)
plt.scatter(v_thin, y_toomre_thin, c='blue', s=3, alpha=0.3, label='thin disk')
plt.scatter(v_thick, y_toomre_thick, c='orange', s=3, alpha=0.3, label='thick disk')
plt.scatter(v_halo1, y_toomre_halo1, c='green', s=3, alpha=0.3, label='halo (TD/H < 0.01)')
plt.ylim(0,600)
plt.xlabel('V (km/s)')
plt.ylabel(r'$\sqrt{U^2+W^2}$ (km/s)')
plt.legend()
plt.savefig('toomre001.png', dpi=300)
plt.show()

fig, ax = plt.subplots(1, 3, dpi=200, figsize=(30,10), sharey=True)
fig.subplots_adjust(wspace=0.1)

ax[0].scatter(v_thin, y_toomre_thin, c='blue', s=3, alpha=0.3, label='thin disk')
ax[0].set_xlabel('V (km/s)', fontsize=12)
ax[0].set_ylabel(r'$\sqrt{U^2+W^2}$ (km/s)')
ax[0].title.set_text('Thin Disk')
for i in range(1,8):
    circle1 = plt.Circle((0,0), 100*i, color='black', linestyle='--', fill=False)
    ax[0].add_artist(circle1)
ax[0].set_xlim(-300,300)
ax[0].set_ylim(0,600)

ax[1].scatter(v_thick, y_toomre_thick, c='orange', s=3, alpha=0.3, label='thick disk')
ax[1].set_xlabel('V (km/s)', fontsize=12)
ax[1].title.set_text('Thick Disk')
for i in range(1,8):
    circle1 = plt.Circle((0,0), 100*i, color='black', linestyle='--', fill=False)
    ax[1].add_artist(circle1)
ax[1].set_xlim(-300,300)
ax[1].set_ylim(0,600)

ax[2].scatter(v_halo1, y_toomre_halo1, c='green', s=3, alpha=0.3, label='halo (TD/H < 0.01)')
ax[2].set_xlabel('V (km/s)', fontsize=12)
ax[2].title.set_text('Halo')
for i in range(1,8):
    circle1 = plt.Circle((0,0), 100*i, color='black', linestyle='--', fill=False)
    ax[2].add_artist(circle1)
ax[2].set_xlim(-300,300)
ax[2].set_ylim(0,600)

plt.savefig('toomredibagi.png')
plt.show()

x_toomre = v_in
y_toomre = np.power((np.power(u_in,2) + np.power(w_in,2)),0.5)

fig, ax = plt.subplots(dpi=200,figsize=(8,8))
for i in range(1,8):
    circle1 = plt.Circle((0,0),100*i,color='black',linestyle='--',fill=False)
    ax.add_artist(circle1)
plt.scatter(x_toomre,y_toomre,s=3,alpha=0.3)
plt.ylim(0,600)
plt.xlabel('V (km/s)')
plt.ylabel(r'$\sqrt{U^2+W^2}$ (km/s)')
plt.savefig('halo_li.png')
plt.show()'''
#===========================================================================#
#                                Scatter Plots                              #
#===========================================================================#

plt.figure(dpi=200, figsize=(12,18))

plt.subplot(321)
plt.scatter(u_in, v_in,s=3)
plt.xlabel('U (km/s)',fontsize=12)
plt.ylabel('V (km/s)',fontsize=12)

plt.subplot(322)
plt.scatter(v_in,uuvv_in,s=3)
plt.xlabel('V (km/s)',fontsize=12)
plt.ylabel(r'$\sqrt{U^2+2V^2}$ (km/s)',fontsize=12)

plt.subplot(323)
plt.scatter(vaz_in,vdele_in,s=3)
plt.xlabel(r'$V_{az}$ (km/s)',fontsize=12)
plt.ylabel(r'$V_{\Delta E}$ (km/s)',fontsize=12)

plt.subplot(324)
plt.scatter(lz_in,lperp_in,s=3)
plt.xlabel(r'$L_{z}$ (kpc km/s)',fontsize=12)
plt.ylabel(r'$L_{\perp}$ (kpc km/s)',fontsize=12)

Ein = E_in/ 1000
plt.subplot(325)
plt.scatter(lz_in,Ein,s=3)
plt.xlabel(r'$L_{z}$ (kpc km/s)',fontsize=12)
plt.ylabel(r'$E ({\times}10^3) (km/s)^{2}$',fontsize=12)

plt.savefig('thickscattersemwa.png')
plt.show()
'''
plt.scatter(v_halo,uuvv_halo,s=3, alpha=0.3)
#plt.xlabel('U (km/s)')
#plt.ylabel('V (km/s)')
plt.xlabel('V (km/s)')
plt.ylabel(r'$\sqrt{U^2+2V^2}$ (km/s)')
plt.title('Halo', fontsize=12)
plt.savefig('halo-uuvv.png')
plt.show()

plt.scatter(v_halo,uuvv_halo,s=3, alpha=0.3)
plt.xlabel('V (km/s)')
plt.ylabel(r'$\sqrt{U^2+2V^2}$ (km/s)')
plt.title('Halo', fontsize=12)
plt.savefig('halo-uuvv.png')
plt.show()
'''
#===========================================================================#
#                                   Wavelet                                 #
#===========================================================================#

def mexico (rr, scale):
    return (2 - rr**2 / (scale**2)) * np.exp(-rr**2 / (2*scale**2))

def koefisien (xstart, xend, ystart, yend, dx, dy, x, y, a):
    xgrid = np.arange(xstart, xend, dx)
    ygrid = np.arange(ystart, yend, dy)
    nxgrid= len(xgrid)
    nygrid= len(ygrid)
    coef = np.zeros((nxgrid, nygrid))
    
    for i in range (len(x)-1):
        for p in range (nxgrid):
            for q in range (nygrid):
                xx = xgrid[p] - x[i]
                yy = ygrid[q] - y[i]
                r = np.sqrt(xx*xx + yy*yy)
                coef[q][p] += mexico(r, a)
    coef[coef<0]=0
    return coef
'''
u_halo_mc = [[] for i in range (5)]
v_halo_mc = [[] for i in range (5)]
w_halo_mc = [[] for i in range (5)]
uuvv_halo_mc = [[] for i in range (5)]

u1_halo_mc = [[] for i in range (5)]
v1_halo_mc = [[] for i in range (5)]
w1_halo_mc = [[] for i in range (5)]

u2_halo_mc = [[] for i in range (5)]
v2_halo_mc = [[] for i in range (5)]
w2_halo_mc = [[] for i in range (5)]

vaz1_halo_mc = [[] for i in range (5)]
vaz2_halo_mc = [[] for i in range (5)]
vdele1_halo_mc = [[] for i in range (5)]
vdele2_halo_mc = [[] for i in range (5)]
#u_thick_mc = [[] for i in range (5)]
#v_thick_mc = [[] for i in range (5)]
for i in range (5):
    u_halo_mc[i] = np.random.normal(0, 160, len(u_in))
    v_halo_mc[i] = np.random.normal(-220, 90, len(v_in))
    w_halo_mc[i] = np.random.normal(0, 90, len(w_in))
    uuvv_halo_mc[i] = np.power((np.power(u_halo_mc[i],2) + 2*np.power(v_halo_mc[i],2)),0.5)
    
    u1_halo_mc[i] = np.random.normal(0, 160, len(vaz1_in))
    v1_halo_mc[i] = np.random.normal(-220, 90, len(vaz1_in))
    w1_halo_mc[i] = np.random.normal(0, 90, len(vaz1_in))
    
    u2_halo_mc[i] = np.random.normal(0, 160, len(vaz2_in))
    v2_halo_mc[i] = np.random.normal(-220, 90, len(vaz2_in))
    w2_halo_mc[i] = np.random.normal(0, 90, len(vaz2_in))
    
    vaz1_halo_mc[i] = np.power((np.power((v1_halo_mc[i] + 232.8),2) + np.power(w1_halo_mc[i],2)),0.5)
    vaz2_halo_mc[i] = -np.power((np.power((v2_halo_mc[i] + 232.8),2) + np.power(w2_halo_mc[i],2)),0.5)
    vdele1_halo_mc[i] = np.power((np.power(u1_halo_mc[i],2) + 2*np.power((vaz1_halo_mc[i] - 232.8),2)),0.5)
    vdele2_halo_mc[i] = np.power((np.power(u2_halo_mc[i],2) + 2*np.power((vaz2_halo_mc[i] + 232.8),2)),0.5)

#    u_thick_mc[i] = np.random.normal(0, 67, len(u_thick))
#    v_thick_mc[i] = np.random.normal(-46, 38, len(v_thick))

def koreksi (xstart, xend, ystart, yend, dx, dy, x, y, x_mc, y_mc, a):
    m = koefisien(xstart, xend, ystart, yend, dx, dy, x, y, a)
    n1= koefisien(xstart, xend, ystart, yend, dx, dy, x_mc[0], y_mc[0], a)
    n2= koefisien(xstart, xend, ystart, yend, dx, dy, x_mc[1], y_mc[1], a)
    n3= koefisien(xstart, xend, ystart, yend, dx, dy, x_mc[2], y_mc[2], a)
    n4= koefisien(xstart, xend, ystart, yend, dx, dy, x_mc[3], y_mc[3], a)
    n5= koefisien(xstart, xend, ystart, yend, dx, dy, x_mc[4], y_mc[4], a)
    koefisien1 = m - n1
    koefisien2 = m - n2
    koefisien3 = m - n3
    koefisien4 = m - n4
    koefisien5 = m - n5
    ratarata = (np.array(koefisien1) + np.array(koefisien2) +
                np.array(koefisien3) + np.array(koefisien4) + 
                np.array(koefisien5))/5
    ratarata[ratarata<0] = 0
    return ratarata

start_time = time.time()

xgrid = np.arange(0,300,10)
ygrid = np.arange(200,500,10)         
X, Y = np.meshgrid(xgrid, ygrid)

xgrid1 = np.arange(-300,0, 10)
ygrid1 = np.arange(200,500, 10)        
X1, Y1 = np.meshgrid(xgrid1, ygrid1)

#z2 = koefisien(-400,200,0,600,5,5,v_in,uuvv_in,2)
#z4 = koefisien(-400,200,0,600,5,5,v_in,uuvv_in,4)
#z6 = koefisien(-400,200,0,600,5,5,v_in,uuvv_in,6)
z8 = koefisien(0,300,200,500,10,10,vaz1_in,vdele1_in,10)
#z8c = koreksi(0,600,0,600,10,10,vaz1_in,vdele1_in,vaz1_halo_mc,vdele1_halo_mc,10)
z9 = koefisien(-300,0,200,500,10,10,vaz2_in,vdele2_in,10)
#z9c = koreksi(-600,0,0,600,10,10,vaz2_in,vdele2_in,vaz2_halo_mc,vdele2_halo_mc,10)
#z8cc = koreksi(-300,300,0,600,10,10,vaz_in,vdele_in,v_halo_mc,uuvv_halo_mc,10)
#z10= koefisien(-400,200,0,600,5,5,v_in,uuvv_in,10)
#z12= koefisien(-400,200,0,600,5,5,v_in,uuvv_in,12)

plt.figure(0)
plt.contourf(X,Y,koefisien (-200,200,0,400,10,10,u_halo1,v_halo1,8),cmap='Reds')
plt.title('Halo, a = 8 km/s')
plt.xlabel('U (km/s)')
plt.ylabel('V (km/s)')
plt.colorbar()
plt.savefig('Halo-10-8.png')

fig, ax = plt.subplots(3, 2, dpi=300, figsize=(20,30), sharey=True)
fig.subplots_adjust(wspace=0.1)
fig.suptitle('Halo(Li, 2017)')

ax[0].contourf(X,Y,z2,cmap='Reds')
ax[0].set_ylabel(r'$\sqrt{U^2+2V^2}$ (km/s)')
ax[0].title.set_text('a = 2 km/s')

ax[1].contourf(X,Y,z4,cmap='Reds')
ax[1].title.set_text('a = 4 km/s')

ax[2].contourf(X,Y,z6,cmap='Reds')
ax[2].set_ylabel(r'$\sqrt{U^2+2V^2}$ (km/s)')
ax[2].title.set_text('a = 6 km/s')

ax[3].contourf(X,Y,z8,cmap='Reds')
ax[3].title.set_text('a = 8 km/s')

ax[4].contourf(X,Y,z10,cmap='Reds')
ax[4].set_xlabel('V (km/s)')
ax[4].set_ylabel(r'$\sqrt{U^2+2V^2}$ (km/s)')
ax[4].title.set_text('a = 10 km/s')

ax[5].contourf(X,Y,z12,cmap='Reds')
ax[5].set_xlabel('V (km/s)')
ax[5].set_ylabel(r'$\sqrt{U^2+2V^2}$ (km/s)')
ax[5].title.set_text('a = 12 km/s')

plt.savefig('cekskalahaloin.png')
plt.show()

plt.figure(dpi=200, figsize=(12,18))

plt.subplot(321)
plt.contourf(X,Y,z2,cmap='Reds')
plt.xlabel('V (km/s)')
plt.ylabel(r'$\sqrt{U^2+2V^2}$ (km/s)')
#plt.xlabel('U (km/s)')
#plt.ylabel('V (km/s)')
plt.title('a = 2 km/s', fontsize=12)
plt.colorbar()

plt.subplot(322)
plt.contourf(X,Y,z4,cmap='Reds')
plt.xlabel('V (km/s)')
plt.ylabel(r'$\sqrt{U^2+2V^2}$ (km/s)')
plt.title('a = 4 km/s', fontsize=12)
plt.colorbar()

plt.subplot(323)
plt.contourf(X,Y,z6,cmap='Reds')
plt.xlabel('V (km/s)')
plt.ylabel(r'$\sqrt{U^2+2V^2}$ (km/s)')
plt.title('a = 6 km/s', fontsize=12)
plt.colorbar()

plt.subplot(324)
plt.contourf(X,Y,z8,cmap='Reds')
plt.xlabel('V (km/s)')
plt.ylabel(r'$\sqrt{U^2+2V^2}$ (km/s)')
plt.title('a = 8 km/s', fontsize=12)
plt.colorbar()

plt.subplot(325)
plt.contourf(X,Y,z10,cmap='Reds')
plt.xlabel('V (km/s)')
plt.ylabel(r'$\sqrt{U^2+2V^2}$ (km/s)')
plt.title('a = 10 km/s', fontsize=12)
plt.colorbar()

plt.subplot(326)
plt.contourf(X,Y,z12,cmap='Reds')
plt.xlabel('V (km/s)')
plt.ylabel(r'$\sqrt{U^2+2V^2}$ (km/s)')
plt.title('a = 12 km/s', fontsize=12)
plt.colorbar()

plt.savefig('cek4haloin.png')
plt.show()

plt.figure(dpi=200)#, figsize=(12,12))

plt.subplot(121)
plt.contourf(X1,Y1,z9,cmap='Reds')
plt.xlabel(r'$V_{az}$ (km/s)', fontsize=14)
plt.ylabel(r'$V_{\Delta E}$ (km/s)', fontsize=14)
plt.colorbar()

plt.subplot(122)
plt.contourf(X,Y,z8,cmap='Reds')
plt.xlabel(r'$V_{az}$ (km/s)', fontsize=14)
plt.ylabel(r'$V_{\Delta E}$ (km/s)', fontsize=14)
plt.colorbar()
plt.show()

plt.subplot(223)
plt.contourf(X1,Y1,z9c,cmap='Reds')
plt.xlabel(r'$V_{az}$ (km/s)', fontsize=14)
plt.ylabel(r'$V_{\Delta E}$ (km/s)', fontsize=14)
plt.colorbar()

plt.subplot(224)
plt.contourf(X,Y,z8c,cmap='Reds')
plt.xlabel(r'$V_{az}$ (km/s)', fontsize=14)
plt.ylabel(r'$V_{\Delta E}$ (km/s)', fontsize=14)
plt.colorbar()

plt.savefig('kelupaan.png')
plt.show()

plt.figure(dpi=200, figsize=(12,12))
plt.contourf(X,Y,z8,cmap='Reds')
#plt.xlabel(r'$V_{az}$ (km/s)', fontsize=14)
#plt.ylabel(r'$V_{\Delta E}$ (km/s)', fontsize=14)
plt.colorbar()
plt.show()
'''
elapsed_time=(time.time()-start_time)/60
print(elapsed_time,'menit')