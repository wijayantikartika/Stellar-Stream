from astropy.io import fits
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
#=======================================================================#
#                              Data                                     #
#=======================================================================#

hdul = fits.open('gaia_lamost_1arcsec_galcoor_uvw.fits')

hdr = hdul[1].header
data = hdul[1].data

#untuk menghitung kecepatan
RA = data.field(34)
DEC = data.field(35)
feh = data.field(40)
rv = data.field(42)
p = data.field(54)
pmRA = data.field(56)
pmDEC = data.field(58)
rvgaia = data.field(71)
bujur = data.field(80)
lintang = data.field(81)
U = data.field(83)
V = data.field(84)
W = data.field(85)

ra = RA.byteswap().newbyteorder()
dec = DEC.byteswap().newbyteorder()
feperh = feh.byteswap().newbyteorder()
RVel = rv.byteswap().newbyteorder()
par = p.byteswap().newbyteorder()
pmra = pmRA.byteswap().newbyteorder()
pmdec = pmDEC.byteswap().newbyteorder()
rv_gaia = rvgaia.byteswap().newbyteorder()
bujur = bujur.byteswap().newbyteorder()
lintang = lintang.byteswap().newbyteorder()
U = U.byteswap().newbyteorder()
V = V.byteswap().newbyteorder()
W = W.byteswap().newbyteorder()

#menghitung galat dan seleksi data
RV_err = data.field(43)
rverr = data.field(72)
errp = data.field(55)
pmra_error = data.field(57)
pmdec_error = data.field(59)

RV_err = RV_err.byteswap().newbyteorder()
rverr = rverr.byteswap().newbyteorder()
errp = errp.byteswap().newbyteorder()
pmra_error = pmra_error.byteswap().newbyteorder()
pmdec_error = pmdec_error.byteswap().newbyteorder()

#seleksi fotometri
gmeanflux = data.field(61)
gmeanerr = data.field(62)
gmeanmag = data.field(63)
bpmeanflux = data.field(64)
bpmeanerr = data.field(65)
bpmeanmag = data.field(66)
rpmeanflux = data.field(67)
rpmeanerr = data.field(68)
rpmeanmag = data.field(69)
bp_rp = data.field(70)

gflux = gmeanflux.byteswap().newbyteorder()
gerr = gmeanerr.byteswap().newbyteorder()
gmag = gmeanmag.byteswap().newbyteorder()
bpflux = bpmeanflux.byteswap().newbyteorder()
bperr = bpmeanerr.byteswap().newbyteorder()
bpmag = bpmeanmag.byteswap().newbyteorder()
rpflux = rpmeanflux.byteswap().newbyteorder()
rperr = rpmeanerr.byteswap().newbyteorder()
rpmag = rpmeanmag.byteswap().newbyteorder()
bprp = bp_rp.byteswap().newbyteorder()

#========================================================================#
#                               Data Selection                           #
#========================================================================#

data = {'RA': ra,'DEC':dec, 'l':bujur, 'b':lintang,'[Fe/H]': feperh,
        'RV (lamost)': RVel, 'RV_err': RV_err,'RV (gaia)': rv_gaia, 
        'RVerr gaia': rverr, 'p': par, 'p_err': errp,
        'pmRA': pmra, 'pmDEC': pmdec,'g_flux': gflux, 'g_err': gerr,
        'g_mag': gmag, 'bp_flux': bpflux, 'bp_err': bperr, 'bp_mag': bpmag, 
        'rp_flux': rpflux, 'rp_err': rperr, 'rp_mag': rpmag, 'bp_rp': bprp, 
        'U (topcat)':U, 'V (topcat)':V, 'W (topcat)': W,
        'pmRA_err': pmra_error, 'pmDEC_err':pmdec_error}

df = pd.DataFrame(data)

df = df[np.isfinite(df['[Fe/H]'])]
df = df[np.isfinite(df['RV (lamost)'])]
df = df[np.isfinite(df['p'])]
df = df[np.isfinite(df['pmRA'])]
df = df[np.isfinite(df['pmDEC'])]
df = df[np.isfinite(df['RV (gaia)'])]

df = df[df['p'] > 0]
df = df[df['p_err'] / df['p'] < 0.10]
df = df[df['RV_err'] <= 10]

excessfactor = (df['bp_flux'] + df['rp_flux']) / df['g_flux']
df['excessfactor'] = excessfactor

df = df[df['g_flux'] / df['g_err'] > 50]
df = df[df['rp_flux'] / df['rp_err'] > 20]
df = df[df['bp_flux'] / df['bp_err'] > 20]
df = df[df['excessfactor'] > (1 + 0.015 * (df['bp_rp']) ** 2)]
df = df[df['excessfactor'] < (1.3 + 0.06 * (df['bp_rp']) ** 2)]

#offset
df['p'] = df['p'] - 0.029
df['RV (lamost)'] = df['RV (lamost)'] - 4.7
#==========================================================================#
#                                    (U, V, W)                             #
#==========================================================================#

RA_NGP=np.deg2rad(192.85948)
DEC_NGP=np.deg2rad(27.12825)
t_NGP=np.deg2rad(122.93192)

t1 = [[np.cos(t_NGP),np.sin(t_NGP),0] , [np.sin(t_NGP),-np.cos(t_NGP),0] , 
     [0,0,1]]
t_1 = np.reshape(t1,(-1,3))
t2 = [[-np.sin(DEC_NGP),0,np.cos(DEC_NGP)] , [0,-1,0] , [np.cos(DEC_NGP),0,
       np.sin(DEC_NGP)]]
t_2 = np.reshape(t2,(-1,3))
t3 = [[np.cos(RA_NGP),np.sin(RA_NGP),0] , [np.sin(RA_NGP),-np.cos(RA_NGP),0] , 
       [0,0,1]]
t_3 = np.reshape(t3,(-1,3))
t12 = np.matmul(t_1,t_2)
T = np.matmul(t12,t_3)

alpa = np.deg2rad(df['RA'])
delta = np.deg2rad(df['DEC'])

a00 = np.cos(alpa) * np.cos(delta)
a01 = -np.sin(alpa)
a02 = -np.cos(alpa) * np.sin(delta)

a10 = np.sin(alpa) * np.cos(delta)
a11 = np.cos(alpa)
a12 = -np.sin(alpa) * np.sin(delta)

a20 = np.sin(delta)
a21 = 0
a22 = np.cos(delta)

b11 = T[0][0]*a00 + T[0][1]*a10 + T[0][2]*a20
b12 = T[0][0]*a01 + T[0][1]*a11 + T[0][2]*a21
b13 = T[0][0]*a02 + T[0][1]*a12 + T[0][2]*a22

b21 = T[1][0]*a00 + T[1][1]*a10 + T[1][2]*a20
b22 = T[1][0]*a01 + T[1][1]*a11 + T[1][2]*a21
b23 = T[1][0]*a02 + T[1][1]*a12 + T[1][2]*a22

b31 = T[2][0]*a00 + T[2][1]*a10 + T[2][2]*a20
b32 = T[2][0]*a01 + T[2][1]*a11 + T[2][2]*a21
b33 = T[2][0]*a02 + T[2][1]*a12 + T[2][2]*a22

k = 4.74057
m1 = df['RV (lamost)']
m2 = k * df['pmRA'] / df['p']
m3 = k * df['pmDEC'] / df['p']

u = b11*m1 + b12*m2 + b13*m3 + 11.1
v = b21*m1 + b22*m2 + b23*m3 + 12.24 #+ 232.8
w = b31*m1 + b32*m2 + b33*m3 + 7.25

df['U'] = u
df['V'] = v
df['W'] = w

k=4.74057
c11=b11**2
c12=b12**2
c13=b13**2

c21=b21**2
c22=b22**2
c23=b23**2

c31=b31**2
c32=b32**2
c33=b33**2

n1=df['RV_err']**2
n2=((k/df['p'])**2)*((df['pmRA_err']**2)+(df['pmRA']*df['p_err']/df['p'])**2)
n3=((k/df['p'])**2)*((df['pmDEC_err']**2)+(df['pmDEC']*df['p_err']/df['p'])**2)

konst2=(2*df['pmRA']*df['pmDEC']*(k**2)*(df['p_err']**2))/(df['p']**4)

dispersi_u=np.sqrt((c11*n1)+(c12*n2)+(c13*n3)+(konst2*b12*b13))
dispersi_v=np.sqrt((c21*n1)+(c22*n2)+(c23*n3)+(konst2*b22*b23))
dispersi_w=np.sqrt((c31*n1)+(c32*n2)+(c33*n3)+(konst2*b32*b33))
#=========================================================================#
#                       Galactocenrtic Coordinates                        #
#=========================================================================#

#posisi
r0 = 8.2 #kpc
z0 = 0.025 #kpc

ex = r0 - (1/df['p'])*np.cos(np.deg2rad(df['b']))*np.cos(np.deg2rad(df['l']))
ye = -(1/df['p'])*np.cos(np.deg2rad(df['b']))*np.sin(np.deg2rad(df['l']))
zi = (1/df['p'])*np.sin(np.deg2rad(df['b'])) + z0

#==========================================================================#
#                     Cylindrical and Spherical Velocities                 #
#==========================================================================#

#cylindrical velocities
fi = np.arctan(ye/ex)
ER = np.sqrt(ex**2 + ye**2)
er = np.sqrt(ex**2 + ye**2 + zi**2)
theta = np.arccos(zi/er)

veer = -df['U']*np.cos(fi) -(df['V']+232.8)*np.sin(fi)
vefi = df['U']*np.sin(fi) -(df['V']+232.8)*np.cos(fi)
vezi = -df['W']

veerr = -df['U']*np.sin(theta)*np.cos(fi)-(df['V']+232.8)*np.sin(theta)*np.sin(
        fi)-df['W']*np.cos(theta)
veteta = -df['U']*np.cos(theta)*np.cos(fi)-(df['V']+232.8)*np.cos(theta
            )*np.sin(fi)+df['W']*np.sin(theta)
vefii = df['U']*np.sin(fi) - (df['V']+232.8)*np.cos(fi)

df['V_R'] = veer
df['V_phi'] = vefi
df['V_z'] = vezi

df['V_r'] = veerr
df['V_theta'] = veteta
df['V_phi (spherical)'] = vefii

#===========================================================================#
#         Inclination, Eccentricity, Energy, and Angle Momentum             #
#===========================================================================#
#inklinasi
inc = np.arctan((df['V'] + 232.8)/df['W'])
i = np.rad2deg(inc)
i[i < 0] = i + 360
df['i (deg)'] = i

#eksentrisitas
ecc = np.sqrt((np.power(df['U'],2) + 2*np.power((df['V']+232.8),2))/(
        2*232.8**2))
df['e'] = ecc

#potensial
nuhalo = 173.2 #km/s
ded = 12 #kpc
p_halo = (nuhalo**2) *np.log(1 + (ER/ded)**2 + (zi/ded)**2)

G = 4.3*(10**-6) #kpc Mo^-1 (km/s)^2
Mdisk = 6.3*(10**10) #Mo
ad = 6.5 #kpc
bd = 0.26 #kpc
p_disk = -(G*Mdisk)/ (np.sqrt(ER**2 + (ad + np.sqrt(zi**2 + bd**2))**2))

Mbulge = 2.1*(10**10) #Mo
cb = 0.7 # kpc
p_bulge = -(G*Mbulge)/(er + cb)

p_tot = 0.5*(df['U']**2 + (df['V']+232.8)**2 + df['W']**2
             ) + p_halo + p_disk + p_bulge
df['E (km/s)^2'] = p_tot

#momentum sudut
lx = ye*df['W'] - zi*(df['V']+232.8)
ly = zi*df['U'] - ex*df['W']
lz = ex*(df['V']+232.8) - ye*df['U']

ltot = np.sqrt(lx**2+ly**2+lz**2)
lperp = np.sqrt(lx**2+ly**2)

df['L_z (kpc km/s)'] = lz
df['L_tot (kpc km/s)'] = ltot
df['L_perp (kpc km/s)'] = lperp

vaz = np.sqrt((df['V'] + 232.8)**2 + df['W']**2)
vdele = np.sqrt(df['U']**2 + (2*(vaz - 232.8))**2)
#vaz[i > 180] = -np.sqrt((df['V'] + 232.8)**2 + df['W']**2)
#vdele[i > 180] = np.sqrt(df['U']**2 + (2*(vaz + 232.8))**2)

df['V_az'] = vaz
df['V_delE'] = vdele
#==========================================================================#
#                   Seleksi Komponen Galaksi: Bensby (2014)                #
#==========================================================================#

#urutan list (thin disk (D), thick disk (TD), halo (H))
sigma_u = [35,67,160]
sigma_v = [20,38,90]
sigma_w = [16,35,90]
u_assym = [0,0,0]
v_assym = [-15,-46,-220]
fraksi = [0.85,0.09,0.0015]

kd = 1 / (((2*np.pi) **1.5) *sigma_u[0] *sigma_v[0] *sigma_w[0])
ktd = 1 / (((2*np.pi) **1.5) *sigma_u[1] *sigma_v[1] *sigma_w[1])
kh = 1 / (((2*np.pi) **1.5) *sigma_u[2] *sigma_v[2] *sigma_w[2])

suku1d = -(((df['U']-u_assym[0])**2)/(2*sigma_u[0]**2))
suku2d = -(((df['V']-v_assym[0])**2)/(2*sigma_v[0]**2))
suku3d = -(((df['W'])**2)/(2*sigma_w[0]**2))

suku1td = -(((df['U']-u_assym[1])**2)/(2*sigma_u[1]**2))
suku2td = -(((df['V']-v_assym[1])**2)/(2*sigma_v[1]**2))
suku3td = -(((df['W'])**2)/(2*sigma_w[1]**2))

suku1h = -(((df['U']-u_assym[2])**2)/(2*sigma_u[2]**2))
suku2h = -(((df['V']-v_assym[2])**2)/(2*sigma_v[2]**2))
suku3h = -(((df['W'])**2)/(2*sigma_w[2]**2))

fd = kd * np.exp(suku1d + suku2d + suku3d)
ftd= ktd* np.exp(suku1td + suku2td + suku3td)
fh = kd * np.exp(suku1h + suku2h + suku3h)

#probabilitas relatif
tdd = (fraksi[1] / fraksi[0]) * (ftd/fd)     #TD/D=(XTD/XD)*(fTD/fD)
tdh = (fraksi[1] / fraksi[2]) * (ftd/fh)     #TD/H=(XTD/XH)*(fTD/fH)

mask_thin = (tdd < 0.5)
mask_thick = (tdd > 2)
mask_halo = (tdh < 1)
mask_halo1 = (tdh <0.01)

#===========================================================================#
#                     Halo Selection: Hefan Li (2017)                      #
#==========================================================================#

#seleksi halo lain
#halo
ka_halo = 1/(((2*np.pi)**2)*141*75*85*0.36)
dalam1_halo = -(veerr**2)/(2*141**2)
dalam2_halo = -(veteta**2/(2*75**2))
dalam3_halo = -((vefii)**2)/(2*85**2)
dalam4_halo = -((df['[Fe/H]']+1.43)**2)/(2*0.36**2)
H = ka_halo*np.exp(dalam1_halo+dalam2_halo+dalam3_halo+dalam4_halo)

#thin
ka_thin = 1/(((2*np.pi)**2)*35*22*17*(0.22))
dalam1_thin = -(veerr**2)/(2*35**2)
dalam2_thin = -(veteta**2/(2*22**2))
dalam3_thin = -((vefii+221)**2)/(2*17**2)
dalam4_thin = -((df['[Fe/H]']+0.10)**2)/(2*0.22**2)
D = ka_thin*np.exp(dalam1_thin+dalam2_thin+dalam3_thin+dalam4_thin)

#thick
ka_thick = 1/(((2*np.pi)**2)*63*45*43*(0.24))
dalam1_thick = -(veerr**2)/(2*63**2)
dalam2_thick = -(veteta**2/(2*45**2))
dalam3_thick = -((vefii+185)**2)/(2*43**2)
dalam4_thick = -((df['[Fe/H]']+0.47)**2)/(2*0.24**2)
TD = ka_thick*np.exp(dalam1_thick+dalam2_thick+dalam3_thick+dalam4_thick)

sbx = np.log(H/TD)
sby = np.log(H/D)

sbx = sbx[vefi < -40]
sby = sby[vefi < -40]
sbx = sbx[df['[Fe/H]'] > -1.4]
sby = sby[df['[Fe/H]'] > -1.4]
sbx = sbx[H > D]
sby = sby[H > D]
sbx = sbx[H > TD]
sby = sby[H > TD]

def func(x, m, c):
    return m * x + c
def ln(x, a, c):
    return (a/x) + c
popt, pcov = curve_fit(func, sbx, sby)
#print(popt)

sbxx = sbx[sby > (90/sbx) + popt[1]]
sbyy = sby[sby > (90/sbx) + popt[1]]

haha = np.linspace(0.01, max(sbx), 10000)
'''
plt.figure(dpi=300)
plt.scatter(sbx,sby,s=1)
plt.scatter(sbxx,sbyy,s=1)
plt.xlabel('ln(H/TD)')
plt.ylabel('ln(H/D)')
plt.xlim(0,30)
plt.ylim(0,200)
plt.plot(sbx, func(sbx, *popt), 'k-',linestyle=':')
plt.plot(haha, ln(haha, 90, popt[1]), 'k-', linewidth=1)
plt.savefig('pilihhalo.png')
plt.show()'''
#===========================================================================#
#                                 Saving . . . . . .                        #
#==========================================================================#
df['X'] = ex
df['Y'] = ye
df['Z'] = zi

df['P_tot']= p_halo + p_disk + p_bulge

df['U_err'] = dispersi_u
df['V_err'] = dispersi_v
df['W_err'] = dispersi_w

#export_thin = df[mask_thin].to_csv (r'D:\kuliah\TA\data\thinbaru.csv', 
#                index = None, header=True)
#export_thick = df[(tdd > 2)&(df['[Fe/H]'] > -1.6)&(df['[Fe/H]'] <-0.2)].to_csv ('thickbarunyobalagi.csv', 
#                 index = None, header=True)
#export_halo = df[mask_halo].to_csv ('halobarulagi.csv', 
#                index = None, header=True)
export_halo1 = df[(tdh< 0.01) & (df['[Fe/H]'] < 1)].to_csv ('halo1barulagi.csv', 
                 index = None, header=True)

#df2 = df[(vefi < -40) & (df['[Fe/H]'] > -1.4) & (H > D) & (H > TD) & (sby > 
#         (90/sbx) + popt[1])]

#export_haloin = df2.to_csv (r'D:\kuliah\TA\data\halohalobandungbaru.csv', 
#                            index = None, header=True)
#export_semua = df.to_csv (r'D:\kuliah\TA\data\semua.csv', index = None, 
#                          header=True)
'''
fig, ax = plt.subplots(1, 3, dpi=200, figsize=(30,10))

ax[0].hist(df['U'], bins=100)
ax[0].set_xlabel(r'$U$ (km/s)', fontsize=12)
ax[0].set_ylabel('N')

ax[1].hist(df['V'], bins=100)
ax[1].set_xlabel(r'$V$ (km/s)', fontsize=12)

ax[2].hist(df['W'], bins=100)
#distribution of spherical velocities
ax[2].set_xlabel(r'$W$ (km/s)', fontsize=12)

plt.savefig('kec.png')
plt.show()'''