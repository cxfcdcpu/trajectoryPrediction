import numpy as np
import pandas as pd
import geopandas as gpd
import json
from shapely.geometry import LineString, Point
from multiprocessing import  Pool

def load_data(fname, nrows):
    if nrows != -1:
      df = pd.read_csv(fname, nrows=nrows)
    else:
      df = pd.read_csv(fname)
    df['traj'] = json.loads('[' + df.POLYLINE.str.cat(sep=',') + ']')
    df = df[df.traj.str.len() > 1].copy()
    df['lines'] = gpd.GeoSeries(df.traj.apply(LineString))
    return gpd.GeoDataFrame(df, geometry='lines')

df = load_data('./train.csv',nrows = 1000)
#df.head()

import matplotlib
matplotlib.rcParams['figure.figsize'] = [15,8]
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style="darkgrid")




from sklearn.neighbors import DistanceMetric
metric = DistanceMetric.get_metric('haversine')
R = 6371 # radius of earth in km
dt = 15/3600 # coordinates are reported in 15 second intervals

def haversine(x):
    return metric.pairwise(np.radians(x)[:,::-1])

def velocity_graph(ax, coords):
    n = len(coords)
    dist = haversine(coords)
    interval = dt*np.abs(np.arange(n)[:,None] - np.arange(n)).clip(1)
    vel = R*dist/interval
    sns.heatmap(vel, ax=ax, square=True, robust=True)
    ax.set_title('Velocity (km/h)')
    ax.set_xlabel('GPS Sample Index')
    ax.set_ylabel('GPS Sample Index')


def offset_distances(coords, offset):
    dist = R*haversine(coords)
    dists = np.diag(dist[offset:])
    return dists

def plot_dists(ax, di):
    dists = np.hstack(df.traj.apply(lambda x: offset_distances(x, di)).values)
    sns.distplot(dists, ax=ax, kde=False, bins=300)
    ax.set_title('Time Offset: {} min'.format(di*dt*60))
    ax.set_xlabel('Distance (km)')
    ax.set_ylabel('Count')
    
def dist_sequence(coords):
    n = len(coords)
    dist = R*metric.pairwise(np.radians(coords)).ravel()
    offsets = (np.arange(n)[:,None] - np.arange(n)).ravel()
    return pd.DataFrame([offsets[offsets>0]*dt*60,dist[offsets>0]], index=['time_offset', 'distance']).T
    
dist_ungrouped = pd.concat(df.traj.apply(dist_sequence).values).set_index('time_offset')
#print(dist_ungrouped)
dists = np.sqrt((dist_ungrouped**2).groupby('time_offset').mean()/2)

#print(dists)

def fit_rational(x,y,w=1):
    ws = np.sqrt(w)
    (a,b),_,_,_ = np.linalg.lstsq(np.column_stack([x,-y])*ws[:,None], x*y*ws, rcond=None)
    return  (a,b)

coeffs = fit_rational(dists.index.values, dists.distance.values, ((1+np.arange(len(dists)))/(1+len(dists)))**-3)

print(coeffs)

def likelihood(coords, ab):
    n = len(coords)
    a,b = coeffs
    dist = R*metric.pairwise(np.radians(coords))
    time = dt*60*np.abs(np.arange(n)[:,None] - np.arange(n))
    sigma = a*time/(time + b) + np.eye(n)
    lr = -0.5*(dist**2/sigma**2).sum(axis=1)
    return lr

def norm_lr(lr):
    return (lr-lr.max())/len(lr)

def plot_likelihood(ax, coords):
    lr = likelihood(coords,coeffs)
    lr = norm_lr(lr)
    sns.lineplot(x=np.arange(len(coords)),y=lr, ax=ax);
    ax.set_xlabel('Sequence ID')
    ax.set_ylabel('Normalized Likelihood')
import itertools   
    
n = int(df.traj.str.len().quantile(0.9))
thresh = -2
#invalid = df.traj.apply(lambda t: (norm_lr(likelihood(t,coeffs)) < thresh)[:n].tolist() + [False]*(n-len(t))).values.tolist()

final = load_data('./train.csv',nrows = -1)


def parallelize_dataframe(df, func, n_cores=4):
    df_split = np.array_split(df, n_cores)
    pool = Pool(n_cores)
    br = pool.map(func, df_split)
    out = np.concatenate(br).ravel()
    pool.close()
    pool.join()
    return out
    
def find_bad_routes(df):
    br = df.traj.apply(lambda t: (norm_lr(likelihood(t,coeffs)) < thresh).any()).values
    return br

#print("Routes with invalid points: {} / 100".format(bad_routes.sum()))

bad_routes = parallelize_dataframe(final, find_bad_routes, 70)

#print(bad_routes)
    
final[~bad_routes].to_csv('./goodData.csv', mode = 'w', columns=['POLYLINE'], index = False)

#df[~bad_routes].lines.plot(figsize=[15,15]);
#plt.show()
