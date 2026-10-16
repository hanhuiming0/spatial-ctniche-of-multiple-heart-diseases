import os
PROJECT_ROOT = Path(os.environ.get('PROJECT_ROOT', '.'))

"""Editable vector boundaries of threshold-defined infarct/fibrotic regions."""
from pathlib import Path
import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from scipy.spatial import cKDTree

ROOT=Path(str(PROJECT_ROOT / Path('')))
HERE=ROOT/"3-结果与分析"/"3-空间转录组"/"4-细胞轨迹"
SCORE=ROOT/"3-结果与分析"/"3-空间转录组"/"3-hotspot"/"结果文件"/"基因集打分"
H5AD=ROOT/"2-数据"/"h5ad"
OUT=HERE/"阈值分组外围曲线_3样本"
CFG=[
 dict(sample="ACH006",disease="ICM-FZ",group="ICM",score="Fibrotic Score",region="Fibrotic core",threshold=1.8,fill="#C65D57",line="#7F1D2D"),
 dict(sample="ACH0011",disease="AMI-IZ",group="AMI",score="Ischemic Score",region="Infarct core",threshold=1.0,fill="#D96B62",line="#8B1E2D"),
 dict(sample="HCM1225D",disease="HCM",group="HCM",score="Fibrotic Score",region="Fibrotic core",threshold=1.0,fill="#C65D57",line="#7F1D2D")]
mpl.rcParams.update({"pdf.fonttype":42,"ps.fonttype":42,"font.family":"Arial","savefig.transparent":True})

def load(c):
 s=c["sample"]; z=pd.read_csv(SCORE/c["group"]/f"{s}_spatial.Fibrotic_info.txt",sep="\t",index_col=0)
 a=ad.read_h5ad(H5AD/f"{s}-processed_adata_vis.h5ad",backed="r"); obs=pd.Index(a.obs_names.astype(str)); ids=[obs,pd.Index([f"{s}_{x}" for x in obs]),pd.Index([x.split("_",1)[-1] for x in obs])]
 names=next((x for x in ids if len(x.intersection(z.index.astype(str)))),None)
 if names is None: raise ValueError(f"{s}: unmatched spot IDs")
 d=pd.DataFrame(np.asarray(a.obsm["spatial"])[:,:2],index=names,columns=["x","y"])
 if a.file is not None: a.file.close()
 d=d.join(z[[c["score"]]],how="inner").rename(columns={c["score"]:"score"}).dropna(); d["positive"]=d.score>c["threshold"]
 return d

def triangulate(d):
 xy=d[["x","y"]].to_numpy(); nn=cKDTree(xy).query(xy,k=2)[0][:,1]; spacing=float(np.median(nn[nn>0])); t=mtri.Triangulation(d.x,d.y); p=xy[t.triangles]
 longest=np.maximum.reduce([np.linalg.norm(p[:,0]-p[:,1],axis=1),np.linalg.norm(p[:,1]-p[:,2],axis=1),np.linalg.norm(p[:,2]-p[:,0],axis=1)])
 t.set_mask(longest>2.25*spacing); return t,spacing

def panel(ax,d,c,legend=True):
 n,p=d.loc[~d.positive],d.loc[d.positive]; ax.scatter(n.x,n.y,s=5,c="#D9D9D9",edgecolors="none",rasterized=False,zorder=1); ax.scatter(p.x,p.y,s=6.2,c=c["fill"],edgecolors="none",rasterized=False,zorder=2)
 t,spacing=triangulate(d); binary=d.positive.astype(float).to_numpy()
 if binary.min()<.5<binary.max(): ax.tricontour(t,binary,levels=[.5],colors=[c["line"]],linewidths=[1.15],zorder=3)
 k,N=int(d.positive.sum()),len(d); ax.set_title(f"{c['sample']}  |  {c['disease']}\n{c['score']} > {c['threshold']:g}  ({k}/{N}, {100*k/N:.1f}%)",fontsize=9.2,pad=7)
 ax.set_aspect("equal"); ax.invert_yaxis(); ax.set_axis_off(); m=1.2*spacing; ax.set_xlim(d.x.min()-m,d.x.max()+m); ax.set_ylim(d.y.max()+m,d.y.min()-m)
 if legend:
  h=[Line2D([0],[0],marker="o",ls="none",ms=4.5,mfc=c["fill"],mec="none",label=f"{c['region']} (positive)"),Line2D([0],[0],marker="o",ls="none",ms=4.5,mfc="#D9D9D9",mec="none",label="Threshold-negative"),Line2D([0],[0],color=c["line"],lw=1.15,label="Threshold boundary")]
  ax.legend(handles=h,loc="lower center",bbox_to_anchor=(.5,-.035),frameon=False,fontsize=7.5,handlelength=1.5)

def main():
 OUT.mkdir(parents=True,exist_ok=True); loaded=[]; rows=[]
 for c in CFG:
  d=load(c); loaded.append((d,c)); fig,ax=plt.subplots(figsize=(4.2,4.5)); panel(ax,d,c); fig.tight_layout(pad=.5); stem=f"{c['sample']}_{c['region'].replace(' ','_')}_threshold_boundary"; fig.savefig(OUT/f"{stem}.pdf",format="pdf",bbox_inches="tight"); fig.savefig(OUT/f"{stem}.svg",format="svg",bbox_inches="tight"); plt.close(fig)
  rows.append(dict(sample=c["sample"],disease=c["disease"],score=c["score"],threshold=c["threshold"],n_spots=len(d),n_threshold_positive=int(d.positive.sum()),positive_percentage=100*d.positive.mean()))
 fig,axes=plt.subplots(1,3,figsize=(11.4,4)); [panel(ax,d,c,False) for ax,(d,c) in zip(axes,loaded)]; h=[Line2D([0],[0],marker="o",ls="none",ms=5,mfc="#C65D57",mec="none",label="Threshold-positive"),Line2D([0],[0],marker="o",ls="none",ms=5,mfc="#D9D9D9",mec="none",label="Threshold-negative"),Line2D([0],[0],color="#7F1D2D",lw=1.15,label="Threshold boundary")]; fig.legend(handles=h,loc="lower center",ncol=3,frameon=False,fontsize=8,bbox_to_anchor=(.5,-.01)); fig.subplots_adjust(left=.01,right=.99,top=.88,bottom=.13,wspace=.08); stem="ACH006_ACH0011_HCM1225D_threshold_boundaries_combined"; fig.savefig(OUT/f"{stem}.pdf",format="pdf",bbox_inches="tight"); fig.savefig(OUT/f"{stem}.svg",format="svg",bbox_inches="tight"); plt.close(fig)
 pd.DataFrame(rows).to_csv(OUT/"threshold_boundary_summary.csv",index=False,encoding="utf-8-sig"); print(OUT)

if __name__=="__main__": main()