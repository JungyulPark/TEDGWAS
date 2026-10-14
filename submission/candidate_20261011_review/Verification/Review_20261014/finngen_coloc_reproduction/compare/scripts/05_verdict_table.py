import pandas as pd, numpy as np
d=pd.read_csv('results/coloc_deltas_vs_codex.csv')
t=pd.read_csv('results/TSHR_snpset_decomposition.csv'); t=t[t.snpset=='EUR_ldref_subset']
m=pd.read_csv('results/mr_recomputed_vs_codex.csv')
rows=[]
for g in ['TSHR','IGF1R','CTLA4']:
  for ep in ['FinnGen_GO','E4_GRAVES_STRICT','E4_GRAVES_OPHT_STRICT']:
    for p in ['1e-05','5e-06','1e-06']:
      a=d[(d.variant=='A_composite')&(d.gene==g)&(d.endpoint==ep)&np.isclose(d.p12,float(p))].iloc[0]
      b=d[(d.variant=='B_nfe_primary')&(d.gene==g)&(d.endpoint==ep)&np.isclose(d.p12,float(p))].iloc[0]
      r=dict(kind='coloc',gene=g,endpoint=ep,p12=p,codex_H2=a.H2_codex,codex_H3=a.H3_codex,codex_H4=a.H4_codex,A_H2=a.H2,A_H3=a.H3,A_H4=a.H4,B_H2=b.H2,B_H3=b.H3,B_H4=b.H4,
             A_dH4=a.dH4,B_dH4=b.dH4,A_max_dPP=max(abs(a[c]) for c in ['dH0','dH1','dH2','dH3','dH4']),B_max_dPP=max(abs(b[c]) for c in ['dH0','dH1','dH2','dH3','dH4']))
      if g=='TSHR':
        x=t[(t.endpoint==ep)&np.isclose(t.p12,float(p))].iloc[0]; r.update(TSHR_1000Gset_H4=x.H4,TSHR_1000Gset_dH4=x.dH4_vs_codex,TSHR_1000Gset_dH3=x.dH3_vs_codex)
      rows.append(r)
    for sc in ['sept_1000G_scale','gnomAD_composite','gnomAD_NFE']:
      x=m[(m.gene==g)&(m.endpoint==ep)&(m.scale==sc)].iloc[0]
      rows.append(dict(kind='MR_'+sc,gene=g,endpoint=ep,k=x.k,OR=x.OR,CI_lo=x.lo_q,CI_hi=x.hi_q,P=x.p,codex_OR=x.codex_OR,codex_lo=x.codex_lo,codex_hi=x.codex_hi,codex_P=x.codex_p,d_beta=x.d_beta,rel_d_P=x.rel_d_p))
pd.DataFrame(rows).to_csv('results/final_verdict_table.csv',index=False); print(len(rows))
