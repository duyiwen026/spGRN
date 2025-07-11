#########                  
data = st.Read10X(data_dir)
data.var_names_make_unique()
st.add.image(adata=data,
             imgpath=data_dir+"spatial/tissue_hires_nobg.png",
             library_id="V1_Breast_Cancer_Block_A_Section_1", visium=True)

# Basic normalisation #
st.pp.filter_genes(data, min_cells=3)
st.pp.normalize_total(data) # NOTE: no log1p                  
                  
# Adding the label transfer results,  #
spot_mixtures = pd.read_csv('/data/stlearn/endspot25.csv', index_col=0)
labels = spot_mixtures.loc[:,'predicted.id'].values.astype(str)
spot_mixtures = spot_mixtures.drop(['predicted.id','prediction.score.max'],
                                   axis=1)
spot_mixtures.columns = [col.replace('prediction.score.', '')
                         for col in spot_mixtures.columns]

# Note the format! #
print(labels)
print(spot_mixtures)
print('Spot mixture order correct?: ',
      np.all(spot_mixtures.index.values==data.obs_names.values)) # Check is in correct order

# NOTE: using the same key in data.obs & data.uns
data.obs['cell_type'] = labels # Adding the dominant cell type labels per spot
data.obs['cell_type'] = data.obs['cell_type'].astype('category')
data.uns['cell_type'] = spot_mixtures # Adding the cell type scores

st.pl.cluster_plot(data, use_label='cell_type')                  
                  
# Loading the LR databases available within stlearn (from NATMI)
lrs = st.tl.cci.load_lrs(['connectomeDB2020_lit'], species='human')
print(len(lrs))

# Running the analysis #
st.tl.cci.run(data, lrs,
                  min_spots = 20, #Filter out any LR pairs with no scores for less than min_spots
                  distance=None, # None defaults to spot+immediate neighbours; dis-tance=0 for within-spot mode
                  n_pairs=200, # Number of random pairs to generate; low as example, recommend ~10,000
                  n_cpus=14, # Number of CPUs for parallel. If None, detects & use all available.
                  )

