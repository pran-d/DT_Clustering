# DTAg Clustering
The proposed algorithm involves a two-step process: geometric Delaunay-based filtering and MSC-based merging. 
- The first step utilizes the inherent properties of Delaunay triangles to distinguish between intra-cluster and inter-cluster triangles, thus identifying preliminary clusters.
- The second step refines the initial clustering by merging the primitive clusters based on their similarity, which is defined using the features of their Minimum Spanning Circles.

DTAg's hyperparameters are tunable in both steps, providing flexibility in achieving desired clustering outcomes. The provision of two convergence criteria allows for intuitive and unambiguous tuning as well as higher user control over the result.
The two versions of the algorithm have been provided in `DTAg_Cascading_Merges.py` and `DTAg_k_merges.py`. 

The evaluation and scoring code is included in `genie_clust.py`. Note that the algorithm code saves the results to a .csv file and is compared to the ground truth data which must be provided in the same directory as another .csv file.
