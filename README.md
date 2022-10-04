# BPL_CHO_Fedbatch

This example of cultivation of CHO culture using fedbatch technique is in laboratory scale. The model describes by-products (lactate and ammonia) formation at over-feeding based on a publication. The original model is extended to describe also recombinant protein production. Simulation is done using an FMU from Bioprocess Library *for* Modelica. Below a diagramwith a typical simulation that you will get at the end of the Jupyter notebook.

![](Fig4_BPL_CHO_Fedbatch.png)

You see in the diagram typical aspects of a CHO Fedbatch process for a certain type of recombinant protein production. The solid line shows the results of constant substrate feeding. The dashed line shows the improved result on recombinant protein production by decreasing the feed rate somewhat at time 100 hours. This behaviour is typical for CHO processes where the recombinant protein production is negatively affected by cell growth.
 
 [start BPL notebook](https://colab.research.google.com/github/janpeter19/BPL_CHO_Fedbatch/blob/main/BPL_CHO_Fedbatch_colab.ipynb)
 
Work in progress - stay tuned!

License information:
* The binary-file with extension FMU is shared under the permissive MIT-license
* The other files are shared under the GPL 3.0 license
