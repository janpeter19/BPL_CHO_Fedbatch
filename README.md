# BPL_CHO_Fedbatch

This example of cultivation of CHO culture using fedbatch technique is in laboratory scale. The model describes by-products (lactate and ammonia) formation at over-feeding based on a publication. The original model is extended to describe also recombinant protein production. Simulation is done using an FMU from Bioprocess Library *for* Modelica. Below a diagramwith a typical simulation that you will get at the end of the Jupyter notebook.

![](Fig4_BPL_CHO_Fedbatch.png)

You see in the diagram typical aspects of a CHO Fedbatch process for a certain type of recombinant protein production. 

* The solid line shows the results of constant substrate feeding. The dashed line shows the improved result on recombinant protein production by decreasing the feed rate somewhat at time 100 hours. This behaviour is typical for CHO processes where the recombinant protein production is negatively affected by cell growth.
 
* In the diagrams to the right for specific glucose and glutamine uptake, the red line shows the bottleneck of metabolism. 

You start up the notebook in Colab by pressing here
[start BPL notebook](https://colab.research.google.com/github/janpeter19/BPL_CHO_Fedbatch/blob/main/BPL_CHO_Fedbatch_colab_me.ipynb)

Then you in the menu choose Runtime/Run all. The installation takes just a few minutes. The subsequent execution of the simulations of microbial growth take just a second or so. 

You can continue in the notebook and make new simulations and follow the examples given. Here are many things to explore!

Note that:
* The script occassionaly get stuck during installation. Then just close the notebook and start from scratch.
* Runtime warnings are at the moment silenced. The main reason is that we run with an older combination of PyFMI and Python that bring depracation warnings of little interest. 
* Remember, you need to have a google-account!

Just to be clear, no installation is done at your local computer.

Work in progress - stay tuned!

License information:
* The binary-file with extension FMU is shared under the permissive MIT-license
* The other files are shared under the GPL 3.0 license
