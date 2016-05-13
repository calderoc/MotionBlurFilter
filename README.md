#Python Modules and Jupyter (IPython) Notebooks Reproducing Illustrative Examples in:
##Motion Blur Filtering: A statistical approach for extracting confinement forces and diffusivity from a single blurred trajectory. 
===================================================

Set of Python 2.7 modules and notebooks illustrating steps for using code.  Examples reproduce Figs. 3 and 4 in article:


```bibtex
@article{PhysRevE.93.053303,
  title = {Motion blur filtering: A statistical approach for extracting confinement forces and diffusivity from a single blurred trajectory},
  author = {Calderon, Christopher P.},
  journal = {Phys. Rev. E},
  volume = {93},
  issue = {5},
  pages = {053303},
  numpages = {15},
  year = {2016},
  month = {May},
  publisher = {American Physical Society},
  doi = {10.1103/PhysRevE.93.053303},
  url = {http://link.aps.org/doi/10.1103/PhysRevE.93.053303}
}

```

```
Subdirectory contents:
* src --- folder containing python modules, notebooks and test data.  

For users who do not wish to install Jupyter or IPython Notebook, 
we have provided HTML versions of the notebooks to illustrate use of methods 
and functions provided in the Python modules.
```






If you utilize these scripts in your work, please cite the Phys. Rev. E article cited above.  


Requirements
------------

* [**IPython**](https://ipython.org/install.html)
or [**Jupyter**](http://jupyter.readthedocs.org/en/latest/install.html#new-to-python-and-jupyter)
* [**Numpy**](http://www.numpy.org/)
* [**Scipy**](http://www.scipy.org/)
* [**Matplotlib**](http://matplotlib.org/)

Routines tested on following  operating systems: Mac OS X 10.10 and Ubuntu 14.04.4 LTS. Tested using package versions: IPython 3.2, Jupyter 4.1, Numpy 1.10, Scipy 0.15.1, Matplotlib 1.4.3.


Running the Notebooks
-------

Three notebooks are provided.  Two of the notebooks provided reproduce Figures 3 and 4 in the 2016 PRE article. The notebook `Example2_ReadInFileIllustration.ipynb` provides a simple illustration of how to read trajectories from a simple plain text file.

For users with IPython v3.*, the terminal command below illustrates how to launch the notebook after dependencies above are installed. 
```bash
cd [LOCAL_GIT_REPO_LOCATION]/src
ipython notebook [FILL_IN_NB_NAME_HERE]
```
e.g., if LOCAL_GIT_REPO_LOCATION is in the current working directory and one wants to 
execute the notebook `Example2_ReadInFileIllustration.ipynb`, the following terminal commands would be used:
```bash
cd ./src
ipython notebook Example2_ReadInFileIllustration.ipynb
```

For users with Jupyter (>= v4.1), replace `ipython` by `jupyter`, e.g. 
```bash
cd ./src
jupyter notebook Example2_ReadInFileIllustration.ipynb
```





