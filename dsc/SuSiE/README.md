## Univariate fine-mapping benchmark

This benchmark is written in DSC with all DSC configurations and module scripts
in SoS Notebook. To edit the notebook please use Jupyter Lab:

```
pip install sos-notebook jupyterlab
conda install nodejs
jupyter labextension install --no-build jupyterlab-toc jupyterlab-sos
```

Then open the notebook with:

```
jupyter-lab finemapping.ipynb
```

Look for the `Contents` vertical tab in JupyterLab's GUI to bring up the table
of contents for navigation.

To run the benchmark, after editing the notebook first save them to module files

```
./export.sos
```

Then run the benchmark from `*.dsc`. To see what are available for example,

```
dsc susie.dsc -h
```

**Do not edit `*.dsc` under this folder, or any scripts under `modules/` folder. 
Always use JupyterLab to edit and export.**
