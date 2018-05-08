## Fine mapping benchmark

This benchmark is written in DSC with all DSC configurations and module scripts
in SoS Notebook. To edit the notebook please use Jupyter Lab:

```
pip install sos-notebook jupyterlab
conda install nodejs
jupyter labextension install jupyterlab-toc
jupyter labextension install jupyterlab-sos
```

Then open the notebook with:

```
jupyter-lab finemapping.ipynb
```

Look for the `Contents` vertical tab in JupyterLab's GUI to bring up the table
of contents for navigation.

To run the benchmark, after editing the notebook first save them to module files

```
./run.sos convert
```

Then run the benchmark from `mnm.sh`. To see what are available for example,

```
dsc mnm.dsc -h
```

**Do not edit `mnm.dsc` or any scripts under `modules/` folder. Always use Jupyter Lab to edit and export.**
