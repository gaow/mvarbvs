## Fine mapping benchmark

This benchmark is written in DSC with all DSC configurations and module scripts
in SoS Notebook. To edit the notebook please use Jupyter Lab:

```
pip install jupyterlab
conda install nodejs
jupyter labextension install jupyterlab-toc
```

Then clone and install an alpha version of `jupyterlab-sos`:

```
git clone https://github.com/vatlab/jupyterlab-sos
cd jupyterlab-sos
./build.sh
```

To run the benchmark, after editing the notebook first save them to module files

```
./run.sos convert
```

Then run the benchmark from `mnm.sh`. To see what are available for example,

```
dsc mnm.sh -h
```
