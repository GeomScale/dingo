This is an example directory where you need to add all your `.mat` or `.json` model files 
your community consists of. 

For example, you may have a community of 2 models, *Escherichia coli* and *Helicobacter pylori*, 
in that case you could get the [`e_coli_core.json`](http://bigg.ucsd.edu/models/e_coli_core) and 
[`iIT341.json`](http://bigg.ucsd.edu/models/iIT341) files from the [BIGG database](http://bigg.ucsd.edu/models)
and then run `dingo`. 

Remember to remove this `README.md` file if you will decide to use this directory. 

Here is an example of how to run `dingo` from your terminal for the case of a community analysis:

```bash=
python -m dingo -cmd  </path_to_community_models_dir/>  --format json
```



