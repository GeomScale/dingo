## using only `dingo`

Make some folders:

```
mkdir polytopes_redundant_removed
mkdir dingo_samples_redundant_removed
mkdir polytopes_redundant_in
mkdir dingo_samples_redundant_in
```

To sample only using `dingo` features (withn and without removing redundant facets), run:

```
./dingo_all_the_way.py  iEC1344_C.json 
```


## `PolyRound` to get a P and then sample with `dingo`

To use the PolyRound transformations: 

again, some directories: 

```
mkdir simpl_transf_polytopes
mkdir polyrounded_polytopes
mkdir dingo_samples_on_polyrounded_polytopes
```

and then first, run the tranformations: 

```
./polyround_preproces.py iEC1344_C.xml
```

> Remember! 
> * For the `PolyRound` transformations, we use the `.xml` format of the model.
> * For our experiments, `dingo` will run only on simplified + transformed polytopes from the `PolyRound` algorithm (that will be saved under the `simpl_transf_polytopes` directory) while, `hopsy` will use simplified + transformed + rounded polytopes from `PolyRound` (will be built under the `polyrounded_polytopes` folder)


and then sample with `dingo`

```
./dingo_on_polyrounded_polyt.py polytope_iEC1344_C.xml.pckl 1000
```

where `1000` stands for the value of ESS required. 


For the case of Recon2, we used [version 2.2](https://www.ebi.ac.uk/biomodels/MODEL1603150001#Files).




## `PolyRound` to get a P and then sample with `hopsy`

To perfmorm sampling using the simplified + transformed + rounded polytopes from PolyRound and `hopsy`:

```
./hopsy_on_polyrounded_polytopes.py polyrounded_polytopes/polytope_iEC1344_C.xml.pckl
```
