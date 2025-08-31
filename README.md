# The Mass Gap in SU(3) Yang-Mills Theory: A Constructive Proof of the bounds and constants

**Author:** Keefe D. Reeves  


REPRODUCIBILITY / BUILD (updated)



Reproduce the test of the bounds (CLI):


```
python -m pip install -e .

pytest -q

```

## ym-bounds report command
```
pytest -q

–beta 6.0 –sym-N 8

–eta0 0.05 –A 3.0 –C 0.2 –steps 20

–tau0 0.4 –sigma-lat 0.045 –a 0.08

–area 1.0 –perimeter 4.0

–perim-scale 1.0 –perim-density 1.0

–mstar 0.3 –pref 1.0 –mp-dps 120

–csv bounds.csv –tex bounds.tex –json bounds.json
```


bounds.json/csv/tex capture inputs, outputs, precision, and git commit.
