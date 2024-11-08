## Adding a model that 'ships' with this package

Step 1. Create a folder here for the model(s). 

Step 2. Add a csv file for the/each new model in the new folder,
each file should have two columns, 'pred.var' and 'coef', e.g.

```
pred.var,coef
"pred.var","coef"
"(Intercept)",0.695507258
"cg00075967",0.12933661
"cg00374717",0.005017857
"cg00864867",1.59976405
"cg00945507",0.056852418
```

Step 3. Ensure that the file does not have a 'Byte Order Mark'
as these will cause problems in some versions of R. 
This can be achieved at the commmand line as follows:

```
vi -c ":set nobomb" -c ":wq" new-model.csv
```

Step 4. Add each new model to the 'models.csv' file in this folder, e.g.

```
name,tissue,target,publication,source,filename
brenner,blood,mortality,10.1038/ncomms14617,Supplementary Figure 1,brenner/coefs.csv
```





