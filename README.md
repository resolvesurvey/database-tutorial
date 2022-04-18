## RESOLVE survey database tutorial

To complete this tutorial:

 * go to https://jupyter.org/try
 * click "Try JupyterLab"
 * close open tabs in the Lab (not necessary, just less confusing)
 * obtain the tutorial: `File` > `Open from URL` > copy+paste `https://raw.githubusercontent.com/resolvesurvey/database-tutorial/master/RESOLVE_DataTutorial.ipynb`.
 * to obtain the RESOLVE data file and the Eckert+ 2015 code, click `File` > `New Console for Notebook`. In the pop-up console, copy+paste and run the following code:
 ```
from js import fetch
url1 = 'https://raw.githubusercontent.com/resolvesurvey/database-tutorial/master/RESOLVE_DR4_prepubl.csv'
contents = await fetch(url1)
contents = await contents.text()
with open('RESOLVE_DR4_prepubl.csv','w') as ff:
    ff.write(contents)
ff.close()

url2 = 'https://raw.githubusercontent.com/keckert7/codes/master/pred_loggs_dist.py'
contents = await fetch(url2)
contents = await contents.text()
with open('pred_loggs_dist.py','w') as ff:
    ff.write(contents)
ff.close()
import scipy
 ```
* close the pop-up console to return to the tutorial.
* you can run or re-run individual cells by clicking on them and typing Ctrl-Enter
