import os



import papermill as pm

filter = ['u','B','g','V','r','i','Y','J','H']

for val in filter:
    pm.execute_notebook(
        "Mark_H0_tripp_alpha.ipynb",
        f"output_notebooks/output_{val}.ipynb",
        parameters={"input_arg": val}
    )



os.system('say "your program has finished."')


sys.exit()
#cal = ['ceph','trgb','sbf','sbfj21','sbfcombined']
cal = ['sbfj21']

for cal in cal:
    filter = ['u','B','g','V','r','i','Y','J','H']
    #filter = ['B','H']
    for filter in filter:        
        file = filter+'_'+cal+'_update3.csv'
        os.system("python Mark_H0_tripp_alpha.ipynb "+file)

os.system('say "your program has finished."')

