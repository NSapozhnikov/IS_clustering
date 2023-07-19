import os

fp = '/mnt/wd/nsap/'
for i in range(1,23):
    os.system(f".././plink --bfile {fp}imp1/merged --chr {i} --make-founders --r2 square --out {fp}in_data1/chr{i}_matrix")
    os.system(f".././plink --bfile {fp}imp1/merged --chr {i} --write-snplist --out {fp}in_data1/chr{i}")
    print(f"{i} element has been plinked...")
os.system(f"rm {fp}in_data1/*.nosex {fp}in_data1/*.log")
print('Done')
