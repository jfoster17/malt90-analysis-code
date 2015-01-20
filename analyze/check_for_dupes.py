import os

year1 = os.listdir("/DATA/MALT_1/MALT90/data/year1/sources")
year2 = os.listdir("/DATA/MALT_1/MALT90/data/year2/sources")
year3 = os.listdir("/DATA/MALT_1/MALT90/data/year3/sources")
year4 = os.listdir("/DATA/MALT_1/MALT90/data/year4/sources")

for source in year4:
    print(source)
    if (source in year1) or (source in year2) or (source in year3):
        print("Dup")
