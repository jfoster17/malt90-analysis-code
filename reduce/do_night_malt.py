#!/usr/bin/env python
# encoding: utf-8
"""
do_night_malt.py fully processes a night of data.
Should be run at the end of a night
after rsyncing data across (if required).
"""
import sys
import preprocess_malt, reduce_malt, ReduceLog
from datetime import date, timedelta

night = sys.argv[1]

dlist = night.split('-')
today = date(int(dlist[0]),int(dlist[1]),int(dlist[2]))
yesterday = today - timedelta(days=1)

print(today.isoformat())
print(yesterday.isoformat())

for date in [today]:

	files_to_process = preprocess_malt.get_new_files(date.isoformat(),in_middle_of_obs = False)
	preprocess_malt.rename_files(files_to_process)

	redlog = ReduceLog.ReduceLog()
	sources = redlog.find_files_with_date(date.isoformat())
	print(sources)
	for one_source in sources:
		reduce_malt.do_reduction(one_source)
