#! /usr/bin/env python
# encoding: utf-8
"""
do_night_malt.py fully processes a night of data.
Should be run at the end of a night.
By default it reduces the current UT date and
UT - 1 day (since some observing sessions run
over midnight UT).

Options
-n : Night -- reduce all files on a give date (YYYY-MM-DD)
-t : Today -- only reduce files from UT today or night as
              specified above.


"""
import sys
import preprocess_malt, reduce_malt, ReduceLog
import datetime
import getopt

def main():
	try:
		opts,args = getopt.getopt(sys.argv[1:], "n:t")
	except getopt.GetoptError,err:
		print(str(err))
		print(__doc__)
		sys.exit(2)
	
	do_yesterday = True
	now_utc = datetime.datetime.utcnow()
        night   = datetime.date(now_utc.year,now_utc.
			       month,now_utc.day).isoformat()

	for o,a in opts:
		if o == "-n":
			night = a
		if o == "t":
			do_yesterday = False

	
	dlist = night.split('-')
	today = datetime.date(int(dlist[0]),int(dlist[1]),int(dlist[2]))
	yesterday = today - datetime.timedelta(days=1)

	daylist = [today]
	if do_yesterday:
		daylist.append(yesterday)
	
	print("### Running do_night_malt.py for the following UT dates: ##")
	for date in daylist:
		print(date.isoformat())

	for date in daylist:
		
		files_to_process = preprocess_malt.get_new_files(date.isoformat(),
							in_middle_of_obs = False)
		preprocess_malt.rename_files(files_to_process)

		redlog = ReduceLog.ReduceLog()
		sources = redlog.find_files_with_date(date.isoformat())
#		print(sources)
		for one_source in sources:
			reduce_malt.do_reduction(one_source)

if __name__ == '__main__':
	main()
