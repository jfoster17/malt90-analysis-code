#!/bin/sh
#
# Copy MOPS spectrometer data for the MALT project.
# This is supposed to run on host 'draco'.
# Example crontab entry - 
#  # run every 15min, starting at 3min past the hour.
#  03,18,33,48 * * * * /bin/sh /path/to/this/file
#
# $URL: https://svn.atnf.csiro.au/operations/TOSS/packages/atca-operations/trunk/data-handling/rsync_data.sh $
# $Id: rsync_data.sh 2567 2010-08-16 04:59:04Z VincentMcIntyre $
# ----------------------------------------------------------------------------
# Config
DEBUG=0
SHOW_LOG=0   # set =1 to display log at end. see also -v option.

DESTINATION="/DATA/MALT_1/MALT90/raw_data/"
LOGFILE="${DESTINATION}/rsync_malt90.log"

SOURCE="kaputar::MALT90"

# ensure we get the expected file permissions at $DESTINATION
umask 002

RSYNC_OPTIONS="--archive --update --no-perms --append --exclude archived --exclude exported"
# --archive  : archive mode; equals -rlptgoD (no -H,-A,-X).
# --update   : skip files that are newer on the receiver.
# --no-perms : don't copy the permssions, use the process' umask.
#              Read the rsync manpage to fully appreciate this option.
# --append   : append data onto incomplete files, don't create new files.
#              This is needed to handle files that are only partly written
#              at the time this job runs.
# --exclude  : on some data disks there are some subdirectories that are
#              used to track which data files had been written to cdrom
#              or tape. These don't need to be archived.
#            

# ----------------------------------------------------------------------------
# Functions
Fail () {
    SHOW_LOG=1
    while test $# -gt 0 ; do echo "$1"; shift; done; exit 1;
} # Fail

CleanUp() {
    test -f "$PIDFILE" && rm -f "$PIDFILE"
    if test -f "$TMPLOG" ; then
	# undo log redirection
	echo "#Exiting at `date +%Y-%m-%dT%H-%M-%S%z` ---------------------------"
        exec 1>&6 6>&- ; exec 2>&7 7>&-

	if test -f "$LOGFILE"; then
            cat "$TMPLOG" >> "$LOGFILE"
	    test "X1" = "X$DEBUG"  && echo "Updated '$LOGFILE'"
	fi

	test "X1" = "X$SHOW_LOG"   && cat "$TMPLOG"
	rm -f "$TMPLOG"
    fi
    
} # CleanUp

# ----------------------------------------------------------------------------
PATH=/usr/bin:/bin; export PATH; unalias -a

# force cleanup on exit
trap "CleanUp" 0

# process commandline arguments
while test $# -gt 0
do
    case "$1" in
	-[hH]) #help
               echo "Usage: $0 [-h*elp] [-v*erbose] [-t*est]"
	       echo "You need to edit the script to set SOURCE and DESTINATION"
	       exit 0
        ;;
	-[tT]) #test mode - no data will be transferred.
               DEBUG=1
        ;;
	-[vV]) #verbose
               SHOW_LOG=1
        ;;
    esac
    shift
done

test 1 -eq $DEBUG && RSYNC_OPTIONS="--dry-run --verbose $RSYNC_OPTIONS"
test 1 -eq $DEBUG && SHOW_LOG="1"

# The log for this invocation doesn't need a fixed location
TMPLOG=`mktemp /tmp/XXXXXX`

# store all output in tmplog, decide at the end if we need to show it to the user
exec 6>&1 ; exec 7>&2; exec 1> "$TMPLOG" 2>&1
echo "#Starting at `date +%Y-%m-%dT%H-%M-%S%z` --------------------------"


# PIDFILE needs to be a predictable path, we can't naively use mktemp.
prog=`basename $0`
PIDFILE="/tmp/${prog}.pid"

# Ditto the long-term log.
test -f "$LOGFILE" || touch -m "$LOGFILE" || Fail "Can't create '$LOGFILE'"

# are we running already?
if test -f "$PIDFILE"
then
    pid=`cat "$PIDFILE"`
    # NB: behaviour of 'ps' is quite variable between platforms
    x=`(ps ax | awk '$1 == "'${pid}'"{print $0}' |grep "$prog") 2>/dev/null`

    # That process is still running, so don't start another.
    if ! test "X" = "X${x}"
    then
	echo "Already running as process $pid, exiting."
        # unset exit trap so we don't kill pidfile of running job
	trap - 0
        exit
    fi
fi
echo $$ >"$PIDFILE" || Fail "Unable to update pidfile '$PIDFILE'"


# check we have what we need to run
# NFS automounting hack
for d in "$SOURCE" "$DESTINATION"
do
    case "X${d}" in
        X/*) ls -d "$d" 1>/dev/null 2>&1
        ;;
    esac
done
sleep 1

test -d "$DESTINATION"  || Fail "$0: no destination directory '$DESTINATION'"

# get on with it
case "X$SOURCE" in
    X/*) # NFS or local filesystem path
	test -d "$SOURCE"       || Fail "$0: no source directory '$SOURCE'"
	nice rsync $RSYNC_OPTIONS "${SOURCE}" "${DESTINATION}" || \
            Fail "$0: failure during rsync"
    ;;

    X*@*::*) # Authenticated rsync - not handled yet.
	Fail "Authenticated rsync is not a supported source ... yet."
    ;;

    X*::*) # Anonymous rsync
	nice rsync $RSYNC_OPTIONS "${SOURCE}" "${DESTINATION}" || \
            Fail "$0: failure during rsync"
    ;;

    X*@*:/*) # Remote login & specific path on remote end. We need ssh
	test -f "$SSHKEY"       || Fail "No ssh key ($SSHKEY) available, aborting"
	nice rsync $RSYNC_OPTIONS -e "$SSHCMD" "${SOURCE}" "${DESTINATION}" || \
            Fail "$0: failure during rsync"
    ;;
esac


exit 0
