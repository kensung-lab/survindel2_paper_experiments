OUTDIR=`dirname $1`
OUTDIR=`dirname $OUTDIR`
spoa $1 > $OUTDIR/consensus_fa/`basename $1`
